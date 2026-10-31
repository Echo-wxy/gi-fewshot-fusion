"""Study 2, final protocol - survival prediction with a clinical anchor, two-stage fusion and diagnostics.

Reuses the prepared data of exp2_model.py (run its `prep` step first) and the same outer splits
(outer 5-fold x seeds 0,1,2, stratified by event).

Inside every outer training fold an inner 3-fold structure is used. Each inner model trains on 2/3 of the outer
training patients (early stopping on validation NLL from a further 20 % split), predicts the remaining 1/3
(giving out-of-fold values for the training patients) and the outer test patients. Outer test predictions are
the average of the three inner models.

Discrete-time head with complementary log-log link: log cumulative hazard in time bin k = b_k + anchor + r(x),
anchor = clinical Cox linear predictor (out-of-fold for training patients) for the "Cox+" methods, 0 otherwise.

Two-stage fusion: unimodal encoders (WSI attention-MIL, omics linear, clinical MLP) are trained with their own
heads, frozen, and only the fusion layer is trained. PORPOISE-style and MCAT-style stay end-to-end.

Usage: python exp2_v2.py [--cohorts STAD,CRC] [--smoke]
"""
import argparse, copy, math, os, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exp2_model as M                                                         # noqa: E402

D, DROP, WD, LR, EPOCHS, PATIENCE, BATCH, INNER = 64, 0.5, 1e-3, 2e-4, 100, 20, 16, 3
LR_HEAD = 1e-3          # stage-2 fusion layers and anchored heads (few parameters, frozen inputs)
FUSIONS = ["Concat", "Attention", "Gated", "DCA"]
ORDER = ["Clinical-Cox", "Cox+zero-residual", "Omics-CoxNet", "Stacked-Cox", "Late-Fusion",
         "WSI-only", "Omics-only", "Concat-Fusion", "Attention-Fusion", "Gated-Fusion", "DCA",
         "PORPOISE-style", "MCAT-style",
         "Cox+WSI", "Cox+Omics", "Cox+Concat", "Cox+Attention", "Cox+Gated", "Cox+DCA"]
PRIMARY = ("Cox+DCA", "Clinical-Cox")
CONTRASTS = [("Cox+DCA", "Clinical-Cox"), ("Cox+DCA", "Cox+Attention"), ("Cox+DCA", "Cox+Concat"),
             ("Cox+DCA", "Cox+Gated"), ("Cox+DCA", "Stacked-Cox"), ("Cox+DCA", "Late-Fusion"),
             ("Cox+WSI", "Clinical-Cox"), ("Cox+Omics", "Clinical-Cox"), ("DCA", "Attention-Fusion"),
             ("DCA", "Concat-Fusion"), ("DCA", "Gated-Fusion")]


# ============================================================================ torch modules
def modules():
    import torch
    import torch.nn as nn

    class ABMIL(nn.Module):
        def __init__(self, d_in):
            super().__init__()
            self.fc = nn.Sequential(nn.Linear(d_in, D), nn.ReLU(), nn.Dropout(DROP))
            self.a, self.b, self.w = nn.Linear(D, 32), nn.Linear(D, 32), nn.Linear(32, 1)

        def forward(self, x, m):
            h = self.fc(x)
            s = self.w(torch.tanh(self.a(h)) * torch.sigmoid(self.b(h))).squeeze(-1).masked_fill(m == 0, -1e4)
            return (torch.softmax(s, 1).unsqueeze(-1) * h).sum(1)

    class Head(nn.Module):
        """cloglog discrete-time head: log cumulative hazard per bin = b_k + anchor + r(z) (r scalar)."""
        def __init__(self, d_in, n_bins, zero=False):
            super().__init__()
            self.b = nn.Parameter(torch.linspace(-3.0, -1.0, n_bins))
            self.r = None if zero else nn.Linear(d_in, 1)

        def forward(self, z, anchor):
            r = 0.0 if self.r is None else self.r(z).squeeze(-1)
            return self.b.unsqueeze(0) + (anchor + r).unsqueeze(-1)            # eta: B x K

    class Fusion(nn.Module):
        def __init__(self, kind):
            super().__init__()
            self.kind = kind
            self.norm = nn.LayerNorm(D)
            if kind == "Concat":
                self.mlp = nn.Sequential(nn.Linear(3 * D, D), nn.ReLU(), nn.Dropout(DROP))
            else:
                self.gate = nn.Linear(D, 1)
            if kind == "Gated":
                self.ctx = nn.Sequential(nn.Linear(2 * D, 32), nn.ReLU(), nn.Linear(32, 1))
            if kind == "DCA":
                self.conf = nn.Sequential(nn.Linear(2 * D, 32), nn.ReLU(), nn.Linear(32, 1))
                self.gamma_raw = nn.Parameter(torch.tensor(0.0))

        def forward(self, z3):                                                # z3: B x 3 x D (frozen embeddings)
            z = self.norm(z3)
            if self.kind == "Concat":
                return self.mlp(z.flatten(1)), {}
            q = self.gate(z).squeeze(-1)
            ex = {}
            if self.kind == "Gated":
                q = q + self.ctx(torch.cat([z, (z.sum(1, keepdim=True) - z) / 2], -1)).squeeze(-1)
            if self.kind == "DCA":
                zi, zj = z.unsqueeze(2), z.unsqueeze(1)
                pair = torch.cat([(zi - zj).expand(-1, -1, 3, -1), (zi * zj).expand(-1, -1, 3, -1)], -1)
                delta = torch.sigmoid(self.conf(pair).squeeze(-1))
                delta = 0.5 * (delta + delta.transpose(1, 2)) * (1 - torch.eye(3, device=z.device)).unsqueeze(0)
                q = q - torch.nn.functional.softplus(self.gamma_raw) * delta.sum(2)
                ex["delta"] = delta
            a = torch.softmax(q, 1)
            ex["alpha"] = a
            return (a.unsqueeze(-1) * z).sum(1), ex

    return ABMIL, Head, Fusion


def cloglog_nll(eta, k, c):
    """eta: B x K log cumulative hazard per bin; k: bin index; c: 1 if censored."""
    import torch
    lam = torch.exp(eta.clamp(-20, 15))
    cumH = torch.cumsum(lam, 1)
    logS = -cumH
    logS_prev = torch.cat([torch.zeros_like(logS[:, :1]), logS], 1).gather(1, k.unsqueeze(1)).squeeze(1)
    log_h = torch.log((-torch.expm1(-lam.gather(1, k.unsqueeze(1)).squeeze(1))).clamp_min(1e-12))
    logS_k = logS.gather(1, k.unsqueeze(1)).squeeze(1)
    return (-(1 - c) * (logS_prev + log_h) - c * logS_k).mean()


def eta_to_S(eta):
    return np.exp(-np.cumsum(np.exp(np.clip(eta, -20, 15)), 1))


# ============================================================================ generic trainer
def fit(step_fn, params, n_fit, val_fn, dev, seed, lr=None):
    """Mini-batch training with early stopping on validation NLL. step_fn(batch_idx) -> loss; val_fn() -> NLL.
    Returns diagnostics."""
    import torch
    opt = torch.optim.Adam(params, lr=lr or LR, weight_decay=WD)
    rng = np.random.default_rng(seed)
    best, best_ep, bad, gnorm, t0 = np.inf, 0, 0, [], time.time()
    state = None
    tr_loss = np.nan
    for ep in range(EPOCHS):
        order = rng.permutation(n_fit)
        losses = []
        for b in range(0, n_fit, BATCH):
            loss = step_fn(order[b:b + BATCH])
            opt.zero_grad(); loss.backward()
            gnorm.append(float(torch.nn.utils.clip_grad_norm_(params, 10.0)))
            opt.step(); losses.append(float(loss.detach()))
        tr_loss = float(np.mean(losses))
        v = val_fn()
        if v < best - 1e-4:
            best, best_ep, bad, state = v, ep, 0, copy.deepcopy([p.detach().clone() for p in params])
        else:
            bad += 1
            if bad >= PATIENCE:
                break
    if state is not None:
        with torch.no_grad():
            for p, s in zip(params, state):
                p.copy_(s)
    return {"best_epoch": best_ep + 1, "stop_epoch": ep + 1, "best_val_nll": best, "last_train_nll": tr_loss,
            "mean_grad_norm": float(np.mean(gnorm)) if gnorm else np.nan, "seconds": time.time() - t0}


# ============================================================================ one inner fold
def run_inner(data, fit_i, va_i, pred_i, anchor, dev, seed, diag, tag):
    """Train everything on fit_i (early stop on va_i); return {method: eta or S for pred_i} plus DCA extras."""
    import torch
    ABMIL, Head, Fusion = modules()
    bags, Xo, Xc, kbin, cens, nb = (data[k] for k in ("bags", "Xo", "Xc", "kbin", "cens", "n_bins"))
    k_t = torch.tensor(kbin, device=dev); c_t = torch.tensor(cens, device=dev, dtype=torch.float32)
    Xo_t = torch.tensor(Xo, device=dev); Xc_t = torch.tensor(Xc, device=dev)
    A_t = torch.tensor(anchor, device=dev, dtype=torch.float32)
    zero = torch.zeros_like(A_t)
    out, ex_out = {}, {}

    def nll_on(fn, idx, anc):
        with torch.no_grad():
            tot = 0.0
            for b in range(0, len(idx), 64):
                bi = idx[b:b + 64]
                tot += float(cloglog_nll(fn(bi), k_t[bi], c_t[bi])) * len(bi)
            return tot / len(idx)

    def predict(fn, idx):
        with torch.no_grad():
            return np.concatenate([fn(idx[b:b + 64]).cpu().numpy() for b in range(0, len(idx), 64)])

    # ---- stage 1: unimodal encoders with their own (unanchored) heads -------------------------------------
    enc = {"wsi": ABMIL(bags[0].shape[1]).to(dev),
           "omics": torch.nn.Sequential(torch.nn.Dropout(DROP), torch.nn.Linear(Xo.shape[1], D)).to(dev),
           "clin": torch.nn.Sequential(torch.nn.Linear(Xc.shape[1], 16), torch.nn.ReLU(), torch.nn.Linear(16, D)).to(dev)}

    def encode(m, bi):
        if m == "wsi":
            xw, mw = M.pad_batch(bags, bi, dev)
            return enc["wsi"](xw, mw)
        return enc[m](Xo_t[bi] if m == "omics" else Xc_t[bi])

    for m in ("wsi", "omics", "clin"):
        torch.manual_seed(seed)
        head = Head(D, nb).to(dev)
        params = list(enc[m].parameters()) + list(head.parameters())
        fn = lambda bi, m=m, head=head: head(encode(m, bi), zero[bi])
        def step(b, fn=fn): enc[m].train(); head.train(); return cloglog_nll(fn(fit_i[b]), k_t[fit_i[b]], c_t[fit_i[b]])
        def val(fn=fn): enc[m].eval(); head.eval(); return nll_on(fn, va_i, zero)
        d = fit(step, params, len(fit_i), val, dev, seed)
        diag.append({**tag, "model": f"stage1-{m}", **d})
        enc[m].eval(); head.eval()
        if m == "wsi":
            out["WSI-only"] = predict(fn, pred_i)
        elif m == "omics":
            out["Omics-only"] = predict(fn, pred_i)
    for e in enc.values():
        for p in e.parameters():
            p.requires_grad_(False)
    n_all = len(bags)
    with torch.no_grad():
        Z = torch.stack([torch.cat([encode(m, np.arange(b, min(b + 32, n_all))) for b in range(0, n_all, 32)])
                         for m in ("wsi", "omics", "clin")], 1)                          # N x 3 x D
    z_mean_fit = Z[fit_i].mean(0, keepdim=True)

    # ---- anchored unimodal heads on frozen embeddings, and Cox + zero residual -------------------------------
    for name, mi in (("Cox+WSI", 0), ("Cox+Omics", 1), ("Cox+zero-residual", None)):
        torch.manual_seed(seed)
        head = Head(D, nb, zero=mi is None).to(dev)
        fn = lambda bi, head=head, mi=mi: head(Z[bi, mi] if mi is not None else Z[bi, 0], A_t[bi])
        def step(b, fn=fn, head=head): head.train(); return cloglog_nll(fn(fit_i[b]), k_t[fit_i[b]], c_t[fit_i[b]])
        def val(fn=fn, head=head): head.eval(); return nll_on(fn, va_i, A_t)
        d = fit(step, list(head.parameters()), len(fit_i), val, dev, seed, lr=LR_HEAD)
        diag.append({**tag, "model": name, **d})
        head.eval()
        out[name] = predict(fn, pred_i)

    # ---- stage 2: fusion layers on frozen embeddings (unanchored and anchored) ---------------------------------
    for kind in FUSIONS:
        for anchored in (False, True):
            torch.manual_seed(seed)
            fus, head = Fusion(kind).to(dev), Head(D, nb).to(dev)
            anc = A_t if anchored else zero
            fn = lambda bi, fus=fus, head=head, anc=anc: head(fus(Z[bi])[0], anc[bi])
            def step(b, fn=fn, fus=fus, head=head):
                fus.train(); head.train(); return cloglog_nll(fn(fit_i[b]), k_t[fit_i[b]], c_t[fit_i[b]])
            def val(fn=fn, fus=fus, head=head, anc=anc): fus.eval(); head.eval(); return nll_on(fn, va_i, anc)
            d = fit(step, list(fus.parameters()) + list(head.parameters()), len(fit_i), val, dev, seed, lr=LR_HEAD)
            name = (f"Cox+{kind}" if anchored else {"Concat": "Concat-Fusion", "Attention": "Attention-Fusion",
                                                   "Gated": "Gated-Fusion", "DCA": "DCA"}[kind])
            diag.append({**tag, "model": name, **d})
            fus.eval(); head.eval()
            out[name] = predict(fn, pred_i)
            if kind in ("DCA", "Attention"):
                with torch.no_grad():
                    f_all, ex = fus(Z[pred_i])
                    a = ex["alpha"].cpu().numpy()
                    rec = {"alpha": a, "entropy": -(a * np.log(a + 1e-12)).sum(1)}
                    if "delta" in ex:
                        rec["delta"] = ex["delta"].cpu().numpy()
                    # modality removal: replace one modality by its training mean, change in NLL on pred_i
                    base = float(cloglog_nll(head(f_all, anc[pred_i]), k_t[pred_i], c_t[pred_i]))
                    for mi, mname in enumerate(("wsi", "omics", "clin")):
                        Zr = Z[pred_i].clone(); Zr[:, mi] = z_mean_fit[:, mi]
                        rec[f"dNLL_remove_{mname}"] = float(cloglog_nll(head(fus(Zr)[0], anc[pred_i]), k_t[pred_i], c_t[pred_i])) - base
                    rec["norms"] = Z[pred_i].norm(dim=-1).mean(0).cpu().numpy()
                    ex_out[name] = rec

    # ---- end-to-end PORPOISE-style and MCAT-style (sigmoid-hazard heads as in exp2_model) ----------------------
    sizes = data["group_sizes"]
    groups = [list(range(sizes[i], sizes[i + 1])) for i in range(len(sizes) - 1)]
    for name in ("PORPOISE-style", "MCAT-style"):
        torch.manual_seed(seed)
        net = M.build_model(name, bags[0].shape[1], Xo.shape[1], Xc.shape[1], nb, groups).to(dev)
        def fwd(bi, net=net):
            xw, mw = M.pad_batch(bags, bi, dev)
            return net(xw, mw, Xo_t[bi], Xc_t[bi])[0]
        def step(b, fwd=fwd, net=net): net.train(); return M.nll_loss(fwd(fit_i[b]), k_t[fit_i[b]], c_t[fit_i[b]])
        def val(fwd=fwd, net=net):
            net.eval()
            with torch.no_grad():
                return float(np.mean([float(M.nll_loss(fwd(va_i[b:b + 32]), k_t[va_i[b:b + 32]], c_t[va_i[b:b + 32]]))
                                      for b in range(0, len(va_i), 32)]))
        d = fit(step, list(net.parameters()), len(fit_i), val, dev, seed)
        diag.append({**tag, "model": name, **d})
        net.eval()
        with torch.no_grad():
            h = np.concatenate([torch.sigmoid(fwd(pred_i[b:b + 32])).cpu().numpy() for b in range(0, len(pred_i), 32)])
        out[name] = ("S", np.cumprod(1 - h, 1))
    return out, ex_out


# ============================================================================ outer loop
def cox_fit_predict(Xc, time_, event, tr, te, grid=None):
    from lifelines import CoxPHFitter
    df = pd.DataFrame(Xc, columns=[f"c{i}" for i in range(Xc.shape[1])])
    keep = [c for c in df.columns if df.iloc[tr][c].std() > 0]
    df = df[keep]; df["T"], df["E"] = time_, event
    cph = CoxPHFitter(penalizer=0.1).fit(df.iloc[tr], "T", "E")
    X = df.iloc[te].drop(columns=["T", "E"])
    lp = cph.predict_log_partial_hazard(X).to_numpy()
    S = cph.predict_survival_function(X, times=grid).to_numpy().T if grid is not None else None
    return lp, S


def coxnet_fit_predict(Xo, time_, event, tr, te, grid, inner_splits):
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sksurv.util import Surv
    y = Surv.from_arrays(event.astype(bool), time_)
    base = CoxnetSurvivalAnalysis(l1_ratio=0.5, n_alphas=30, alpha_min_ratio=0.05, max_iter=20000).fit(Xo[tr], y[tr])
    alphas = base.alphas_
    scores = np.zeros(len(alphas))
    for itr, ite in inner_splits:
        mdl = CoxnetSurvivalAnalysis(l1_ratio=0.5, alphas=alphas, max_iter=20000).fit(Xo[itr], y[itr])
        for j, a in enumerate(mdl.alphas_):
            scores[j] += M.harrell_c(time_[ite], event[ite], mdl.predict(Xo[ite], alpha=a))
    a_best = alphas[int(np.argmax(scores))]
    mdl = CoxnetSurvivalAnalysis(l1_ratio=0.5, alphas=[a_best], fit_baseline_model=True, max_iter=20000).fit(Xo[tr], y[tr])
    risk = mdl.predict(Xo[te])
    fns = mdl.predict_survival_function(Xo[te])
    g = np.clip(grid, fns[0].domain[0], fns[0].domain[1])
    S = np.stack([f(g) for f in fns])
    return risk, S, a_best


def main_loop(args, log):
    import torch
    from sklearn.model_selection import StratifiedKFold, train_test_split
    from lifelines import CoxPHFitter
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    seeds = [0] if args.smoke else [0, 1, 2]
    for c in args.cohorts.split(","):
        od = M.res_dir(args) / f"v2_{c}_{args.encoder}_cap{args.cap}"
        od.mkdir(parents=True, exist_ok=True)
        tab, hm, mut, cna = M.load_prep(args, c)
        time_, event = tab.time.to_numpy(float), tab.event.to_numpy(int)
        n = len(tab)
        log(f"[v2 {c}] {n} patients, {event.sum()} deaths, encoder {args.encoder}, tile cap {args.cap}")
        preds, fms, diag, mech, weights, inner_rows = [], [], [], [], [], []
        for seed in seeds:
            bags, _ = M.load_bags(args, c, args.encoder, list(tab.patient), args.cap, seed)
            skf = StratifiedKFold(n_splits=args.folds, shuffle=True, random_state=seed)
            for f, (tr_all, te) in enumerate(skf.split(np.zeros(n), event)):
                t0 = time.time()
                edges = M.make_bins(time_, event, tr_all)
                tau = float(np.quantile(time_[tr_all], 0.8))
                grid = np.linspace(np.quantile(time_[tr_all], 0.1), tau, 30)
                parts = M.omics_matrix(tab, hm, mut, cna, tr_all)
                keys = ("hallmark", "mutation", "cna", "burden")
                Xo = np.concatenate([parts[k] for k in keys], 1).astype(np.float32)
                sizes = list(np.cumsum([0] + [parts[k].shape[1] for k in keys]))
                Xc = M.clin_matrix(tab, tr_all)
                inner = list(StratifiedKFold(n_splits=INNER, shuffle=True, random_state=1000 + seed)
                             .split(np.zeros(len(tr_all)), event[tr_all]))
                inner_abs = [(tr_all[a], tr_all[b]) for a, b in inner]
                # clinical anchor: out-of-fold for training patients, full-training fit for test patients
                anchor = np.zeros(n)
                for itr, ite in inner_abs:
                    anchor[ite] = cox_fit_predict(Xc, time_, event, itr, ite)[0]
                lp_te, S_cox = cox_fit_predict(Xc, time_, event, tr_all, te, grid)
                anchor[te] = lp_te
                data = dict(bags=bags, Xo=Xo, Xc=Xc, kbin=M.to_bin(time_, edges), cens=1 - event,
                            n_bins=len(edges) - 1, group_sizes=sizes)
                te_acc, oof = {}, {}
                ex_acc = {}
                for g, (itr, ite) in enumerate(inner_abs):
                    fit_i, va_i = train_test_split(itr, test_size=0.2, stratify=event[itr], random_state=seed + g)
                    pred_i = np.concatenate([ite, te])
                    tag = {"seed": seed, "fold": f, "inner": g}
                    out, ex = run_inner(data, fit_i, va_i, pred_i, anchor, dev, seed * 100 + f * 10 + g, diag, tag)
                    for m, v in out.items():
                        S = v[1] if isinstance(v, tuple) else eta_to_S(v)
                        oof.setdefault(m, np.full((n, S.shape[1]), np.nan))[ite] = S[:len(ite)]
                        te_acc.setdefault(m, []).append(S[len(ite):])
                    for m, rec in ex.items():
                        ex_acc.setdefault(m, []).append(rec)
                S_te = {m: np.mean(v, 0) for m, v in te_acc.items()}
                risk = {m: -S.sum(1) for m, S in S_te.items()}
                Sg = {m: M.surv_at(S, edges, grid) for m, S in S_te.items()}
                risk["Clinical-Cox"], Sg["Clinical-Cox"] = lp_te, S_cox
                # omics CoxNet (alpha chosen on the inner folds)
                r_cn, S_cn, a_cn = coxnet_fit_predict(Xo, time_, event, tr_all, te, grid, inner_abs)
                risk["Omics-CoxNet"], Sg["Omics-CoxNet"] = r_cn, S_cn
                # stacked Cox: clinical covariates + out-of-fold unimodal risks
                r_oof = {m: -oof[m][tr_all].sum(1) for m in ("WSI-only", "Omics-only")}
                mu = {m: (r_oof[m].mean(), r_oof[m].std() + 1e-8) for m in r_oof}
                def stack_X(idx, rr):
                    return np.column_stack([Xc[idx]] + [(rr[m] - mu[m][0]) / mu[m][1] for m in ("WSI-only", "Omics-only")])
                Xs_tr = stack_X(tr_all, r_oof)
                Xs_te = stack_X(te, {m: risk[m] for m in ("WSI-only", "Omics-only")})
                Xs = np.zeros((n, Xs_tr.shape[1])); Xs[tr_all] = Xs_tr; Xs[te] = Xs_te
                risk["Stacked-Cox"], Sg["Stacked-Cox"] = cox_fit_predict(Xs, time_, event, tr_all, te, grid)
                # late fusion: average of within-fold ranks of WSI-only, Omics-only and clinical Cox
                risk["Late-Fusion"] = np.mean([pd.Series(risk[m]).rank(pct=True) for m in ("WSI-only", "Omics-only", "Clinical-Cox")], 0)
                Sg["Late-Fusion"] = None
                for m, r in risk.items():
                    mt = M.fold_metrics(time_, event, r, time_[tr_all], event[tr_all], te, Sg.get(m), grid, tau)
                    fms.append({"seed": seed, "fold": f, "method": m, **mt})
                    rr = pd.Series(r).rank(pct=True).to_numpy()
                    for i, p in enumerate(te):
                        preds.append({"seed": seed, "fold": f, "patient": tab.patient[p], "method": m,
                                      "risk": r[i], "risk_rank_in_fold": rr[i]})
                for i, p in enumerate(tr_all):
                    inner_rows.append({"seed": seed, "fold": f, "role": "train", "patient": tab.patient[p],
                                       "risk_wsi": r_oof["WSI-only"][i], "risk_omics": r_oof["Omics-only"][i]})
                for i, p in enumerate(te):
                    inner_rows.append({"seed": seed, "fold": f, "role": "test", "patient": tab.patient[p],
                                       "risk_wsi": risk["WSI-only"][i], "risk_omics": risk["Omics-only"][i]})
                for m, recs in ex_acc.items():                                  # mechanism, test part only
                    nte = len(te)
                    a = np.mean([r["alpha"][-nte:] for r in recs], 0)
                    ent = np.mean([r["entropy"][-nte:] for r in recs], 0)
                    mech.append({"seed": seed, "fold": f, "method": m, "alpha_wsi": a[:, 0].mean(),
                                 "alpha_omics": a[:, 1].mean(), "alpha_clin": a[:, 2].mean(), "entropy": ent.mean(),
                                 **{k: np.mean([r[k] for r in recs]) for k in recs[0] if k.startswith("dNLL")},
                                 **{f"norm_{nm}": np.mean([r["norms"][i] for r in recs]) for i, nm in enumerate(("wsi", "omics", "clin"))}})
                    for i, p in enumerate(te):
                        w = {"seed": seed, "fold": f, "method": m, "patient": tab.patient[p],
                             "alpha_wsi": a[i, 0], "alpha_omics": a[i, 1], "alpha_clin": a[i, 2]}
                        if "delta" in recs[0]:
                            dl = np.mean([r["delta"][-nte:] for r in recs], 0)
                            w.update({"delta_wsi_omics": dl[i, 0, 1], "delta_wsi_clin": dl[i, 0, 2], "delta_omics_clin": dl[i, 1, 2]})
                        weights.append(w)
                cd = [x["harrell_c"] for x in fms if x["seed"] == seed and x["fold"] == f and x["method"] in ("Cox+DCA", "Clinical-Cox")]
                log(f"[v2 {c}] seed {seed} fold {f + 1}/{args.folds} done, {time.time() - t0:.0f}s; "
                    f"Cox C={cd[1] if len(cd) > 1 else float('nan'):.3f}，Cox+DCA C={cd[0]:.3f}，CoxNet alpha={a_cn:.3g}")
        pr, fm = pd.DataFrame(preds), pd.DataFrame(fms)
        pr.to_csv(od / "oof_predictions.csv", index=False); fm.to_csv(od / "fold_metrics.csv", index=False)
        pd.DataFrame(diag).to_csv(od / "diagnostics.csv", index=False)
        pd.DataFrame(mech).to_csv(od / "fusion_mechanism.csv", index=False)
        pd.DataFrame(weights).to_csv(od / "fusion_weights.csv", index=False)
        pd.DataFrame(inner_rows).to_csv(od / "unimodal_risks_inner.csv", index=False)
        summarize(tab, pr, fm, pd.DataFrame(diag), pd.DataFrame(mech), od, args, log, c)


# ============================================================================ summary
def summarize(tab, pr, fm, dg, mech, od, args, log, c):
    time_, event = tab.time.to_numpy(float), tab.event.to_numpy(int)
    pooled = pr.groupby(["method", "patient"])["risk_rank_in_fold"].mean().unstack("method").reindex(tab.patient)
    n_boot = 200 if args.smoke else 2000
    rows = []
    for m in [x for x in ORDER if x in fm.method.unique()]:
        g = fm[fm.method == m]
        r = {"method": m}
        for k in ("harrell_c", "uno_c", "ibs", "brier_12m", "brier_24m"):
            v = g[k].dropna() if k in g else pd.Series(dtype=float)
            r[k + "_mean"], r[k + "_sd"] = (v.mean(), v.std()) if len(v) else (np.nan, np.nan)
        r["pooled_c"] = M.harrell_c(time_, event, pooled[m].to_numpy())
        rows.append(r)
    s = pd.DataFrame(rows); s.to_csv(od / "summary.csv", index=False)
    con = []
    for a, b in CONTRASTS:
        if a in pooled and b in pooled:
            d, lo, hi, p = M.bootstrap_diff(time_, event, pooled[a].to_numpy(), pooled[b].to_numpy(), n_boot)
            con.append({"contrast": f"{a} − {b}", "diff_c": d, "ci_lo": lo, "ci_hi": hi, "share_boot_le_0": p})
    cn = pd.DataFrame(con); cn.to_csv(od / "contrasts.csv", index=False)
    L = [f"# Experiment 2, protocol v2 — {c} (encoder {args.encoder}, tiles <= {args.cap})",
         f"Patients {len(tab)}, deaths {int(event.sum())}, censoring {1 - event.mean():.0%}; outer {args.folds}-fold x "
         f"{fm.seed.nunique()} seeds; inner {INNER}-fold ensembles.", "",
         "| Method | Harrell C (mean ± SD) | Uno C | IBS | Brier 12 m | Brier 24 m | Pooled C |", "|---|---|---|---|---|---|---|"]
    fmt = lambda v: "—" if pd.isna(v) else f"{v:.3f}"
    for _, r in s.iterrows():
        L.append(f"| {r.method} | {r.harrell_c_mean:.3f} ± {r.harrell_c_sd:.3f} | {fmt(r.uno_c_mean)} | {fmt(r.ibs_mean)} | "
                 f"{fmt(r.brier_12m_mean)} | {fmt(r.brier_24m_mean)} | {r.pooled_c:.3f} |")
    L += ["", f"Primary contrast (pre-specified): {PRIMARY[0]} − {PRIMARY[1]}, pooled out-of-fold C, patient-level bootstrap.", "",
          "| Contrast | ΔC | 95% CI | Share of bootstrap ≤ 0 |", "|---|---|---|---|"]
    for _, r in cn.iterrows():
        L.append(f"| {r.contrast} | {r.diff_c:+.3f} | [{r.ci_lo:+.3f}, {r.ci_hi:+.3f}] | {r.share_boot_le_0:.3f} |")
    if len(dg):
        g = dg.groupby("model")[["best_epoch", "stop_epoch", "best_val_nll", "last_train_nll", "mean_grad_norm", "seconds"]].median()
        L += ["", "Training diagnostics (median over all inner models):", "",
              "| Model | Best epoch | Stop epoch | Best val NLL | Last train NLL | Grad norm | Seconds |", "|---|---|---|---|---|---|---|"]
        for mname, r in g.iterrows():
            L.append(f"| {mname} | {r.best_epoch:.0f} | {r.stop_epoch:.0f} | {r.best_val_nll:.3f} | {r.last_train_nll:.3f} | "
                     f"{r.mean_grad_norm:.2f} | {r.seconds:.1f} |")
    if len(mech):
        g = mech.groupby("method").mean(numeric_only=True)
        cols = [x for x in g.columns if x not in ("seed", "fold")]
        L += ["", "Fusion mechanism (mean over folds; dNLL_remove_* = rise in test NLL when that modality is replaced by its training mean):", "",
              "| Method | " + " | ".join(cols) + " |", "|---" * (len(cols) + 1) + "|"]
        for mname, r in g.iterrows():
            L.append(f"| {mname} | " + " | ".join(f"{r[x]:.3f}" for x in cols) + " |")
    (od / "summary.md").write_text("\n".join(L), encoding="utf-8")
    log("\n".join(L))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohorts", default="STAD,CRC")
    ap.add_argument("--encoder", default="uni")
    ap.add_argument("--cap", type=int, default=500)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    global EPOCHS, PATIENCE
    if args.smoke:
        args.folds, args.cap = 3, 100
        EPOCHS, PATIENCE = 3, 3
    log = M.Log(M.res_dir(args))
    log(f"===== exp2 v2 cohorts={args.cohorts} encoder={args.encoder} cap={args.cap} smoke={args.smoke} =====")
    main_loop(args, log)
    log("===== done =====")


if __name__ == "__main__":
    main()
