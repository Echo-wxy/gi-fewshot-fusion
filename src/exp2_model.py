"""Study 2, shared data preparation, encoders and metrics (also runs the MSI/hypermutation tile check, task A) - pathology tiles + omics + clinical, conflict-aware fusion for overall survival.

Steps (each reuses the previous step's files):
  prep   : align patients across tiles / survival / expression / CNA / mutation; Hallmark ssGSEA (all genes);
           mutation and CNA matrices; clinical table  -> results/exp2/prep/<COHORT>_*.{csv,npz}
  taskA  : sanity check - attention-MIL predicts the Kather MSIMUT label from tiles (5-fold AUC)
  taskB  : main survival experiment, outer 5-fold x 3 seeds, every learned preprocessing step fitted on the
           outer-training patients only; 12 methods incl. conflict-aware fusion (DCA); metrics + paired bootstrap
  all    : prep, taskA, taskB
Usage: python exp2_model.py all [--cohorts STAD,CRC] [--encoder uni] [--smoke]
"""
import argparse, json, math, os, re, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("GIFS_ROOT", str(Path.home() / "gifs_data")))
STUDY = {"STAD": "stad_tcga_pan_can_atlas_2018", "CRC": "coadread_tcga_pan_can_atlas_2018"}
NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
          "In_Frame_Ins", "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}
METHODS = ["Clinical-Cox", "WSI-ABMIL", "Omics-MLP", "Clin+Omics-MLP", "Clin+WSI", "Concat-Fusion",
           "Late-Fusion", "Attention-Fusion", "Gated-Fusion(param-matched)", "PORPOISE-style", "MCAT-style",
           "DCA (ours)"]
OURS = "DCA (ours)"
ABLATIONS = ["DCA-shuffled-conflict", "DCA-random-conflict"]     # conflict-score ablations (batch 2)
ALL_METHODS = METHODS + ABLATIONS
DCA_LIKE = (OURS,) + tuple(ABLATIONS)


# ============================================================================ utilities
class Log:
    def __init__(self, d):
        d.mkdir(parents=True, exist_ok=True)
        self.f = open(d / "log.txt", "a", encoding="utf-8")

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        self.f.write(time.strftime("[%H:%M:%S] ") + s + "\n"); self.f.flush()


def res_dir(args):
    d = ROOT / "results" / ("exp2_smoke" if args.smoke else "exp2")
    d.mkdir(parents=True, exist_ok=True)
    return d


def feat_dir():
    """Tile features always live in results/exp2/features (written by exp2_extract.py), also for smoke runs."""
    return ROOT / "results" / "exp2" / "features"


def read_cbio(path, **kw):
    return pd.read_csv(path, sep="\t", comment="#", low_memory=False, **kw)


def set_seed(s):
    import random, torch
    random.seed(s); np.random.seed(s); torch.manual_seed(s)


# ============================================================================ step: prep
def load_hallmark(rd, log):
    """MSigDB Hallmark 2020 gene sets (50). Downloaded once through gseapy and cached as json."""
    cache = rd / "prep" / "hallmark_2020.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    import gseapy
    lib = gseapy.get_library("MSigDB_Hallmark_2020")
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(lib), encoding="utf-8")
    log(f"[prep] Hallmark gene sets downloaded and cached: {len(lib)}")
    return lib


def ssgsea(expr, gene_sets, alpha=0.25, min_size=10):
    """Single-sample GSEA (Barbie et al. 2009): per sample, genes ranked by expression; enrichment score =
    sum over positions of (weighted in-set ECDF - out-of-set ECDF). Purely per-sample, so no information
    crosses between patients. expr: genes x samples DataFrame (log2 scale)."""
    genes = expr.index.to_numpy()
    X = expr.to_numpy(dtype=np.float64)
    n_g = len(genes)
    order = np.argsort(-X, axis=0, kind="stable")                    # descending rank per sample
    rank_w = np.empty_like(X)
    for j in range(X.shape[1]):
        r = np.empty(n_g)
        r[order[:, j]] = np.arange(n_g, 0, -1)                          # top gene gets the largest rank value
        rank_w[:, j] = r
    out, names = [], []
    gidx = {g: i for i, g in enumerate(genes)}
    for name, gs in gene_sets.items():
        idx = np.array(sorted({gidx[g] for g in gs if g in gidx}))
        if len(idx) < min_size:
            continue
        member = np.zeros(n_g, bool); member[idx] = True
        scores = np.empty(X.shape[1])
        for j in range(X.shape[1]):
            o = order[:, j]
            m = member[o]
            w = np.where(m, rank_w[o, j] ** alpha, 0.0)
            p_in = np.cumsum(w) / w.sum()
            p_out = np.cumsum(~m) / (~m).sum()
            scores[j] = (p_in - p_out).sum()
        out.append(scores); names.append(name)
    return pd.DataFrame(np.array(out).T, index=expr.columns, columns=names)


def stage_group(s):
    if not isinstance(s, str):
        return "missing"
    m = re.search(r"STAGE\s+(IV|III|II|I)", s.upper())
    return m.group(1) if m else "missing"


def step_prep(args, log):
    import h5py
    rd = res_dir(args)
    pdir = rd / "prep"; pdir.mkdir(parents=True, exist_ok=True)
    hallmark = load_hallmark(rd, log)
    for c in args.cohorts.split(","):
        out_tab = pdir / f"{c}_patients.csv"
        if out_tab.exists() and not args.force:
            log(f"[prep {c}] exists, skipping"); continue
        d = ROOT / "data" / "tcga" / STUDY[c]
        pat = read_cbio(d / "data_clinical_patient.txt")
        smp = read_cbio(d / "data_clinical_sample.txt")
        # tiles present (from the main encoder's feature file; all encoders hold the same tiles)
        h = h5py.File(feat_dir() / f"{c}_{args.encoder}.h5", "r")
        tp = pd.Series(h["patient"][:].astype(str))
        zips = pd.Series(h["source_zip"][:].astype(str))
        tile_n = tp.value_counts()
        msimut = pd.Series(zips.str.contains("MSIMUT").to_numpy()).groupby(tp.to_numpy()).max()
        h.close()
        # omics
        expr = read_cbio(d / "data_mrna_seq_v2_rsem.txt")
        expr = expr.dropna(subset=["Hugo_Symbol"]).drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
        expr = expr[[x for x in expr.columns if x.startswith("TCGA-")]]
        expr.columns = [x[:12] for x in expr.columns]
        expr = expr.loc[:, ~expr.columns.duplicated()]
        cna = read_cbio(d / "data_cna.txt")
        cna = cna.dropna(subset=["Hugo_Symbol"]).drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
        cna = cna[[x for x in cna.columns if x.startswith("TCGA-")]]
        cna.columns = [x[:12] for x in cna.columns]
        cna = cna.loc[:, ~cna.columns.duplicated()]
        mut = read_cbio(d / "data_mutations.txt", usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"])
        mut["pid"] = mut["Tumor_Sample_Barcode"].astype(str).str[:12]
        mut_pat = set(mut["pid"])
        mut = mut[mut["Variant_Classification"].isin(NONSYN)]
        # cohort
        pat = pat.dropna(subset=["OS_MONTHS", "OS_STATUS"])
        keep = (set(pat.PATIENT_ID) & set(tile_n.index) & set(expr.columns) & set(cna.columns) & mut_pat)
        keep = sorted(keep)
        log(f"[prep {c}] patients with all modalities: {len(keep)}")
        P = pat.set_index("PATIENT_ID").loc[keep]
        S = smp.drop_duplicates("PATIENT_ID").set_index("PATIENT_ID").reindex(keep)
        tab = pd.DataFrame({
            "patient": keep,
            "time": P["OS_MONTHS"].astype(float).clip(lower=0.1).to_numpy(),
            "event": P["OS_STATUS"].astype(str).str.startswith("1").astype(int).to_numpy(),
            "age": pd.to_numeric(P["AGE"], errors="coerce").to_numpy(),
            "male": (P["SEX"].astype(str).str.upper() == "MALE").astype(int).to_numpy(),
            "stage": [stage_group(s) for s in P["AJCC_PATHOLOGIC_TUMOR_STAGE"]],
            "subtype": P["SUBTYPE"].astype(str).to_numpy() if "SUBTYPE" in P else "nan",
            "tmb": pd.to_numeric(S["TMB_NONSYNONYMOUS"], errors="coerce").to_numpy(),
            "msi_mantis": pd.to_numeric(S["MSI_SCORE_MANTIS"], errors="coerce").to_numpy(),
            "msi_sensor": pd.to_numeric(S["MSI_SENSOR_SCORE"], errors="coerce").to_numpy(),
            "n_tiles": tile_n.reindex(keep).to_numpy(),
            "kather_msimut": msimut.reindex(keep).astype(int).to_numpy(),
        })
        tab.to_csv(out_tab, index=False)
        # Hallmark ssGSEA on all genes (log2(x+1))
        E = np.log2(expr[keep].astype(float).clip(lower=0) + 1)
        hm = ssgsea(E, hallmark)
        hm.to_csv(pdir / f"{c}_hallmark.csv")
        # mutation matrix (genes mutated in >= 1% of patients; final top-k chosen inside each training fold)
        mm = pd.crosstab(mut["pid"], mut["Hugo_Symbol"]).reindex(keep).fillna(0).clip(upper=1).astype(np.int8)
        mm = mm.loc[:, mm.mean(0) >= 0.01]
        cn = cna[keep].T.astype(float).fillna(0)
        cn = cn.loc[:, cn.std(0) > 0]
        np.savez_compressed(pdir / f"{c}_omics.npz", mut=mm.to_numpy(), mut_genes=mm.columns.to_numpy(str),
                            cna=cn.to_numpy(np.float32), cna_genes=cn.columns.to_numpy(str),
                            patients=np.array(keep))
        log(f"[prep {c}] {len(keep)} patients, {int(tab.event.sum())} deaths; {hm.shape[1]} Hallmark sets; "
            f"{mm.shape[1]} mutation genes; {cn.shape[1]} CNA genes; MSIMUT {int(tab.kather_msimut.sum())}")


def load_prep(args, c):
    pdir = res_dir(args) / "prep"
    tab = pd.read_csv(pdir / f"{c}_patients.csv")
    hm = pd.read_csv(pdir / f"{c}_hallmark.csv", index_col=0).reindex(tab.patient)
    om = np.load(pdir / f"{c}_omics.npz", allow_pickle=False)
    assert list(om["patients"]) == list(tab.patient)
    hm, mut, cna = hm.to_numpy(np.float32), om["mut"].astype(np.float32), om["cna"].astype(np.float32)
    if args.smoke:                       # rehearsal: 80 patients, chosen once with a fixed seed
        keep = np.sort(np.random.default_rng(0).choice(len(tab), min(80, len(tab)), replace=False))
        tab, hm, mut, cna = tab.iloc[keep].reset_index(drop=True), hm[keep], mut[keep], cna[keep]
    return tab, hm, mut, cna


def load_bags(args, c, enc, patients, cap, seed):
    """Per patient: up to `cap` tile features (random subset with a fixed seed). Returns list of float16 arrays
    and the tile indices used (saved so that runs can be repeated)."""
    import h5py
    h = h5py.File(feat_dir() / f"{c}_{enc}.h5", "r")
    tp = h["patient"][:].astype(str)
    feats = h["feats"][:]
    h.close()
    rng = np.random.default_rng(seed)
    idx_by = pd.Series(np.arange(len(tp))).groupby(tp).apply(np.array).to_dict()
    bags, used = [], {}
    for p in patients:
        idx = idx_by[p]
        if cap and len(idx) > cap:
            idx = np.sort(rng.choice(idx, cap, replace=False))
        bags.append(feats[idx].astype(np.float16)); used[p] = idx
    return bags, used


# ============================================================================ fold-level feature building
def clin_matrix(tab, tr):
    age = tab.age.to_numpy(float)
    med, sd = np.nanmedian(age[tr]), np.nanstd(age[tr]) + 1e-6
    a = np.where(np.isnan(age), med, age)
    cols = [(a - med) / sd, np.isnan(age).astype(float), tab.male.to_numpy(float)]
    for st in ["I", "II", "III", "IV", "missing"]:
        cols.append((tab.stage == st).to_numpy(float))
    return np.stack(cols, 1).astype(np.float32)


def omics_matrix(tab, hm, mut, cna, tr, k_mut=50, k_cna=100):
    """Every selection / scaling fitted on training patients `tr` only."""
    def zs(x):
        mu, sd = np.nanmean(x[tr], 0), np.nanstd(x[tr], 0) + 1e-6
        x = np.where(np.isnan(x), mu, x)
        return (x - mu) / sd
    top_m = np.argsort(-mut[tr].mean(0), kind="stable")[:k_mut]
    top_c = np.argsort(-cna[tr].var(0), kind="stable")[:k_cna]
    sc = np.stack([np.log1p(tab.tmb.to_numpy(float)), tab.msi_mantis.to_numpy(float),
                   tab.msi_sensor.to_numpy(float)], 1)
    parts = {"hallmark": zs(hm), "mutation": mut[:, top_m], "cna": zs(cna[:, top_c]), "burden": zs(sc)}
    return parts


def make_bins(time, event, tr, n_bins=4):
    q = np.quantile(time[tr][event[tr] == 1], np.linspace(0, 1, n_bins + 1))
    q[0], q[-1] = 0.0, np.inf
    return q


def to_bin(time, edges):
    return np.clip(np.searchsorted(edges, time, side="right") - 1, 0, len(edges) - 2)


# ============================================================================ models
def build_model(name, d_wsi, d_omics, d_clin, n_bins, groups):
    import torch
    import torch.nn as nn

    def mlp(i, o, h=128, p=0.25):
        return nn.Sequential(nn.Linear(i, h), nn.ReLU(), nn.Dropout(p), nn.Linear(h, o))

    class ABMIL(nn.Module):        # gated attention MIL (Ilse et al. 2018)
        def __init__(self, d, h=128):
            super().__init__()
            self.fc = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(0.25))
            self.a, self.b, self.w = nn.Linear(h, 64), nn.Linear(h, 64), nn.Linear(64, 1)

        def forward(self, x, m):     # x: B x N x d, m: B x N (1 = real tile)
            h = self.fc(x)
            s = self.w(torch.tanh(self.a(h)) * torch.sigmoid(self.b(h))).squeeze(-1)
            s = s.masked_fill(m == 0, -1e4)
            att = torch.softmax(s, 1)
            return (att.unsqueeze(-1) * h).sum(1)

    D = 128

    class Net(nn.Module):
        def __init__(self):
            super().__init__()
            self.kind = name
            if name in ("WSI-ABMIL", "Clin+WSI", "Concat-Fusion", "Attention-Fusion", "Gated-Fusion(param-matched)",
                        "PORPOISE-style") + DCA_LIKE:
                self.wsi = ABMIL(d_wsi, D)
            if name in ("Omics-MLP", "Clin+Omics-MLP", "Concat-Fusion", "Attention-Fusion",
                        "Gated-Fusion(param-matched)", "PORPOISE-style") + DCA_LIKE:
                self.om = mlp(d_omics, D)
            if name in ("Clin+WSI", "Concat-Fusion", "Attention-Fusion", "Gated-Fusion(param-matched)") + DCA_LIKE:
                self.cl = mlp(d_clin, D, h=32)
            if name == "Clin+Omics-MLP":
                self.om = mlp(d_omics + d_clin, D)
            if name in ("Attention-Fusion", "Gated-Fusion(param-matched)") + DCA_LIKE:
                self.norm = nn.LayerNorm(D)
                self.gate = nn.Linear(D, 1)
            if name in DCA_LIKE:
                # pairwise conflict score from [z_i - z_j, z_i * z_j]; penalty strength gamma >= 0 via softplus
                self.conf = nn.Sequential(nn.Linear(2 * D, 64), nn.ReLU(), nn.Linear(64, 1))
                self.gamma_raw = nn.Parameter(torch.tensor(0.0))
            if name == "Gated-Fusion(param-matched)":
                # same extra capacity as the conflict MLP, but sees [z_m, mean of the other modalities] only
                self.ctx = nn.Sequential(nn.Linear(2 * D, 64), nn.ReLU(), nn.Linear(64, 1))
            if name == "PORPOISE-style":    # gated modality embeddings + Kronecker (outer-product) fusion
                self.gp, self.go = nn.Linear(2 * D, D), nn.Linear(2 * D, D)
                self.kp, self.ko = nn.Linear(D, 32), nn.Linear(D, 32)
                self.post = nn.Sequential(nn.Linear(33 * 33, D), nn.ReLU(), nn.Dropout(0.25))
            if name == "MCAT-style":        # genomic-group tokens attend over tiles (co-attention), then pooled
                self.tok = nn.ModuleList([mlp(len(g), D, h=64) for g in groups])
                self.tile_proj = nn.Linear(d_wsi, D)
                self.coattn = nn.MultiheadAttention(D, 4, batch_first=True, dropout=0.1)
                self.pool_w = nn.Linear(D, 1)
            in_head = 2 * D if name == "MCAT-style" else (3 * D if name == "Concat-Fusion" else
                                                          (2 * D if name == "Clin+WSI" else D))
            self.head = nn.Sequential(nn.Linear(in_head, 64), nn.ReLU(), nn.Dropout(0.25), nn.Linear(64, n_bins))

        def forward(self, xw, mw, xo, xc):
            extra = {}
            if self.kind == "WSI-ABMIL":
                f = self.wsi(xw, mw)
            elif self.kind == "Omics-MLP":
                f = self.om(xo)
            elif self.kind == "Clin+Omics-MLP":
                f = self.om(torch.cat([xo, xc], 1))
            elif self.kind == "Clin+WSI":
                f = torch.cat([self.wsi(xw, mw), self.cl(xc)], 1)
            elif self.kind == "Concat-Fusion":
                f = torch.cat([self.wsi(xw, mw), self.om(xo), self.cl(xc)], 1)
            elif self.kind in ("Attention-Fusion", "Gated-Fusion(param-matched)") + DCA_LIKE:
                z = self.norm(torch.stack([self.wsi(xw, mw), self.om(xo), self.cl(xc)], 1))   # B x 3 x D
                q = self.gate(z).squeeze(-1)                                                   # B x 3
                if self.kind == "Gated-Fusion(param-matched)":
                    others = (z.sum(1, keepdim=True) - z) / 2
                    q = q + self.ctx(torch.cat([z, others], -1)).squeeze(-1)
                if self.kind in DCA_LIKE:
                    zi, zj = z.unsqueeze(2), z.unsqueeze(1)                                    # B x 3 x 1 x D
                    pair = torch.cat([(zi - zj).expand(-1, -1, 3, -1), (zi * zj).expand(-1, -1, 3, -1)], -1)
                    delta = torch.sigmoid(self.conf(pair).squeeze(-1))                         # B x 3 x 3
                    delta = 0.5 * (delta + delta.transpose(1, 2))                              # symmetric
                    eye = torch.eye(3, device=z.device).unsqueeze(0)
                    delta = delta * (1 - eye)
                    if self.kind == "DCA-shuffled-conflict":       # conflict scores of other patients in the batch
                        delta = delta[torch.randperm(delta.shape[0], device=delta.device)]
                    elif self.kind == "DCA-random-conflict":       # conflict scores unrelated to the patient
                        delta = torch.rand_like(delta) * (1 - eye)
                        delta = 0.5 * (delta + delta.transpose(1, 2))
                    gamma = torch.nn.functional.softplus(self.gamma_raw)
                    q = q - gamma * delta.sum(2)
                    extra["delta"] = delta
                alpha = torch.softmax(q, 1)
                extra["alpha"] = alpha
                f = (alpha.unsqueeze(-1) * z).sum(1)
            elif self.kind == "PORPOISE-style":
                hp, ho = self.wsi(xw, mw), self.om(xo)
                hp = hp * torch.sigmoid(self.gp(torch.cat([hp, ho], 1)))
                ho = ho * torch.sigmoid(self.go(torch.cat([hp, ho], 1)))
                one = torch.ones(hp.shape[0], 1, device=hp.device)
                a, b = torch.cat([self.kp(hp), one], 1), torch.cat([self.ko(ho), one], 1)
                f = self.post(torch.bmm(a.unsqueeze(2), b.unsqueeze(1)).flatten(1))
            elif self.kind == "MCAT-style":
                toks = torch.stack([t(xo[:, g]) for t, g in zip(self.tok, groups_t)], 1)       # B x G x D
                tiles = self.tile_proj(xw)
                att, _ = self.coattn(toks, tiles, tiles, key_padding_mask=(mw == 0))
                w = torch.softmax(self.pool_w(att).squeeze(-1), 1)
                f = torch.cat([(w.unsqueeze(-1) * att).sum(1), toks.mean(1)], 1)
            return self.head(f), extra

    groups_t = groups
    return Net()


def nll_loss(logits, k, c):
    """Discrete-time survival NLL (Zadeh & Schmid 2021); c = 1 if censored."""
    import torch
    h = torch.sigmoid(logits).clamp(1e-6, 1 - 1e-6)
    S = torch.cumprod(1 - h, 1)
    S_pad = torch.cat([torch.ones_like(S[:, :1]), S], 1)
    s_prev = S_pad.gather(1, k.unsqueeze(1)).squeeze(1)
    h_k = h.gather(1, k.unsqueeze(1)).squeeze(1)
    s_k = S.gather(1, k.unsqueeze(1)).squeeze(1)
    unc = -(1 - c) * (torch.log(s_prev) + torch.log(h_k))
    cen = -c * torch.log(s_k)
    return (unc + cen).mean()


def pad_batch(bags, idx, dev):
    import torch
    n = max(len(bags[i]) for i in idx)
    d = bags[idx[0]].shape[1]
    x = np.zeros((len(idx), n, d), np.float16); m = np.zeros((len(idx), n), np.float32)
    for r, i in enumerate(idx):
        x[r, :len(bags[i])] = bags[i]; m[r, :len(bags[i])] = 1
    return torch.tensor(x, device=dev).float(), torch.tensor(m, device=dev)


def harrell_c(time, event, risk):
    from lifelines.utils import concordance_index
    return concordance_index(time, -np.asarray(risk), event)


def train_deep(name, data, tr, va, te, args, seed, dev):
    """Train one deep model on `tr`, early-stop on validation C-index (`va`), predict `te`."""
    import copy, torch
    set_seed(seed)
    bags, Xo, Xc, kbin, cens, groups, n_bins = (data[k] for k in ("bags", "Xo", "Xc", "kbin", "cens", "groups", "n_bins"))
    model = build_model(name, bags[0].shape[1], Xo.shape[1], Xc.shape[1], n_bins, groups).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-4)
    Xo_t, Xc_t = torch.tensor(Xo, device=dev), torch.tensor(Xc, device=dev)
    k_t, c_t = torch.tensor(kbin, device=dev), torch.tensor(cens, device=dev, dtype=torch.float32)
    needs_wsi = hasattr(model, "wsi") or name == "MCAT-style"

    def predict(ids):
        model.eval()
        outs, extras = [], []
        with torch.no_grad():
            for b in range(0, len(ids), args.batch):
                bi = ids[b:b + args.batch]
                xw, mw = pad_batch(bags, bi, dev) if needs_wsi else (None, None)
                lo, ex = model(xw, mw, Xo_t[bi], Xc_t[bi])
                outs.append(torch.sigmoid(lo).cpu().numpy())
                extras.append({k: v.cpu().numpy() for k, v in ex.items()})
        h = np.concatenate(outs)
        S = np.cumprod(1 - h, 1)
        ex = {k: np.concatenate([e[k] for e in extras]) for k in (extras[0] if extras else {})}
        return S, ex

    best, best_state, bad = -1, copy.deepcopy(model.state_dict()), 0
    rng = np.random.default_rng(seed)
    for ep in range(args.epochs):
        model.train()
        order = rng.permutation(tr)
        for b in range(0, len(order), args.batch):
            bi = order[b:b + args.batch]
            xw, mw = pad_batch(bags, bi, dev) if needs_wsi else (None, None)
            lo, _ = model(xw, mw, Xo_t[bi], Xc_t[bi])
            loss = nll_loss(lo, k_t[bi], c_t[bi])
            opt.zero_grad(); loss.backward(); opt.step()
        S_va, _ = predict(va)
        cidx = harrell_c(data["time"][va], data["event"][va], -S_va.sum(1))
        if cidx > best:
            best, best_state, bad = cidx, copy.deepcopy(model.state_dict()), 0
        else:
            bad += 1
            if bad >= args.patience:
                break
    model.load_state_dict(best_state)
    S_te, ex = predict(te)
    return S_te, ex


# ============================================================================ metrics
def surv_at(S, edges, grid):
    """Step survival from bin-level survival S (n x K) evaluated at times `grid` (end-of-bin convention)."""
    ends = edges[1:].copy(); ends[-1] = np.inf
    pos = np.searchsorted(ends, grid, side="left")              # bin containing each grid time
    S_pad = np.concatenate([np.ones((S.shape[0], 1)), S], 1)
    return S_pad[:, np.clip(pos, 0, S.shape[1])]


def fold_metrics(time, event, risk, tr_time, tr_event, te, Sgrid, grid, tau):
    from sksurv.metrics import concordance_index_ipcw, integrated_brier_score, brier_score
    from sksurv.util import Surv
    y_tr = Surv.from_arrays(tr_event.astype(bool), tr_time)
    y_te = Surv.from_arrays(event[te].astype(bool), time[te])
    out = {"harrell_c": harrell_c(time[te], event[te], risk)}
    try:
        out["uno_c"] = concordance_index_ipcw(y_tr, y_te, risk, tau=tau)[0]
    except Exception:
        out["uno_c"] = np.nan
    if Sgrid is not None:
        ok = (grid > time[te].min()) & (grid < time[te].max())
        g = grid[ok]
        try:
            out["ibs"] = integrated_brier_score(y_tr, y_te, Sgrid[:, ok], g) if len(g) >= 2 else np.nan
        except Exception:
            out["ibs"] = np.nan
        for m in (12, 24, 36):
            try:
                j = int(np.argmin(np.abs(grid - m)))
                out[f"brier_{m}m"] = (brier_score(y_tr, y_te, Sgrid[:, [j]], [grid[j]])[1][0]
                                      if grid[0] <= m <= grid[-1] else np.nan)
            except Exception:
                out[f"brier_{m}m"] = np.nan
    return out


# ============================================================================ step: task A
def step_taskA(args, log):
    import torch
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    rd = res_dir(args)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    rows = []
    for c in args.cohorts.split(","):
        tab, *_ = load_prep(args, c)
        bags, _ = load_bags(args, c, args.encoder, list(tab.patient), args.cap, 0)
        y = tab.kather_msimut.to_numpy()
        skf = StratifiedKFold(n_splits=args.folds, shuffle=True, random_state=0)
        oof = np.zeros(len(y))
        for f, (tr, te) in enumerate(skf.split(np.zeros(len(y)), y)):
            set_seed(f)
            net = build_model("WSI-ABMIL", bags[0].shape[1], 1, 1, 1, []).to(dev)
            opt = torch.optim.Adam(net.parameters(), lr=args.lr, weight_decay=1e-4)
            pos_w = torch.tensor((y[tr] == 0).sum() / max(1, (y[tr] == 1).sum()), device=dev)
            rng = np.random.default_rng(f)
            for ep in range(max(3, args.epochs // 2)):
                net.train()
                order = rng.permutation(tr)
                for b in range(0, len(order), args.batch):
                    bi = order[b:b + args.batch]
                    xw, mw = pad_batch(bags, bi, dev)
                    lo, _ = net(xw, mw, None, None)
                    loss = torch.nn.functional.binary_cross_entropy_with_logits(
                        lo.squeeze(1), torch.tensor(y[bi], device=dev, dtype=torch.float32), pos_weight=pos_w)
                    opt.zero_grad(); loss.backward(); opt.step()
            net.eval()
            with torch.no_grad():
                for b in range(0, len(te), args.batch):
                    bi = te[b:b + args.batch]
                    xw, mw = pad_batch(bags, bi, dev)
                    oof[bi] = net(xw, mw, None, None)[0].squeeze(1).cpu().numpy()
            rows.append({"cohort": c, "fold": f, "auc": roc_auc_score(y[te], oof[te]) if len(set(y[te])) > 1 else np.nan})
        rows.append({"cohort": c, "fold": "pooled", "auc": roc_auc_score(y, oof)})
        log(f"[taskA {c}] MSIMUT {args.folds}-fold AUC (pooled out-of-fold) = {rows[-1]['auc']:.3f}")
    pd.DataFrame(rows).to_csv(rd / f"taskA_msimut_auc_{args.encoder}.csv", index=False)


# ============================================================================ step: task B
def bootstrap_diff(time, event, r_a, r_b, n_boot, seed=0):
    rng = np.random.default_rng(seed)
    n = len(time); diffs = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        if event[i].sum() < 2:
            continue
        diffs.append(harrell_c(time[i], event[i], r_a[i]) - harrell_c(time[i], event[i], r_b[i]))
    d = np.array(diffs)
    return d.mean(), np.percentile(d, 2.5), np.percentile(d, 97.5), float((d <= 0).mean())


def step_taskB(args, log):
    import torch
    from sklearn.model_selection import StratifiedKFold, train_test_split
    from lifelines import CoxPHFitter
    rd = res_dir(args)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    seeds = [0] if args.smoke else [0, 1, 2]
    methods = [m for m in ALL_METHODS if m in args.methods.split(",")] if args.methods else list(METHODS)
    for c in args.cohorts.split(","):
        od = rd / f"taskB_{c}_{args.encoder}_cap{args.cap}" if not args.tag else rd / f"taskB_{c}_{args.encoder}_cap{args.cap}_{args.tag}"
        od.mkdir(parents=True, exist_ok=True)
        tab, hm, mut, cna = load_prep(args, c)
        time_, event = tab.time.to_numpy(float), tab.event.to_numpy(int)
        n = len(tab)
        log(f"[taskB {c}] {n} patients, {event.sum()} deaths, encoder {args.encoder}, tile cap {args.cap}, {len(methods)} methods")
        pred_rows, fm_rows, extra_rows = [], [], []
        for seed in seeds:
            bags, used = load_bags(args, c, args.encoder, list(tab.patient), args.cap, seed)
            skf = StratifiedKFold(n_splits=args.folds, shuffle=True, random_state=seed)
            for f, (tr_all, te) in enumerate(skf.split(np.zeros(n), event)):
                t0 = time.time()
                tr, va = train_test_split(tr_all, test_size=0.2, stratify=event[tr_all], random_state=seed)
                edges = make_bins(time_, event, tr_all)
                tau = float(np.quantile(time_[tr_all], 0.8))
                grid = np.linspace(np.quantile(time_[tr_all], 0.1), tau, 30)
                parts = omics_matrix(tab, hm, mut, cna, tr_all)
                Xo = np.concatenate([parts[k] for k in ("hallmark", "mutation", "cna", "burden")], 1).astype(np.float32)
                sizes = np.cumsum([0] + [parts[k].shape[1] for k in ("hallmark", "mutation", "cna", "burden")])
                groups = [list(range(sizes[i], sizes[i + 1])) for i in range(4)]
                Xc = clin_matrix(tab, tr_all)
                data = dict(bags=bags, Xo=Xo, Xc=Xc, kbin=to_bin(time_, edges), cens=1 - event,
                            groups=groups, n_bins=len(edges) - 1, time=time_, event=event)
                risks, Sg = {}, {}
                for m in methods:
                    if m == "Late-Fusion":
                        continue
                    if m == "Clinical-Cox":
                        df = pd.DataFrame(Xc, columns=[f"c{i}" for i in range(Xc.shape[1])])
                        df = df.loc[:, df.iloc[tr_all].std() > 0]
                        df["T"], df["E"] = time_, event
                        cph = CoxPHFitter(penalizer=0.1).fit(df.iloc[tr_all], "T", "E")
                        X_te = df.iloc[te].drop(columns=["T", "E"])
                        risks[m] = cph.predict_partial_hazard(X_te).to_numpy()
                        Sg[m] = cph.predict_survival_function(X_te, times=grid).to_numpy().T
                        continue
                    S, ex = train_deep(m, data, tr, va, te, args, seed * 100 + f, dev)
                    risks[m] = -S.sum(1)
                    Sg[m] = surv_at(S, edges, grid)
                    if m == OURS and ex:
                        for i, p in enumerate(te):
                            extra_rows.append({"seed": seed, "fold": f, "patient": tab.patient[p],
                                               "alpha_wsi": ex["alpha"][i, 0], "alpha_omics": ex["alpha"][i, 1],
                                               "alpha_clin": ex["alpha"][i, 2], "delta_wsi_omics": ex["delta"][i, 0, 1],
                                               "delta_wsi_clin": ex["delta"][i, 0, 2], "delta_omics_clin": ex["delta"][i, 1, 2]})
                if "Late-Fusion" in methods:
                    have = [k for k in ("WSI-ABMIL", "Omics-MLP", "Clinical-Cox") if k in risks]
                    if have:
                        rk = [pd.Series(risks[k]).rank(pct=True).to_numpy() for k in have]
                        risks["Late-Fusion"] = np.mean(rk, 0)
                        Sg["Late-Fusion"] = None
                for m, r in risks.items():
                    mt = fold_metrics(time_, event, r, time_[tr_all], event[tr_all], te, Sg.get(m), grid, tau)
                    fm_rows.append({"seed": seed, "fold": f, "method": m, **mt})
                    rr = pd.Series(r).rank(pct=True).to_numpy()            # within-fold rank, for pooling
                    for i, p in enumerate(te):
                        pred_rows.append({"seed": seed, "fold": f, "patient": tab.patient[p], "method": m,
                                          "risk": r[i], "risk_rank_in_fold": rr[i]})
                ours_c = [x["harrell_c"] for x in fm_rows if x["seed"] == seed and x["fold"] == f and x["method"] == OURS]
                log(f"[taskB {c}] seed {seed} fold {f + 1}/{args.folds} done, {time.time() - t0:.0f}s"
                    + (f"；DCA C={ours_c[0]:.3f}" if ours_c else ""))
        pr, fm = pd.DataFrame(pred_rows), pd.DataFrame(fm_rows)
        pr.to_csv(od / "oof_predictions.csv", index=False)
        fm.to_csv(od / "fold_metrics.csv", index=False)
        if extra_rows:
            pd.DataFrame(extra_rows).to_csv(od / "dca_weights.csv", index=False)
        pd.DataFrame([{"patient": p, "tiles_used": len(v)} for p, v in used.items()]).to_csv(od / "tiles_used_seed_last.csv", index=False)
        summarize_taskB(tab, pr, fm, od, args, log, c)


def summarize_taskB(tab, pr, fm, od, args, log, c):
    time_, event = tab.time.to_numpy(float), tab.event.to_numpy(int)
    pos = {p: i for i, p in enumerate(tab.patient)}
    # pooled out-of-fold risk: mean over seeds of the within-fold rank
    pooled = pr.groupby(["method", "patient"])["risk_rank_in_fold"].mean().unstack("method")
    pooled = pooled.reindex(tab.patient)
    n_boot = 200 if args.smoke else 2000
    rows = []
    for m in [x for x in ALL_METHODS if x in fm.method.unique()]:
        g = fm[fm.method == m]
        r = {"method": m}
        for k in ("harrell_c", "uno_c", "ibs", "brier_12m", "brier_24m", "brier_36m"):
            if k in g:
                v = g[k].dropna()
                r[k + "_mean"], r[k + "_sd"] = (v.mean(), v.std()) if len(v) else (np.nan, np.nan)
        r["pooled_c"] = harrell_c(time_, event, pooled[m].to_numpy())
        if m != OURS and OURS in pooled:
            d, lo, hi, p = bootstrap_diff(time_, event, pooled[OURS].to_numpy(), pooled[m].to_numpy(), n_boot)
            r.update({"ours_minus_this_c": d, "diff_ci_lo": lo, "diff_ci_hi": hi, "p_one_sided": p})
        rows.append(r)
    s = pd.DataFrame(rows)
    s.to_csv(od / "summary.csv", index=False)
    lines = [f"# Task B — {c}, encoder {args.encoder}, tiles per patient <= {args.cap}",
             f"Patients {len(tab)}, deaths {int(event.sum())}; outer {args.folds}-fold CV x {fm.seed.nunique()} seeds.", "",
             "| Method | Harrell C (mean ± SD over folds) | Uno C | IBS | Pooled C | DCA − this (95% bootstrap CI) |",
             "|---|---|---|---|---|---|"]
    for _, r in s.iterrows():
        diff = (f"{r.ours_minus_this_c:+.3f} [{r.diff_ci_lo:+.3f}, {r.diff_ci_hi:+.3f}]"
                if "ours_minus_this_c" in r and pd.notna(r.get("ours_minus_this_c")) else "—")
        ibs = f"{r.ibs_mean:.3f}" if pd.notna(r.get("ibs_mean", np.nan)) else "—"
        lines.append(f"| {r.method} | {r.harrell_c_mean:.3f} ± {r.harrell_c_sd:.3f} | {r.uno_c_mean:.3f} | {ibs} | "
                     f"{r.pooled_c:.3f} | {diff} |")
    lines += ["", "IBS: lower is better. Late-Fusion has no survival curve, so no IBS.",
              "Pooled C uses each patient's out-of-fold risk rank averaged over seeds."]
    (od / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    log("\n".join(lines))


# ============================================================================ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["prep", "taskA", "taskB", "all"])
    ap.add_argument("--cohorts", default="STAD,CRC")
    ap.add_argument("--encoder", default="uni")
    ap.add_argument("--cap", type=int, default=500, help="max tiles per patient (0 = all)")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--patience", type=int, default=10)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--lr", type=float, default=2e-4)
    ap.add_argument("--methods", default="", help="comma-separated subset of methods")
    ap.add_argument("--tag", default="", help="suffix for the task B output folder")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if args.smoke:
        args.folds, args.epochs, args.patience, args.cap = 3, 3, 3, 100
    log = Log(res_dir(args))
    log(f"===== exp2 {args.step} cohorts={args.cohorts} encoder={args.encoder} cap={args.cap} smoke={args.smoke} =====")
    steps = ["prep", "taskA", "taskB"] if args.step == "all" else [args.step]
    for s in steps:
        {"prep": step_prep, "taskA": step_taskA, "taskB": step_taskB}[s](args, log)
    log("===== done =====")


if __name__ == "__main__":
    main()
