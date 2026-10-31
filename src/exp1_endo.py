"""Study 1 - few-shot recognition of unseen-label endoscopic lesions (HyperKvasir).

Steps (run in order; each step reuses the previous step's output files):
  extract  : frozen DINOv2 ViT-B/14 features for every labeled image
  dedup    : exact hash + perceptual hash + embedding similarity -> near-duplicate clusters
  run      : episodic evaluation of all methods on all class splits, shots and both sampling modes
  figures  : summary tables, accuracy-vs-shot plot, occlusion attribution figure
  all      : the four steps above
Add --smoke for a few-minute rehearsal on a small subset (results go to results/exp1_smoke).

Usage:  python exp1_endo.py all [--smoke]
Root folder: env var GIFS_ROOT (default ~/gifs_data).
"""
import argparse, collections, csv, hashlib, json, math, os, random, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(os.environ.get("GIFS_ROOT", str(Path.home() / "gifs_data")))
ENDO = ROOT / "data" / "endoscopy"

ALIAS = {"hemorroids": "hemorrhoids", "normal-cecum": "cecum", "normal-pylorus": "pylorus",
         "normal-z-line": "z-line", "oesophagitis-a": "esophagitis-a", "oesophagitis-b-d": "esophagitis-b-d",
         "polyp": "polyps", "short-segment-barretts": "barretts-short-segment"}
MIN_IMAGES = 30                                    # classes with fewer images are dropped (4 classes)
MAIN_TEST = ["ulcerative-colitis-grade-0-1", "barretts", "barretts-short-segment",
             "ulcerative-colitis-grade-3", "esophagitis-b-d"]
UC_FAMILY = ["ulcerative-colitis-grade-0-1", "ulcerative-colitis-grade-1",
             "ulcerative-colitis-grade-2", "ulcerative-colitis-grade-3"]
ESO_FAMILY = ["barretts", "barretts-short-segment", "esophagitis-a", "esophagitis-b-d"]
SHOTS = [1, 5, 10]
MODES = ["scene", "strict", "none"]      # episode sampling: avoid scene-level / strict near-duplicates / plain random
N_QUERY = 15
METHODS = ["NC-L2N", "SimpleShot-CL2N", "LogReg-CL2N", "NonEpisodic-Head+NC", "ProtoNet-Head (ours)"]


# ----------------------------------------------------------------------------- utilities
class Logger:
    def __init__(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.f = open(path, "a", encoding="utf-8")

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        self.f.write(time.strftime("[%H:%M:%S] ") + s + "\n")
        self.f.flush()


MASK_OVERLAY = False        # set by --mask-overlay: black out on-screen text and the scope-position inset


def out_dir(smoke, masked=None):
    masked = MASK_OVERLAY if masked is None else masked
    d = ROOT / "results" / (("exp1_smoke" if smoke else "exp1") + ("_masked" if masked else ""))
    d.mkdir(parents=True, exist_ok=True)
    return d


def set_seed(s):
    random.seed(s)
    np.random.seed(s)
    try:
        import torch
        torch.manual_seed(s)
    except ImportError:
        pass


def list_images():
    items = []
    for r, _, fs in os.walk(ENDO):
        cls = os.path.basename(r).lower()
        cls = ALIAS.get(cls, cls)
        for f in sorted(fs):
            if f.lower().endswith((".jpg", ".jpeg", ".png")):
                items.append((os.path.join(r, f), cls))
    items.sort()
    return items


# ----------------------------------------------------------------------------- step 1: features
def build_backbone(random_weights=False):
    import timm
    if random_weights:   # container tests only: tiny randomly initialised ViT, same code path otherwise
        return timm.create_model("vit_tiny_patch16_224", pretrained=False, num_classes=0).eval()
    m = timm.create_model("vit_base_patch14_dinov2.lvd142m", pretrained=True, num_classes=0, img_size=224)
    return m.eval()


def load_image(path):
    """Open an endoscopy frame. With --mask-overlay, two corner regions are blacked out: the upper-left text block
    (date/time/settings) and the lower-left scope-position inset (green box) that many HyperKvasir frames carry."""
    from PIL import Image, ImageDraw
    im = Image.open(path).convert("RGB")
    if MASK_OVERLAY:
        w, h = im.size
        d = ImageDraw.Draw(im)
        d.rectangle([0, 0, int(0.28 * w), int(0.35 * h)], fill=(0, 0, 0))
        d.rectangle([0, int(0.60 * h), int(0.33 * w), h], fill=(0, 0, 0))
    return im


def make_transform():
    from torchvision import transforms as T
    return T.Compose([T.Resize(224), T.CenterCrop(224), T.ToTensor(),
                      T.Normalize((0.485, 0.456, 0.406), (0.229, 0.224, 0.225))])


def step_extract(args, log):
    import torch
    from PIL import Image
    od = out_dir(args.smoke)
    fpath = od / "features_dinov2.npz"
    if fpath.exists() and not args.force:
        log(f"[extract] feature file exists, skipping: {fpath}")
        return
    items = list_images()
    if not items:
        sys.exit(f"[extract] no endoscopic images found under {ENDO}")
    if args.smoke:                                    # keep at most 40 images per class
        per = collections.defaultdict(list)
        for p, c in items:
            per[c].append((p, c))
        rng = random.Random(0)
        items = sorted(x for c in per for x in rng.sample(per[c], min(40, len(per[c]))))
    log(f"[extract] {len(items)} images, {len(set(c for _, c in items))} classes")
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model = build_backbone(args.random_backbone).to(dev)
    tf = make_transform()
    feats, md5s = [], []
    bs, t0 = 64, time.time()
    for i in range(0, len(items), bs):
        batch = items[i:i + bs]
        imgs = []
        for p, _ in batch:
            raw = open(p, "rb").read()
            md5s.append(hashlib.md5(raw).hexdigest())
            imgs.append(tf(load_image(p)))
        x = torch.stack(imgs).to(dev)
        with torch.no_grad(), torch.autocast(dev, dtype=torch.float16, enabled=dev == "cuda"):
            feats.append(model(x).float().cpu().numpy())
        if (i // bs) % 20 == 0:
            log(f"[extract] {min(i + bs, len(items))}/{len(items)}  {time.time() - t0:.0f}s")
    feats = np.concatenate(feats).astype(np.float32)
    np.savez_compressed(fpath, feats=feats, labels=np.array([c for _, c in items]),
                        paths=np.array([os.path.relpath(p, ENDO) for p, _ in items]), md5=np.array(md5s))
    log(f"[extract] done: {feats.shape}, {time.time() - t0:.0f}s -> {fpath}")


def load_features(smoke):
    d = np.load(out_dir(smoke) / "features_dinov2.npz", allow_pickle=False)
    return d["feats"], d["labels"].astype(str), d["paths"].astype(str), d["md5"].astype(str)


# ----------------------------------------------------------------------------- step 2: near-duplicates
class UnionFind:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def step_dedup(args, log):
    """Two levels of near-duplicate grouping (same class only; cross-class pairs are reported, never merged):
    strict = the same frame re-saved / almost identical;
    scene  = frames of the same scene from one examination (camera moved a little). HyperKvasir has no patient or
             examination IDs, so scene-level grouping is the conservative proxy used for the main analysis."""
    import imagehash
    from PIL import Image
    od = out_dir(args.smoke)
    feats, labels, paths, md5 = load_features(args.smoke)
    n = len(labels)
    log(f"[dedup] {n} images")
    ph = []
    for i, p in enumerate(paths):
        ph.append(imagehash.phash(Image.open(ENDO / p).convert("RGB"), hash_size=8))
        if i % 2000 == 0:
            log(f"[dedup] perceptual hash {i}/{n}")
    bits = np.array([h.hash.flatten() for h in ph], dtype=np.uint8)          # n x 64
    z = feats / np.linalg.norm(feats, axis=1, keepdims=True)
    T, C, SH, SC = args.phash_thr, args.cos_thr, args.scene_ham, args.scene_cos
    uf_strict, uf_scene = UnionFind(n), UnionFind(n)
    border, cross, n_exact, n_strict, n_scene = [], 0, 0, 0, 0
    for i in range(n):
        j = np.arange(i + 1, n)
        if len(j) == 0:
            break
        same = labels[j] == labels[i]
        ham = (bits[j] != bits[i]).sum(1)
        cos = z[j] @ z[i]
        ex = md5[j] == md5[i]
        strict = ex | ((ham <= T) & (cos >= C - 0.05)) | ((cos >= C) & (ham <= 2 * T))
        scene = strict | ((ham <= SH) & (cos >= SC)) | ((cos >= SC + 0.03) & (ham <= SH + 8))
        n_exact += int(ex.sum()); n_strict += int((strict & same).sum()); n_scene += int((scene & same).sum())
        cross += int((scene & ~same).sum())
        for k in j[strict & same]:
            uf_strict.union(i, int(k))
        for k in j[scene & same]:
            uf_scene.union(i, int(k))
        bmask = same & ~scene & (ham <= SH + 6) & (cos >= SC - 0.05)
        for k in j[bmask][:3]:
            border.append((i, int(k), int((bits[k] != bits[i]).sum()), float(z[k] @ z[i])))
    cl_strict = np.array([uf_strict.find(i) for i in range(n)])
    cl_scene = np.array([uf_scene.find(i) for i in range(n)])
    with open(od / "dedup_clusters.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["path", "label", "cluster", "cluster_scene"])
        for row in zip(paths, labels, cl_strict, cl_scene):
            w.writerow(row)
    lines = [f"strict (re-saved or almost identical frame): identical file; or pHash <= {T} and cosine >= {C - 0.05:.2f}; or cosine >= {C} and pHash <= {2 * T}",
             f"scene (neighboring frames of one examination): strict; or pHash <= {SH} and cosine >= {SC}; or cosine >= {SC + 0.03:.2f} and pHash <= {SH + 8}",
             f"within-class pairs: identical {n_exact}, strict {n_strict}, scene {n_scene}; cross-class scene-level pairs {cross} (reported, not merged)",
             f"{'class':32s} {'images':>6s} {'strict':>8s} {'max':>5s} {'scene':>8s} {'max':>5s}"]
    for c in sorted(set(labels), key=lambda c: (labels == c).sum()):
        m = labels == c
        c1, c2 = collections.Counter(cl_strict[m]), collections.Counter(cl_scene[m])
        warn = "  <-- largest scene cluster exceeds 20% of the class" if max(c2.values()) > 0.2 * m.sum() and m.sum() >= 10 else ""
        lines.append(f"{c:32s} {int(m.sum()):6d} {len(c1):8d} {max(c1.values()):5d} {len(c2):8d} {max(c2.values()):5d}{warn}")
    (od / "dedup_report.txt").write_text("\n".join(lines), encoding="utf-8")
    for s_ in lines:
        log("[dedup] " + s_)
    rng = random.Random(0)
    border = rng.sample(border, min(24, len(border)))
    if border:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        rows = len(border)
        fig, ax = plt.subplots(rows, 2, figsize=(5.2, 2.5 * rows))
        ax = np.atleast_2d(ax)
        for r, (i, k, hd, cs) in enumerate(border):
            for col, idx in enumerate((i, k)):
                ax[r, col].imshow(Image.open(ENDO / paths[idx]).convert("RGB").resize((200, 160)))
                ax[r, col].axis("off")
            ax[r, 0].set_title(f"{r + 1}. {labels[i]}   ham={hd}  cos={cs:.3f}", fontsize=9, loc="left")
        plt.tight_layout()
        fig.savefig(od / "dedup_borderline_pairs.png", dpi=110)
        plt.close(fig)
        log(f"[dedup] borderline pairs (within-class pairs just below the scene threshold): {od / 'dedup_borderline_pairs.png'}")
    log(f"[dedup] done -> {od / 'dedup_clusters.csv'}")


def load_clusters(smoke, paths):
    """Return {'none': None, 'strict': array, 'scene': array}."""
    f = out_dir(smoke) / "dedup_clusters.csv"
    if not f.exists() and MASK_OVERLAY:
        f = out_dir(smoke, masked=False) / "dedup_clusters.csv"   # grouping is a property of the images, reuse it
    m1, m2 = {}, {}
    with open(f, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if "cluster_scene" not in row:
                sys.exit("[run] dedup_clusters.csv is outdated; rerun the dedup step first")
            m1[row["path"]] = int(row["cluster"])
            m2[row["path"]] = int(row["cluster_scene"])
    return {"none": None, "strict": np.array([m1[p] for p in paths]), "scene": np.array([m2[p] for p in paths])}


# ----------------------------------------------------------------------------- splits & episodes
def make_splits(classes, counts, n_random, seed=0):
    ok = sorted(c for c in classes if counts[c] >= MIN_IMAGES)
    splits = {
        "main_rare5": MAIN_TEST,
        "family_UC4": UC_FAMILY,
        "family_esophagus4": ESO_FAMILY,
    }
    rng = random.Random(seed)
    for r in range(n_random):
        splits[f"random5_{r + 1}"] = sorted(rng.sample(ok, 5))
    out = {}
    for name, test in splits.items():
        test = [c for c in test if c in ok]
        train = [c for c in ok if c not in test]
        out[name] = (train, test)
    return out


def _split_one_class(rng, idx, k, q, clusters, dedup, tries=300):
    """Pick k support and q query images of one class. With dedup, no near-duplicate cluster is shared
    between support and query (a cluster touched by the support is never used for the query)."""
    idx = list(idx)
    if not dedup:
        if len(idx) < k + q:
            return None
        rng.shuffle(idx)
        return idx[:k], idx[k:k + q]
    groups = collections.defaultdict(list)
    for i in idx:
        groups[clusters[i]].append(i)
    keys = list(groups)
    for _ in range(tries):
        rng.shuffle(keys)
        s, qq = [], []
        for g in keys:
            members = groups[g][:]
            rng.shuffle(members)
            if len(s) < k:
                s += members[:k - len(s)]
            elif len(qq) < q:
                qq += members[:q - len(qq)]
            else:
                break
        if len(s) == k and len(qq) == q:
            return s, qq
    return None


def sample_episode(rng, idx_by_class, test_classes, k, q, clusters, dedup):
    """Return support / query indices and labels, or None if some class cannot supply them."""
    sup, qry, ys, yq = [], [], [], []
    for ci, c in enumerate(test_classes):
        got = _split_one_class(rng, idx_by_class[c], k, q, clusters, dedup)
        if got is None:
            return None
        s, qq = got
        sup += s; qry += qq; ys += [ci] * len(s); yq += [ci] * len(qq)
    return (np.array(sup, dtype=int), np.array(qry, dtype=int), np.array(ys, dtype=int), np.array(yq, dtype=int))


def episode_budget(idx_by_class, test_classes, k, clusters, dedup):
    """Largest query count per class (<= N_QUERY) that every test class can supply."""
    rng = random.Random(12345)
    for q in range(N_QUERY, 0, -1):
        if all(_split_one_class(rng, idx_by_class[c], k, q, clusters, dedup) is not None for c in test_classes):
            return q
    return 0


# ----------------------------------------------------------------------------- methods
def l2n(x):
    return x / (np.linalg.norm(x, axis=1, keepdims=True) + 1e-8)


def nc_predict(zs, ys, zq, n_way, temp=10.0):
    protos = np.stack([zs[ys == c].mean(0) for c in range(n_way)])
    d = ((zq[:, None, :] - protos[None]) ** 2).sum(-1)
    logit = -d * temp
    logit -= logit.max(1, keepdims=True)
    p = np.exp(logit)
    return p / p.sum(1, keepdims=True)


def logreg_predict(zs, ys, zq, n_way):
    from sklearn.linear_model import LogisticRegression
    clf = LogisticRegression(C=1.0, max_iter=1000)
    clf.fit(zs, ys)
    p = np.zeros((len(zq), n_way))
    p[:, clf.classes_] = clf.predict_proba(zq)
    return p


def head_module(d_in, d_hid=512):
    """Residual adapter on SimpleShot (centred + L2-normalised) features. The last layer starts at zero, so an
    untrained head is exactly SimpleShot; training can only move away from it if the training classes say so."""
    import torch.nn as nn

    class Residual(nn.Module):
        def __init__(self):
            super().__init__()
            self.mlp = nn.Sequential(nn.Linear(d_in, d_hid), nn.ReLU(), nn.Dropout(0.1), nn.Linear(d_hid, d_in))
            nn.init.zeros_(self.mlp[-1].weight)
            nn.init.zeros_(self.mlp[-1].bias)

        def forward(self, x):
            return x + self.mlp(x)

    return Residual()


def split_meta_val(train_classes, seed, n_val=4):
    """Hold out a few TRAINING classes as meta-validation (model selection never sees the test classes)."""
    rng = random.Random(seed)
    tc = sorted(train_classes)
    val = sorted(rng.sample(tc, min(n_val, max(0, len(tc) - 5))))
    return [c for c in tc if c not in val], val


def make_val_episodes(y, val_classes, seed, n=100, shot=5, query=10):
    rng = np.random.default_rng(10_000 + seed)
    idx_by = {c: np.where(y == c)[0] for c in val_classes}
    eps = []
    for _ in range(n):
        s, q = [], []
        for c in val_classes:
            pick = rng.choice(idx_by[c], shot + query, replace=len(idx_by[c]) < shot + query)
            s.append(pick[:shot]); q.append(pick[shot:])
        eps.append((np.concatenate(s), np.concatenate(q)))
    return eps, shot, query


def val_accuracy(head, z, val_eps, dev):
    eps, shot, query = val_eps
    if not eps:
        return 0.0
    e = embed(head, z, dev)
    way = len(eps[0][0]) // shot
    accs = []
    for s, q in eps:
        protos = e[s].reshape(way, shot, -1).mean(1)
        pred = (-((e[q][:, None] - protos[None]) ** 2).sum(-1)).argmax(1)
        accs.append((pred == np.repeat(np.arange(way), query)).mean())
    return float(np.mean(accs))


def train_protonet(x, y, train_classes, seed, episodes, dev, way=10, shot=5, query=10, check_every=250):
    import copy
    import torch
    import torch.nn.functional as F
    set_seed(seed)
    fit_c, val_c = split_meta_val(train_classes, seed)
    val_eps = make_val_episodes(y, val_c, seed) if val_c else ([], 5, 10)
    head = head_module(x.shape[1]).to(dev)
    log_t = torch.nn.Parameter(torch.tensor(math.log(10.0), device=dev))
    opt = torch.optim.Adam(list(head.parameters()) + [log_t], lr=1e-4, weight_decay=1e-4)
    xt = torch.tensor(x, device=dev)
    idx_by = {c: np.where(y == c)[0] for c in fit_c}
    rng = np.random.default_rng(seed)
    way = min(way, len(fit_c))
    head.eval()
    best, best_state = val_accuracy(head, x, val_eps, dev), copy.deepcopy(head.state_dict())   # episode 0 = SimpleShot
    base = best
    for ep in range(1, episodes + 1):
        head.train()
        cls = rng.choice(fit_c, way, replace=False)
        s_idx, q_idx = [], []
        for c in cls:
            pick = rng.choice(idx_by[c], shot + query, replace=len(idx_by[c]) < shot + query)
            s_idx.append(pick[:shot]); q_idx.append(pick[shot:])
        zs = F.normalize(head(xt[np.concatenate(s_idx)]), dim=-1).view(way, shot, -1)
        zq = F.normalize(head(xt[np.concatenate(q_idx)]), dim=-1)
        logits = -torch.cdist(zq, zs.mean(1)) ** 2 * log_t.exp()
        loss = F.cross_entropy(logits, torch.arange(way, device=dev).repeat_interleave(query))
        opt.zero_grad(); loss.backward(); opt.step()
        if ep % check_every == 0 or ep == episodes:
            head.eval()
            acc = val_accuracy(head, x, val_eps, dev)
            if acc > best:
                best, best_state = acc, copy.deepcopy(head.state_dict())
    head.load_state_dict(best_state)
    head.eval()
    return head, best, base


def train_nonepisodic(x, y, train_classes, seed, epochs, dev):
    import copy
    import torch
    import torch.nn.functional as F
    set_seed(seed)
    fit_c, val_c = split_meta_val(train_classes, seed)
    val_eps = make_val_episodes(y, val_c, seed) if val_c else ([], 5, 10)
    cmap = {c: i for i, c in enumerate(fit_c)}
    keep = np.isin(y, fit_c)
    xt = torch.tensor(x[keep], device=dev)
    yt = torch.tensor([cmap[c] for c in y[keep]], device=dev)
    head = head_module(x.shape[1]).to(dev)
    clf = torch.nn.Linear(x.shape[1], len(fit_c)).to(dev)
    opt = torch.optim.Adam(list(head.parameters()) + list(clf.parameters()), lr=1e-4, weight_decay=1e-4)
    counts = torch.bincount(yt, minlength=len(fit_c)).float()
    w_sample = (1.0 / counts)[yt]                 # class-balanced sampling
    n = len(yt)
    head.eval()
    best, best_state = val_accuracy(head, x, val_eps, dev), copy.deepcopy(head.state_dict())
    base = best
    for _ in range(epochs):
        head.train()
        order = torch.multinomial(w_sample, n, replacement=True)
        for b in range(0, n, 256):
            bi = order[b:b + 256]
            logits = clf(F.normalize(head(xt[bi]), dim=-1) * 10.0)
            loss = F.cross_entropy(logits, yt[bi])
            opt.zero_grad(); loss.backward(); opt.step()
        head.eval()
        acc = val_accuracy(head, x, val_eps, dev)
        if acc > best:
            best, best_state = acc, copy.deepcopy(head.state_dict())
    head.load_state_dict(best_state)
    head.eval()
    return head, best, base


def embed(head, x, dev):
    import torch
    import torch.nn.functional as F
    with torch.no_grad():
        return F.normalize(head(torch.tensor(x, device=dev)), dim=-1).cpu().numpy()


def metrics(p, yq, n_way):
    from sklearn.metrics import f1_score, roc_auc_score
    pred = p.argmax(1)
    acc = float((pred == yq).mean())
    f1 = float(f1_score(yq, pred, average="macro", labels=list(range(n_way)), zero_division=0))
    try:
        auc = float(roc_auc_score(yq, p, multi_class="ovr", average="macro", labels=list(range(n_way))))
    except ValueError:
        auc = float("nan")
    return acc, f1, auc


# ----------------------------------------------------------------------------- step 3: episodic evaluation
def step_run(args, log):
    import torch
    od = out_dir(args.smoke)
    feats, labels, paths, _ = load_features(args.smoke)
    clmodes = load_clusters(args.smoke, paths)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    counts = collections.Counter(labels)
    splits = make_splits(sorted(counts), counts, n_random=1 if args.smoke else 5)
    if args.splits:
        splits = {k: v for k, v in splits.items() if k in args.splits.split(",")}
    n_ep = 20 if args.smoke else args.episodes
    seeds = [0] if args.smoke else [0, 1, 2]
    train_eps = 100 if args.smoke else 3000
    epochs = 2 if args.smoke else 30
    (od / "splits.json").write_text(json.dumps({k: {"train": v[0], "test": v[1]} for k, v in splits.items()},
                                               indent=1, ensure_ascii=False), encoding="utf-8")
    rows, qrows, trainlog = [], [], []
    t0 = time.time()
    for sname, (train_c, test_c) in splits.items():
        if len(test_c) < 2:
            log(f"[run] {sname}: not enough novel classes, skipped"); continue
        n_way = len(test_c)
        idx_by = {c: np.where(labels == c)[0] for c in test_c}
        base = feats[np.isin(labels, train_c)]
        mu = base.mean(0, keepdims=True)                      # SimpleShot centering uses meta-train classes only
        z_l2 = l2n(feats)
        z_cl2 = l2n(feats - mu)
        heads = {}
        for seed in seeds:
            for kind, fn, arg in (("proto", train_protonet, train_eps), ("nonep", train_nonepisodic, epochs)):
                h, best, base = fn(z_cl2, labels, train_c, seed, arg, dev)
                heads[(kind, seed)] = h
                trainlog.append([sname, kind, seed, round(base * 100, 2), round(best * 100, 2)])
        emb = {k: embed(h, z_cl2, dev) for k, h in heads.items()}
        log(f"[run] {sname}: {n_way}-way, {len(train_c)} base classes, adapters trained in {time.time() - t0:.0f}s")
        for k in SHOTS:
            # one query count for all three sampling modes, so the modes differ only in near-duplicate control
            budgets = {md: episode_budget(idx_by, test_c, k, clmodes[md], md != "none") for md in MODES}
            # a sampling mode that can only supply a handful of query images (e.g. a class whose frames collapse
            # into one or two scene clusters) is skipped instead of dragging every mode down to that query count
            feasible = [b_ for b_ in budgets.values() if b_ >= 5] or [b_ for b_ in budgets.values() if b_ >= 1]
            q = min(feasible) if feasible else 0
            for dedup in MODES:
                clusters = clmodes[dedup]
                if q < 1 or budgets[dedup] < q:
                    log(f"[run] {sname} {k}-shot {dedup}: too few images for this sampling mode (at most {budgets[dedup]} queries per class), skipped")
                    continue
                rng = random.Random(1000 + k)               # identical episodes for every method
                skipped = 0
                for e in range(n_ep):
                    got = sample_episode(rng, idx_by, test_c, k, q, clusters, dedup != "none")
                    if got is None:
                        skipped += 1
                        continue
                    s, qi, ys, yq = got
                    preds = {
                        "NC-L2N": [nc_predict(z_l2[s], ys, z_l2[qi], n_way)],
                        "SimpleShot-CL2N": [nc_predict(z_cl2[s], ys, z_cl2[qi], n_way)],
                        "LogReg-CL2N": [logreg_predict(z_cl2[s], ys, z_cl2[qi], n_way)],
                        "NonEpisodic-Head+NC": [nc_predict(emb[("nonep", sd)][s], ys, emb[("nonep", sd)][qi], n_way)
                                                for sd in seeds],
                        "ProtoNet-Head (ours)": [nc_predict(emb[("proto", sd)][s], ys, emb[("proto", sd)][qi], n_way)
                                                 for sd in seeds],
                    }
                    for mname, plist in preds.items():
                        trained = mname in ("NonEpisodic-Head+NC", "ProtoNet-Head (ours)")
                        for sd, p in zip(seeds if trained else [-1], plist):
                            acc, f1, auc = metrics(p, yq, n_way)
                            rows.append([sname, dedup, k, q, e, mname, sd, acc, f1, auc])
                            if sname == "main_rare5":
                                pred = p.argmax(1)
                                for ii, yy, pp in zip(qi, yq, pred):
                                    qrows.append([dedup, k, mname, sd, e, int(ii), int(clmodes["scene"][ii]), int(pp == yy)])
                log(f"[run] {sname} {dedup} {k}-shot ({q} queries per class) done"
                    + (f", {skipped} episodes could not be sampled" if skipped else "") + f", elapsed {time.time() - t0:.0f}s")
    with open(od / "episodes.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["split", "sampling", "shot", "query_per_class", "episode", "method", "seed", "acc", "f1", "auc"])
        w.writerows(rows)
    with open(od / "head_training.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["split", "head", "seed", "metaval_acc_untrained_SimpleShot", "metaval_acc_selected"])
        w.writerows(trainlog)
    with open(od / "main_query_predictions.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["sampling", "shot", "method", "seed", "episode", "image_index", "cluster", "correct"])
        w.writerows(qrows)
    log(f"[run] done: {len(rows)} rows -> {od / 'episodes.csv'}, {time.time() - t0:.0f}s")


# ----------------------------------------------------------------------------- step 4: tables & figures
def cluster_bootstrap(df_q, n_boot=1000, seed=0):
    """95% CI of accuracy resampling near-duplicate clusters (each query image's correctness pooled per cluster)."""
    g = df_q.groupby("cluster")["correct"].agg(["sum", "count"]).to_numpy()
    rng = np.random.default_rng(seed)
    stats = []
    for _ in range(n_boot):
        pick = g[rng.integers(0, len(g), len(g))]
        stats.append(pick[:, 0].sum() / pick[:, 1].sum())
    return np.percentile(stats, [2.5, 97.5])


def step_figures(args, log):
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    od = out_dir(args.smoke)
    df = pd.read_csv(od / "episodes.csv")
    if df.empty:
        log("[figures] episodes.csv is empty (no episode could be sampled); check log.txt"); return
    # seed-average first, so every episode counts once per method
    ep = df.groupby(["split", "sampling", "shot", "method", "episode"])[["acc", "f1", "auc"]].mean().reset_index()
    out = []
    for (sp, dd, k, m), g in ep.groupby(["split", "sampling", "shot", "method"]):
        r = {"split": sp, "sampling": dd, "shot": k, "method": m, "episodes": len(g)}
        for col in ["acc", "f1", "auc"]:
            v = g[col].dropna().to_numpy() * 100
            r[col] = v.mean()
            r[col + "_ci95"] = 1.96 * v.std(ddof=1) / math.sqrt(len(v)) if len(v) > 1 else float("nan")
        out.append(r)
    summ = pd.DataFrame(out)
    # spread across the random class splits (between-split variability)
    rnd = summ[summ.split.str.startswith("random5_")]
    if len(rnd):
        between = rnd.groupby(["sampling", "shot", "method"])["acc"].agg(["mean", "std", "count"]).reset_index()
        between.to_csv(od / "summary_random_splits.csv", index=False)
    # cluster bootstrap CI for the main split
    qp = od / "main_query_predictions.csv"
    if qp.exists():
        qdf = pd.read_csv(qp)
        cb = []
        for (dd, k, m), g in qdf.groupby(["sampling", "shot", "method"]):
            lo, hi = cluster_bootstrap(g, n_boot=200 if args.smoke else 1000)
            cb.append({"split": "main_rare5", "sampling": dd, "shot": k, "method": m,
                       "acc_cluster_ci_lo": lo * 100, "acc_cluster_ci_hi": hi * 100})
        summ = summ.merge(pd.DataFrame(cb), on=["split", "sampling", "shot", "method"], how="left")
    summ.to_csv(od / "summary.csv", index=False)
    wide = ep.pivot_table(index=["split", "sampling", "shot", "episode"], columns="method", values="acc")
    ours = "ProtoNet-Head (ours)"
    pd_rows = []
    for (sp, dd, k), g in wide.groupby(level=[0, 1, 2]):
        for m in METHODS:
            if m == ours or m not in g or ours not in g:
                continue
            d = (g[ours] - g[m]).dropna().to_numpy() * 100
            if len(d) > 1:
                h = 1.96 * d.std(ddof=1) / math.sqrt(len(d))
                pd_rows.append({"split": sp, "sampling": dd, "shot": k, "versus": m, "diff_pp": d.mean(),
                                "ci_lo": d.mean() - h, "ci_hi": d.mean() + h, "episodes": len(d),
                                "ours_better_share": float((d > 0).mean())})
    paired = pd.DataFrame(pd_rows)
    paired.to_csv(od / "paired_differences.csv", index=False)
    # readable markdown table: main split, dedup on
    lines = ["# Experiment 1 summary (main split = 5 rare lesion classes; episodes avoid scene-level near-duplicates)", "",
             "| Method | Shot | Accuracy % (95% CI) | Macro-F1 % | Macro-AUC % | Cluster-bootstrap CI |",
             "|---|---|---|---|---|---|"]
    main_mode = next((md for md in MODES if ((summ.split == "main_rare5") & (summ.sampling == md)).any()), "none")
    mm = summ[(summ.split == "main_rare5") & (summ.sampling == main_mode)]
    lines[0] += f" — sampling used in this table: {main_mode}"
    for k in SHOTS:
        for m in METHODS:
            r = mm[(mm.shot == k) & (mm.method == m)]
            if len(r):
                r = r.iloc[0]
                cbs = f"[{r.get('acc_cluster_ci_lo', float('nan')):.1f}, {r.get('acc_cluster_ci_hi', float('nan')):.1f}]"
                lines.append(f"| {m} | {k} | {r.acc:.1f} ± {r.acc_ci95:.1f} | {r.f1:.1f} | {r.auc:.1f} | {cbs} |")
    pm = paired[(paired.split == "main_rare5") & (paired.sampling == main_mode)] if len(paired) else paired
    if len(pm):
        lines += ["", "Paired difference, ours minus each method (percentage points, same episodes, 95% CI):", "",
                  "| Shot | Versus | Difference | 95% CI | Episodes where ours is better |", "|---|---|---|---|---|"]
        for _, r in pm.sort_values(["shot", "versus"]).iterrows():
            lines.append(f"| {r.shot} | {r.versus} | {r.diff_pp:+.1f} | [{r.ci_lo:+.1f}, {r.ci_hi:+.1f}] | {r.ours_better_share * 100:.0f}% |")
    lines += ["", "Full tables: summary.csv (all splits and sampling modes), paired_differences.csv"]
    (od / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    log("[figures] " + "\n".join(lines))
    # accuracy vs shot, main split
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), sharey=True)
    titles = {"scene": "Scene-level near-duplicates separated", "strict": "Strict near-duplicates separated",
              "none": "Plain random episodes"}
    for ax, dd in zip(axes, ["scene", "strict", "none"]):
        s = summ[(summ.split == "main_rare5") & (summ.sampling == dd)]
        for m in METHODS:
            g = s[s.method == m].sort_values("shot")
            if len(g):
                ax.errorbar(g.shot, g.acc, yerr=g.acc_ci95, marker="o", capsize=3, label=m)
        ax.set_xticks(SHOTS); ax.set_xlabel("Support images per class")
        ax.set_title(titles[dd], fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Accuracy (%)")
    axes[2].legend(fontsize=7, loc="lower right")
    plt.tight_layout()
    fig.savefig(od / "fig_accuracy_vs_shot.png", dpi=300)
    plt.close(fig)
    log(f"[figures] figure: {od / 'fig_accuracy_vs_shot.png'}")
    step_occlusion(args, log)


def step_occlusion(args, log):
    """Occlusion attribution of the 'query -> predicted prototype' score of our head (main split, 5-shot)."""
    import torch
    import torch.nn.functional as F
    from PIL import Image
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    od = out_dir(args.smoke)
    feats, labels, paths, _ = load_features(args.smoke)
    clusters = load_clusters(args.smoke, paths)["scene"]
    counts = collections.Counter(labels)
    train_c, test_c = make_splits(sorted(counts), counts, 0)["main_rare5"]
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    mu = feats[np.isin(labels, train_c)].mean(0, keepdims=True)
    z_cl2 = l2n(feats - mu)
    head = train_protonet(z_cl2, labels, train_c, 0, 100 if args.smoke else 3000, dev)[0]
    backbone = build_backbone(args.random_backbone).to(dev)
    tf = make_transform()
    idx_by = {c: np.where(labels == c)[0] for c in test_c}
    k = 5
    dedup = episode_budget(idx_by, test_c, k, clusters, True) >= 1
    q = episode_budget(idx_by, test_c, k, clusters, dedup)
    got = sample_episode(random.Random(7), idx_by, test_c, k, min(q, 15), clusters, dedup) if q >= 1 else None
    if got is None:
        log("[occlusion] no episode could be sampled, skipping occlusion maps"); return
    s, qi, ys, yq = got
    zs = embed(head, z_cl2[s], dev)
    protos = np.stack([zs[ys == c].mean(0) for c in range(len(test_c))])
    zq = embed(head, z_cl2[qi], dev)
    pred = (-((zq[:, None] - protos[None]) ** 2).sum(-1)).argmax(1)
    correct = [i for i in range(len(qi)) if pred[i] == yq[i]][:3]
    wrong = [i for i in range(len(qi)) if pred[i] != yq[i]][:3]
    chosen = [(i, True) for i in correct] + [(i, False) for i in wrong]
    if not chosen:
        log("[occlusion] no cases to plot"); return
    grid, patch = 14, 16
    fig, ax = plt.subplots(len(chosen), 2, figsize=(5, 2.4 * len(chosen)))
    ax = np.atleast_2d(ax)
    pt = torch.tensor(protos, device=dev)
    for r, (i, ok) in enumerate(chosen):
        img = tf(load_image(ENDO / paths[qi[i]]))
        batch = [img]
        for gy in range(grid):
            for gx in range(grid):
                im = img.clone()
                im[:, gy * patch:(gy + 1) * patch, gx * patch:(gx + 1) * patch] = 0.0
                batch.append(im)
        x = torch.stack(batch).to(dev)
        with torch.no_grad(), torch.autocast(dev, dtype=torch.float16, enabled=dev == "cuda"):
            f = backbone(x).float()
        with torch.no_grad():
            f = F.normalize(f - torch.tensor(mu, device=dev, dtype=f.dtype), dim=-1)
            z = F.normalize(head(f), dim=-1)
            score = -((z - pt[pred[i]]) ** 2).sum(-1)          # score towards the predicted prototype
        drop = (score[0] - score[1:]).cpu().numpy().reshape(grid, grid)
        show = img.permute(1, 2, 0).numpy() * np.array([0.229, 0.224, 0.225]) + np.array([0.485, 0.456, 0.406])
        ax[r, 0].imshow(np.clip(show, 0, 1)); ax[r, 0].axis("off")
        ax[r, 0].set_title(f"true: {test_c[yq[i]]}\npred: {test_c[pred[i]]} ({'correct' if ok else 'wrong'})", fontsize=6)
        ax[r, 1].imshow(np.clip(show, 0, 1))
        lim = float(np.abs(drop).max()) or 1.0
        hm = ax[r, 1].imshow(np.kron(drop, np.ones((patch, patch))), cmap="RdBu_r", vmin=-lim, vmax=lim, alpha=0.55)
        ax[r, 1].axis("off")
        ax[r, 1].set_title("occlusion (red = supports prediction)", fontsize=6)
        fig.colorbar(hm, ax=ax[r, 1], fraction=0.046, pad=0.02)
    plt.tight_layout()
    fig.savefig(od / "fig_occlusion.png", dpi=300)
    plt.close(fig)
    log(f"[occlusion] figure: {od / 'fig_occlusion.png'} (model visualization, not a clinical explanation)")


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["extract", "dedup", "run", "figures", "all"])
    ap.add_argument("--smoke", action="store_true", help="few-minute rehearsal on a small subset")
    ap.add_argument("--force", action="store_true", help="recompute features even if they exist")
    ap.add_argument("--episodes", type=int, default=600)
    ap.add_argument("--phash-thr", type=int, default=6)
    ap.add_argument("--cos-thr", type=float, default=0.97)
    ap.add_argument("--scene-ham", type=int, default=16)
    ap.add_argument("--scene-cos", type=float, default=0.92)
    ap.add_argument("--mask-overlay", action="store_true", help="sensitivity analysis: black out text and inset")
    ap.add_argument("--splits", default="", help="comma-separated subset of splits, e.g. main_rare5")
    ap.add_argument("--random-backbone", action="store_true", help=argparse.SUPPRESS)  # container tests only
    args = ap.parse_args()
    global MASK_OVERLAY
    MASK_OVERLAY = args.mask_overlay
    log = Logger(out_dir(args.smoke) / "log.txt")
    log(f"===== exp1 {args.step} smoke={args.smoke} root={ROOT} =====")
    steps = ["extract", "dedup", "run", "figures"] if args.step == "all" else [args.step]
    if MASK_OVERLAY and "dedup" in steps:
        # near-duplicate groups describe the images themselves, so the masked run reuses the unmasked grouping
        if not (out_dir(args.smoke, masked=False) / "dedup_clusters.csv").exists():
            sys.exit("[mask] run the unmasked analysis first (its dedup_clusters.csv is reused)")
        steps.remove("dedup")
        log("[mask] reusing near-duplicate clusters of the unmasked run")
    for s in steps:
        {"extract": step_extract, "dedup": step_dedup, "run": step_run, "figures": step_figures}[s](args, log)
    log("===== all done =====")


if __name__ == "__main__":
    main()
