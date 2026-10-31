"""Study 2, step 1 - frozen tile features for every Kather TCGA tile (read straight from the zip files).

Encoders: uni (main; no TCGA in pre-training), resnet50 (ImageNet, CLAM-style truncated after layer3, 1024-d),
          phikon (sensitivity analysis; pre-trained on TCGA - must be disclosed).
Output  : <GIFS_ROOT>/results/exp2/features/<COHORT>_<encoder>.h5 with datasets
          feats (N x D, float16), patient (N), tile (N, path inside the zip), source_zip (N)
Each (cohort, encoder) file is written once and skipped on re-runs; an interrupted file is redone from scratch.

Usage: python exp2_extract.py [--encoders uni,resnet50,phikon] [--cohorts STAD,CRC] [--smoke]
"""
import argparse, io, os, re, sys, time, zipfile
from pathlib import Path

import numpy as np

ROOT = Path(os.environ.get("GIFS_ROOT", str(Path.home() / "gifs_data")))
TILES = ROOT / "data" / "wsi_tiles"
PID = re.compile(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})")
ZIPS = {"STAD": ["STAD_TRAIN_MSS.zip", "STAD_TRAIN_MSIMUT.zip", "STAD_TEST_MSS.zip", "STAD_TEST_MSIMUT.zip"],
        "CRC": ["CRC_DX_TRAIN_MSS.zip", "CRC_DX_TRAIN_MSIMUT.zip", "CRC_DX_TEST_MSS.zip", "CRC_DX_TEST_MSIMUT.zip"]}
MEAN, STD = (0.485, 0.456, 0.406), (0.229, 0.224, 0.225)     # ImageNet statistics, used by all three encoders


def log(msg, f=[None]):
    if f[0] is None:
        d = ROOT / "results" / "exp2"
        d.mkdir(parents=True, exist_ok=True)
        f[0] = open(d / "extract_log.txt", "a", encoding="utf-8")
    line = time.strftime("[%H:%M:%S] ") + msg
    print(line, flush=True)
    f[0].write(line + "\n"); f[0].flush()


def list_tiles(cohort, smoke):
    items = []
    for z in ZIPS[cohort]:
        zp = TILES / cohort / z
        if not zp.exists():
            sys.exit(f"missing {zp}")
        names = sorted(n for n in zipfile.ZipFile(zp).namelist() if n.lower().endswith(".png"))
        if smoke:
            by_p = {}
            for n in names:
                by_p.setdefault(PID.search(n).group(1), []).append(n)
            names = [n for p in sorted(by_p)[:6] for n in by_p[p][:20]]     # 6 patients x 20 tiles per zip
        items += [(z, n, PID.search(n).group(1)) for n in names]
    return items


class TileSet:
    """torch Dataset over tiles inside zip files; each worker opens its own zip handles."""

    def __init__(self, cohort, items):
        self.cohort, self.items, self.handles = cohort, items, {}
        import torch
        self.mean = torch.tensor(MEAN).view(3, 1, 1)
        self.std = torch.tensor(STD).view(3, 1, 1)

    def __len__(self):
        return len(self.items)

    def __getitem__(self, i):
        import torch
        from PIL import Image
        z, name, _ = self.items[i]
        if z not in self.handles:
            self.handles[z] = zipfile.ZipFile(TILES / self.cohort / z)
        im = Image.open(io.BytesIO(self.handles[z].read(name))).convert("RGB")
        if im.size != (224, 224):
            im = im.resize((224, 224))
        x = torch.from_numpy(np.asarray(im, dtype=np.uint8).copy()).permute(2, 0, 1).float() / 255.0
        return (x - self.mean) / self.std


def read_token():
    f = Path(__file__).resolve().parent / "hf_token.txt"
    return f.read_text(encoding="utf-8-sig").strip() if f.exists() else os.environ.get("HF_TOKEN", "").strip()


def build_encoder(name, random_weights=False):
    """Return (model, output_dim). random_weights is for container tests only."""
    import torch
    import timm
    if name == "uni":
        m = timm.create_model("vit_large_patch16_224", img_size=224, patch_size=16, init_values=1e-5,
                              num_classes=0, dynamic_img_size=True)
        if not random_weights:
            from huggingface_hub import hf_hub_download
            tok = read_token()
            if not tok:
                sys.exit("UNI requires a Hugging Face access token: set HF_TOKEN or put it in hf_token.txt next to this script")
            path = hf_hub_download("MahmoodLab/UNI", "pytorch_model.bin", endpoint="https://huggingface.co", token=tok)
            m.load_state_dict(torch.load(path, map_location="cpu"), strict=True)
        return m.eval(), 1024
    if name == "resnet50":
        import torchvision
        r = torchvision.models.resnet50(weights=None if random_weights else "IMAGENET1K_V2")
        trunk = torch.nn.Sequential(r.conv1, r.bn1, r.relu, r.maxpool, r.layer1, r.layer2, r.layer3,
                                    torch.nn.AdaptiveAvgPool2d(1), torch.nn.Flatten())
        return trunk.eval(), 1024
    if name == "phikon":
        from transformers import ViTConfig, ViTModel
        if random_weights:
            vit = ViTModel(ViTConfig(), add_pooling_layer=False)
        else:
            vit = ViTModel.from_pretrained("owkin/phikon", add_pooling_layer=False)

        class CLS(torch.nn.Module):
            def __init__(self, v):
                super().__init__()
                self.v = v

            def forward(self, x):
                return self.v(pixel_values=x).last_hidden_state[:, 0]

        return CLS(vit).eval(), 768
    raise ValueError(name)


def extract(cohort, enc_name, args):
    import h5py
    import torch
    out = ROOT / "results" / ("exp2_smoke" if args.smoke else "exp2") / "features"
    out.mkdir(parents=True, exist_ok=True)
    final, tmp = out / f"{cohort}_{enc_name}.h5", out / f"{cohort}_{enc_name}.h5.part"
    if final.exists():
        log(f"[{cohort} {enc_name}] exists, skipping: {final}")
        return
    items = list_tiles(cohort, args.smoke)
    if args.max_tiles:
        items = items[:args.max_tiles]
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model, dim = build_encoder(enc_name, args.random_backbone)
    model = model.to(dev)
    from torch.utils.data import DataLoader
    dl = DataLoader(TileSet(cohort, items), batch_size=args.batch, num_workers=args.workers,
                    persistent_workers=args.workers > 0, pin_memory=dev == "cuda")
    n = len(items)
    log(f"[{cohort} {enc_name}] {n} tiles, {len(set(p for _, _, p in items))} patients, feature dimension {dim}")
    if tmp.exists():
        tmp.unlink()
    t0, done = time.time(), 0
    with h5py.File(tmp, "w") as h:
        feats = h.create_dataset("feats", (n, dim), dtype="float16", chunks=(min(4096, n), dim))
        str_dt = h5py.string_dtype()
        h.create_dataset("patient", data=np.array([p for _, _, p in items], dtype=object), dtype=str_dt)
        h.create_dataset("tile", data=np.array([t for _, t, _ in items], dtype=object), dtype=str_dt)
        h.create_dataset("source_zip", data=np.array([z for z, _, _ in items], dtype=object), dtype=str_dt)
        h.attrs["encoder"], h.attrs["cohort"] = enc_name, cohort
        with torch.no_grad():
            for x in dl:
                x = x.to(dev, non_blocking=True)
                with torch.autocast(dev, dtype=torch.float16, enabled=dev == "cuda"):
                    f = model(x).float()
                feats[done:done + len(f)] = f.cpu().numpy().astype(np.float16)
                done += len(f)
                if (done // args.batch) % 200 == 0 or done == n:
                    el = time.time() - t0
                    log(f"[{cohort} {enc_name}] {done}/{n}  {done / max(el, 1e-6):.0f} tiles/s  "
                        f"about {(n - done) / max(done / max(el, 1e-6), 1e-6) / 60:.0f} min left")
        bad = int((~np.isfinite(feats[:].astype(np.float32))).any(1).sum())
        h.attrs["n_nonfinite_rows"] = bad
    tmp.replace(final)
    log(f"[{cohort} {enc_name}] done in {(time.time() - t0) / 60:.1f} min, non-finite rows {bad} -> {final}")
    del model
    if dev == "cuda":
        torch.cuda.empty_cache()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoders", default="uni,resnet50,phikon")
    ap.add_argument("--cohorts", default="STAD,CRC")
    ap.add_argument("--batch", type=int, default=128)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--random-backbone", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--max-tiles", type=int, default=0, help=argparse.SUPPRESS)
    args = ap.parse_args()
    log(f"===== exp2 extract encoders={args.encoders} cohorts={args.cohorts} smoke={args.smoke} =====")
    for enc in args.encoders.split(","):
        for c in args.cohorts.split(","):
            extract(c, enc, args)
    log("===== feature extraction done =====")


if __name__ == "__main__":
    main()
