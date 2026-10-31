"""Per-class 5-shot accuracy and paired gains on the rare-lesion split (scene-level near-duplicate control)."""
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
R = REPO / "results" / "study1"


def main():
    d = pd.read_csv(R / "dedup_clusters.csv")
    q = pd.read_csv(R / "main_query_predictions.csv.gz")
    q = q[(q.shot == 5) & (q.sampling == "scene")].copy()
    q["label"] = d.label.values[q.image_index.values]
    e = q.groupby(["method", "episode", "label"]).correct.mean().unstack("method")
    rows = []
    for lab, g in e.groupby(level="label"):
        r = {"label": lab}
        for m in g.columns:
            r[m] = 100 * g[m].mean()
        for vs in ["SimpleShot-CL2N", "NonEpisodic-Head+NC"]:
            diff = 100 * (g["ProtoNet-Head (ours)"] - g[vs])
            r["d_" + vs] = diff.mean()
            r["ci_" + vs] = 1.96 * diff.std(ddof=1) / np.sqrt(len(diff))
        rows.append(r)
    out = pd.DataFrame(rows).set_index("label")
    out.to_csv(R / "per_class_5shot.csv")
    print(out.round(1).T.to_string())


if __name__ == "__main__":
    main()
