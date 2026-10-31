"""Print the main result tables from the files in results/ (no GPU or raw data needed)."""
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
S1 = REPO / "results" / "study1"
S2 = REPO / "results" / "study2"
SETTINGS = [("uni", "cap500", "main"), ("resnet50", "cap500", "sensitivity"), ("phikon", "cap500", "sensitivity"),
            ("uni", "cap100", "sensitivity"), ("uni", "cap0", "sensitivity")]


def study1_main():
    s = pd.read_csv(S1 / "summary.csv")
    m = s[(s.split == "main_rare5") & (s.sampling == "scene")]
    t = m.pivot_table(index="method", columns="shot", values="acc").round(1)
    print("\nStudy 1, rare-lesion split, scene-level control: accuracy (%) by number of support images")
    print(t.to_string())
    return t


def study1_robustness():
    p = pd.read_csv(S1 / "paired_differences.csv")
    r = p[(p.split == "main_rare5") & (p.shot == 5)][["sampling", "versus", "diff_pp", "ci_lo", "ci_hi"]].round(1)
    print("\nStudy 1, 5-shot paired gain of the episodic adapter (percentage points)")
    print(r.to_string(index=False))
    return r


def study2_main():
    out = {}
    for c in ["STAD", "CRC"]:
        s = pd.read_csv(S2 / "main" / f"v2_{c}_uni_cap500_summary.csv").set_index("method")
        out[c] = s["pooled_c"].round(3)
    t = pd.DataFrame(out)
    print("\nStudy 2, pooled out-of-fold C-index (UNI features, at most 500 tiles per patient)")
    print(t.to_string())
    return t


def study2_contrasts():
    rows = []
    for enc, cap, sub in SETTINGS:
        for c in ["STAD", "CRC"]:
            k = pd.read_csv(S2 / sub / f"v2_{c}_{enc}_{cap}_contrasts.csv")
            k["contrast"] = k.contrast.str.replace("\u2212", "-")
            for name in ["Cox+DCA - Clinical-Cox", "Cox+WSI - Clinical-Cox"]:
                r = k.set_index("contrast").loc[name]
                rows.append({"setting": f"{enc} {cap}", "cohort": c, "contrast": name,
                             "diff": round(r.diff_c, 3), "ci": f"({r.ci_lo:.3f}, {r.ci_hi:.3f})"})
    t = pd.DataFrame(rows)
    print("\nStudy 2, differences in pooled C-index versus the clinical Cox model (95% bootstrap CI)")
    print(t.to_string(index=False))
    return t


def main():
    study1_main(); study1_robustness(); study2_main(); study2_contrasts()


if __name__ == "__main__":
    main()
