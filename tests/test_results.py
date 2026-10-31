"""Consistency checks on the shipped result files."""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import make_tables  # noqa: E402


def test_study1_main_accuracy_table():
    t = make_tables.study1_main()
    assert set(t.columns) == {1, 5, 10}
    assert t.loc["ProtoNet-Head (ours)", 5] > t.loc["SimpleShot-CL2N", 5]


def test_clinical_cox_identical_across_settings():
    vals = {}
    for enc, cap, sub in make_tables.SETTINGS:
        for c in ["STAD", "CRC"]:
            s = pd.read_csv(REPO / "results" / "study2" / sub / f"v2_{c}_{enc}_{cap}_summary.csv").set_index("method")
            vals.setdefault(c, set()).add(round(float(s.loc["Clinical-Cox", "pooled_c"]), 6))
    assert all(len(v) == 1 for v in vals.values())


def test_zero_residual_matches_clinical_ranking():
    for c in ["STAD", "CRC"]:
        s = pd.read_csv(REPO / "results" / "study2" / "main" / f"v2_{c}_uni_cap500_summary.csv").set_index("method")
        assert abs(s.loc["Cox+zero-residual", "pooled_c"] - s.loc["Clinical-Cox", "pooled_c"]) < 1e-3


def test_contrast_table_complete():
    t = make_tables.study2_contrasts()
    assert len(t) == 5 * 2 * 2
