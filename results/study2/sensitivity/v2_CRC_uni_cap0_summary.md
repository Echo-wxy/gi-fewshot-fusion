# Experiment 2, protocol v2 — CRC (encoder uni, tiles <= 0)
Patients 358, deaths 74, censoring 79%; outer 5-fold x 3 seeds; inner 3-fold ensembles.

| Method | Harrell C (mean ± SD) | Uno C | IBS | Brier 12 m | Brier 24 m | Pooled C |
|---|---|---|---|---|---|---|
| Clinical-Cox | 0.681 ± 0.067 | 0.702 | 0.111 | 0.091 | 0.133 | 0.687 |
| Cox+zero-residual | 0.681 ± 0.067 | 0.702 | 0.113 | 0.094 | 0.132 | 0.687 |
| Omics-CoxNet | 0.496 ± 0.090 | 0.512 | 0.131 | 0.098 | 0.152 | 0.493 |
| Stacked-Cox | 0.698 ± 0.065 | 0.723 | 0.110 | 0.090 | 0.132 | 0.703 |
| Late-Fusion | 0.689 ± 0.100 | 0.715 | — | — | — | 0.701 |
| WSI-only | 0.653 ± 0.078 | 0.682 | 0.120 | 0.097 | 0.137 | 0.672 |
| Omics-only | 0.533 ± 0.107 | 0.535 | 0.122 | 0.098 | 0.142 | 0.535 |
| Concat-Fusion | 0.680 ± 0.093 | 0.701 | 0.118 | 0.096 | 0.135 | 0.700 |
| Attention-Fusion | 0.650 ± 0.086 | 0.667 | 0.119 | 0.096 | 0.136 | 0.674 |
| Gated-Fusion | 0.676 ± 0.100 | 0.699 | 0.119 | 0.096 | 0.135 | 0.706 |
| DCA | 0.657 ± 0.100 | 0.675 | 0.119 | 0.096 | 0.136 | 0.684 |
| PORPOISE-style | 0.613 ± 0.078 | 0.633 | 0.121 | 0.096 | 0.139 | 0.619 |
| MCAT-style | 0.634 ± 0.087 | 0.659 | 0.122 | 0.096 | 0.141 | 0.646 |
| Cox+WSI | 0.717 ± 0.085 | 0.746 | 0.113 | 0.093 | 0.129 | 0.734 |
| Cox+Omics | 0.647 ± 0.085 | 0.666 | 0.116 | 0.095 | 0.136 | 0.656 |
| Cox+Concat | 0.706 ± 0.076 | 0.731 | 0.113 | 0.094 | 0.132 | 0.717 |
| Cox+Attention | 0.689 ± 0.076 | 0.710 | 0.115 | 0.094 | 0.132 | 0.704 |
| Cox+Gated | 0.696 ± 0.083 | 0.722 | 0.115 | 0.095 | 0.132 | 0.710 |
| Cox+DCA | 0.692 ± 0.082 | 0.719 | 0.115 | 0.094 | 0.132 | 0.707 |

Primary contrast (pre-specified): Cox+DCA − Clinical-Cox, pooled out-of-fold C, patient-level bootstrap.

| Contrast | ΔC | 95% CI | Share of bootstrap ≤ 0 |
|---|---|---|---|
| Cox+DCA − Clinical-Cox | +0.020 | [-0.031, +0.074] | 0.231 |
| Cox+DCA − Cox+Attention | +0.002 | [-0.015, +0.020] | 0.383 |
| Cox+DCA − Cox+Concat | -0.011 | [-0.034, +0.014] | 0.807 |
| Cox+DCA − Cox+Gated | -0.003 | [-0.012, +0.005] | 0.775 |
| Cox+DCA − Stacked-Cox | +0.004 | [-0.033, +0.045] | 0.435 |
| Cox+DCA − Late-Fusion | +0.005 | [-0.030, +0.039] | 0.392 |
| Cox+WSI − Clinical-Cox | +0.047 | [-0.007, +0.107] | 0.046 |
| Cox+Omics − Clinical-Cox | -0.031 | [-0.069, +0.005] | 0.954 |
| DCA − Attention-Fusion | +0.010 | [-0.013, +0.032] | 0.184 |
| DCA − Concat-Fusion | -0.015 | [-0.047, +0.017] | 0.831 |
| DCA − Gated-Fusion | -0.022 | [-0.039, -0.006] | 0.998 |

Training diagnostics (median over all inner models):

| Model | Best epoch | Stop epoch | Best val NLL | Last train NLL | Grad norm | Seconds |
|---|---|---|---|---|---|---|
| Attention-Fusion | 9 | 29 | 0.621 | 0.428 | 0.75 | 0.6 |
| Concat-Fusion | 2 | 22 | 0.603 | 0.301 | 1.17 | 0.5 |
| Cox+Attention | 8 | 28 | 0.612 | 0.420 | 0.75 | 0.6 |
| Cox+Concat | 2 | 22 | 0.601 | 0.300 | 1.15 | 0.5 |
| Cox+DCA | 7 | 27 | 0.611 | 0.362 | 0.74 | 0.9 |
| Cox+Gated | 6 | 26 | 0.610 | 0.306 | 0.76 | 0.7 |
| Cox+Omics | 6 | 26 | 0.642 | 0.580 | 0.46 | 0.5 |
| Cox+WSI | 2 | 22 | 0.613 | 0.463 | 0.56 | 0.4 |
| Cox+zero-residual | 99 | 100 | 0.651 | 0.668 | 0.13 | 1.6 |
| DCA | 7 | 27 | 0.630 | 0.344 | 0.73 | 0.9 |
| Gated-Fusion | 6 | 26 | 0.621 | 0.295 | 0.79 | 0.7 |
| MCAT-style | 8 | 28 | 0.659 | 0.291 | 1.98 | 7.7 |
| PORPOISE-style | 6 | 26 | 0.673 | 0.272 | 2.09 | 7.2 |
| stage1-clin | 78 | 98 | 0.642 | 0.652 | 0.32 | 2.1 |
| stage1-omics | 22 | 42 | 0.671 | 0.650 | 1.47 | 0.9 |
| stage1-wsi | 8 | 28 | 0.647 | 0.331 | 1.53 | 7.1 |

Fusion mechanism (mean over folds; dNLL_remove_* = rise in test NLL when that modality is replaced by its training mean):

| Method | alpha_wsi | alpha_omics | alpha_clin | entropy | dNLL_remove_wsi | dNLL_remove_omics | dNLL_remove_clin | norm_wsi | norm_omics | norm_clin |
|---|---|---|---|---|---|---|---|---|---|---|
| Attention-Fusion | 0.293 | 0.336 | 0.371 | 0.854 | -0.005 | 0.006 | 0.024 | 6.173 | 3.937 | 1.880 |
| Cox+Attention | 0.293 | 0.358 | 0.349 | 0.840 | 0.003 | 0.006 | 0.007 | 6.173 | 3.937 | 1.880 |
| Cox+DCA | 0.287 | 0.350 | 0.363 | 0.844 | -0.011 | 0.010 | 0.005 | 6.173 | 3.937 | 1.880 |
| DCA | 0.293 | 0.329 | 0.379 | 0.854 | -0.001 | 0.009 | 0.020 | 6.173 | 3.937 | 1.880 |