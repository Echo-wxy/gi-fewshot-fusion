# Experiment 2, protocol v2 — CRC (encoder uni, tiles <= 500)
Patients 358, deaths 74, censoring 79%; outer 5-fold x 3 seeds; inner 3-fold ensembles.

| Method | Harrell C (mean ± SD) | Uno C | IBS | Brier 12 m | Brier 24 m | Pooled C |
|---|---|---|---|---|---|---|
| Clinical-Cox | 0.681 ± 0.067 | 0.702 | 0.111 | 0.091 | 0.133 | 0.687 |
| Cox+zero-residual | 0.681 ± 0.067 | 0.702 | 0.113 | 0.094 | 0.132 | 0.687 |
| Omics-CoxNet | 0.496 ± 0.090 | 0.512 | 0.131 | 0.098 | 0.152 | 0.493 |
| Stacked-Cox | 0.697 ± 0.066 | 0.721 | 0.110 | 0.090 | 0.132 | 0.702 |
| Late-Fusion | 0.689 ± 0.101 | 0.714 | — | — | — | 0.702 |
| WSI-only | 0.653 ± 0.081 | 0.682 | 0.120 | 0.097 | 0.137 | 0.673 |
| Omics-only | 0.533 ± 0.107 | 0.535 | 0.122 | 0.098 | 0.142 | 0.535 |
| Concat-Fusion | 0.680 ± 0.094 | 0.702 | 0.118 | 0.096 | 0.136 | 0.700 |
| Attention-Fusion | 0.647 ± 0.087 | 0.662 | 0.119 | 0.096 | 0.135 | 0.673 |
| Gated-Fusion | 0.671 ± 0.106 | 0.694 | 0.119 | 0.097 | 0.136 | 0.704 |
| DCA | 0.657 ± 0.104 | 0.675 | 0.119 | 0.097 | 0.136 | 0.685 |
| PORPOISE-style | 0.612 ± 0.077 | 0.630 | 0.121 | 0.096 | 0.140 | 0.621 |
| MCAT-style | 0.628 ± 0.090 | 0.656 | 0.122 | 0.096 | 0.140 | 0.642 |
| Cox+WSI | 0.716 ± 0.086 | 0.744 | 0.112 | 0.093 | 0.129 | 0.733 |
| Cox+Omics | 0.647 ± 0.085 | 0.666 | 0.116 | 0.095 | 0.136 | 0.656 |
| Cox+Concat | 0.701 ± 0.079 | 0.726 | 0.114 | 0.094 | 0.133 | 0.713 |
| Cox+Attention | 0.688 ± 0.077 | 0.711 | 0.114 | 0.094 | 0.132 | 0.705 |
| Cox+Gated | 0.699 ± 0.082 | 0.723 | 0.115 | 0.095 | 0.132 | 0.713 |
| Cox+DCA | 0.690 ± 0.084 | 0.717 | 0.115 | 0.095 | 0.133 | 0.704 |

Primary contrast (pre-specified): Cox+DCA − Clinical-Cox, pooled out-of-fold C, patient-level bootstrap.

| Contrast | ΔC | 95% CI | Share of bootstrap ≤ 0 |
|---|---|---|---|
| Cox+DCA − Clinical-Cox | +0.018 | [-0.032, +0.072] | 0.258 |
| Cox+DCA − Cox+Attention | -0.001 | [-0.017, +0.016] | 0.529 |
| Cox+DCA − Cox+Concat | -0.008 | [-0.031, +0.017] | 0.749 |
| Cox+DCA − Cox+Gated | -0.009 | [-0.018, -0.000] | 0.982 |
| Cox+DCA − Stacked-Cox | +0.003 | [-0.036, +0.045] | 0.466 |
| Cox+DCA − Late-Fusion | +0.002 | [-0.033, +0.037] | 0.460 |
| Cox+WSI − Clinical-Cox | +0.047 | [-0.007, +0.107] | 0.047 |
| Cox+Omics − Clinical-Cox | -0.031 | [-0.069, +0.005] | 0.954 |
| DCA − Attention-Fusion | +0.012 | [-0.011, +0.035] | 0.135 |
| DCA − Concat-Fusion | -0.014 | [-0.046, +0.017] | 0.825 |
| DCA − Gated-Fusion | -0.019 | [-0.036, -0.003] | 0.989 |

Training diagnostics (median over all inner models):

| Model | Best epoch | Stop epoch | Best val NLL | Last train NLL | Grad norm | Seconds |
|---|---|---|---|---|---|---|
| Attention-Fusion | 8 | 28 | 0.623 | 0.418 | 0.75 | 0.6 |
| Concat-Fusion | 2 | 22 | 0.602 | 0.307 | 1.17 | 0.5 |
| Cox+Attention | 7 | 27 | 0.609 | 0.419 | 0.75 | 0.6 |
| Cox+Concat | 2 | 22 | 0.600 | 0.306 | 1.13 | 0.4 |
| Cox+DCA | 7 | 27 | 0.615 | 0.372 | 0.74 | 0.8 |
| Cox+Gated | 6 | 26 | 0.608 | 0.307 | 0.77 | 0.7 |
| Cox+Omics | 6 | 26 | 0.642 | 0.580 | 0.46 | 0.4 |
| Cox+WSI | 2 | 22 | 0.611 | 0.463 | 0.56 | 0.4 |
| Cox+zero-residual | 99 | 100 | 0.651 | 0.668 | 0.13 | 1.5 |
| DCA | 7 | 27 | 0.635 | 0.357 | 0.74 | 0.8 |
| Gated-Fusion | 6 | 26 | 0.623 | 0.295 | 0.78 | 0.7 |
| MCAT-style | 8 | 28 | 0.655 | 0.284 | 1.98 | 3.0 |
| PORPOISE-style | 7 | 27 | 0.683 | 0.257 | 2.13 | 2.7 |
| stage1-clin | 78 | 98 | 0.642 | 0.652 | 0.32 | 1.9 |
| stage1-omics | 22 | 42 | 0.671 | 0.650 | 1.47 | 0.8 |
| stage1-wsi | 8 | 28 | 0.645 | 0.325 | 1.52 | 2.3 |

Fusion mechanism (mean over folds; dNLL_remove_* = rise in test NLL when that modality is replaced by its training mean):

| Method | alpha_wsi | alpha_omics | alpha_clin | entropy | dNLL_remove_wsi | dNLL_remove_omics | dNLL_remove_clin | norm_wsi | norm_omics | norm_clin |
|---|---|---|---|---|---|---|---|---|---|---|
| Attention-Fusion | 0.294 | 0.336 | 0.370 | 0.854 | -0.010 | 0.007 | 0.024 | 6.208 | 3.937 | 1.880 |
| Cox+Attention | 0.293 | 0.356 | 0.351 | 0.837 | -0.001 | 0.007 | 0.009 | 6.208 | 3.937 | 1.880 |
| Cox+DCA | 0.287 | 0.347 | 0.366 | 0.841 | -0.008 | 0.010 | 0.005 | 6.208 | 3.937 | 1.880 |
| DCA | 0.290 | 0.328 | 0.382 | 0.854 | -0.010 | 0.030 | 0.018 | 6.208 | 3.937 | 1.880 |