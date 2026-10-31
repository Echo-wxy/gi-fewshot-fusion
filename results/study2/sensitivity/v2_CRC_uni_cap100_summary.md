# Experiment 2, protocol v2 — CRC (encoder uni, tiles <= 100)
Patients 358, deaths 74, censoring 79%; outer 5-fold x 3 seeds; inner 3-fold ensembles.

| Method | Harrell C (mean ± SD) | Uno C | IBS | Brier 12 m | Brier 24 m | Pooled C |
|---|---|---|---|---|---|---|
| Clinical-Cox | 0.681 ± 0.067 | 0.702 | 0.111 | 0.091 | 0.133 | 0.687 |
| Cox+zero-residual | 0.681 ± 0.067 | 0.702 | 0.113 | 0.094 | 0.132 | 0.687 |
| Omics-CoxNet | 0.496 ± 0.090 | 0.512 | 0.131 | 0.098 | 0.152 | 0.493 |
| Stacked-Cox | 0.695 ± 0.067 | 0.721 | 0.110 | 0.090 | 0.132 | 0.700 |
| Late-Fusion | 0.688 ± 0.100 | 0.714 | — | — | — | 0.703 |
| WSI-only | 0.653 ± 0.079 | 0.684 | 0.121 | 0.097 | 0.138 | 0.673 |
| Omics-only | 0.533 ± 0.107 | 0.535 | 0.122 | 0.098 | 0.142 | 0.535 |
| Concat-Fusion | 0.683 ± 0.093 | 0.706 | 0.118 | 0.096 | 0.136 | 0.702 |
| Attention-Fusion | 0.647 ± 0.086 | 0.664 | 0.119 | 0.097 | 0.136 | 0.675 |
| Gated-Fusion | 0.666 ± 0.100 | 0.692 | 0.119 | 0.097 | 0.137 | 0.696 |
| DCA | 0.653 ± 0.103 | 0.675 | 0.119 | 0.097 | 0.137 | 0.681 |
| PORPOISE-style | 0.610 ± 0.087 | 0.631 | 0.121 | 0.097 | 0.140 | 0.615 |
| MCAT-style | 0.636 ± 0.085 | 0.658 | 0.122 | 0.096 | 0.141 | 0.649 |
| Cox+WSI | 0.714 ± 0.086 | 0.743 | 0.113 | 0.093 | 0.130 | 0.732 |
| Cox+Omics | 0.647 ± 0.085 | 0.666 | 0.116 | 0.095 | 0.136 | 0.656 |
| Cox+Concat | 0.703 ± 0.076 | 0.729 | 0.114 | 0.094 | 0.132 | 0.712 |
| Cox+Attention | 0.685 ± 0.076 | 0.707 | 0.115 | 0.094 | 0.133 | 0.700 |
| Cox+Gated | 0.693 ± 0.083 | 0.719 | 0.115 | 0.095 | 0.132 | 0.709 |
| Cox+DCA | 0.687 ± 0.087 | 0.716 | 0.116 | 0.095 | 0.133 | 0.700 |

Primary contrast (pre-specified): Cox+DCA − Clinical-Cox, pooled out-of-fold C, patient-level bootstrap.

| Contrast | ΔC | 95% CI | Share of bootstrap ≤ 0 |
|---|---|---|---|
| Cox+DCA − Clinical-Cox | +0.013 | [-0.037, +0.065] | 0.319 |
| Cox+DCA − Cox+Attention | -0.000 | [-0.018, +0.017] | 0.507 |
| Cox+DCA − Cox+Concat | -0.012 | [-0.035, +0.011] | 0.847 |
| Cox+DCA − Cox+Gated | -0.009 | [-0.018, +0.000] | 0.970 |
| Cox+DCA − Stacked-Cox | -0.000 | [-0.038, +0.042] | 0.514 |
| Cox+DCA − Late-Fusion | -0.003 | [-0.038, +0.030] | 0.573 |
| Cox+WSI − Clinical-Cox | +0.046 | [-0.008, +0.105] | 0.048 |
| Cox+Omics − Clinical-Cox | -0.031 | [-0.069, +0.005] | 0.954 |
| DCA − Attention-Fusion | +0.006 | [-0.016, +0.028] | 0.283 |
| DCA − Concat-Fusion | -0.021 | [-0.052, +0.012] | 0.909 |
| DCA − Gated-Fusion | -0.014 | [-0.031, +0.002] | 0.966 |

Training diagnostics (median over all inner models):

| Model | Best epoch | Stop epoch | Best val NLL | Last train NLL | Grad norm | Seconds |
|---|---|---|---|---|---|---|
| Attention-Fusion | 8 | 28 | 0.623 | 0.424 | 0.76 | 0.6 |
| Concat-Fusion | 2 | 22 | 0.604 | 0.303 | 1.18 | 0.5 |
| Cox+Attention | 7 | 27 | 0.614 | 0.415 | 0.76 | 0.6 |
| Cox+Concat | 2 | 22 | 0.603 | 0.305 | 1.14 | 0.5 |
| Cox+DCA | 7 | 27 | 0.612 | 0.369 | 0.74 | 0.9 |
| Cox+Gated | 5 | 25 | 0.613 | 0.293 | 0.78 | 0.7 |
| Cox+Omics | 6 | 26 | 0.642 | 0.580 | 0.46 | 0.5 |
| Cox+WSI | 2 | 22 | 0.610 | 0.465 | 0.53 | 0.4 |
| Cox+zero-residual | 99 | 100 | 0.651 | 0.668 | 0.13 | 1.5 |
| DCA | 7 | 27 | 0.630 | 0.353 | 0.75 | 0.9 |
| Gated-Fusion | 6 | 26 | 0.624 | 0.296 | 0.78 | 0.7 |
| MCAT-style | 9 | 29 | 0.657 | 0.276 | 1.98 | 2.0 |
| PORPOISE-style | 7 | 27 | 0.681 | 0.271 | 2.13 | 1.6 |
| stage1-clin | 78 | 98 | 0.642 | 0.652 | 0.32 | 2.0 |
| stage1-omics | 22 | 42 | 0.671 | 0.650 | 1.47 | 0.8 |
| stage1-wsi | 8 | 28 | 0.645 | 0.319 | 1.53 | 1.2 |

Fusion mechanism (mean over folds; dNLL_remove_* = rise in test NLL when that modality is replaced by its training mean):

| Method | alpha_wsi | alpha_omics | alpha_clin | entropy | dNLL_remove_wsi | dNLL_remove_omics | dNLL_remove_clin | norm_wsi | norm_omics | norm_clin |
|---|---|---|---|---|---|---|---|---|---|---|
| Attention-Fusion | 0.291 | 0.334 | 0.376 | 0.845 | -0.005 | 0.006 | 0.024 | 6.075 | 3.937 | 1.880 |
| Cox+Attention | 0.283 | 0.357 | 0.360 | 0.837 | 0.005 | 0.007 | 0.007 | 6.075 | 3.937 | 1.880 |
| Cox+DCA | 0.284 | 0.349 | 0.367 | 0.839 | -0.002 | 0.008 | 0.005 | 6.075 | 3.937 | 1.880 |
| DCA | 0.289 | 0.330 | 0.381 | 0.854 | -0.003 | 0.007 | 0.019 | 6.075 | 3.937 | 1.880 |