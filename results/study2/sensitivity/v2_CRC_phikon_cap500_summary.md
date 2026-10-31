# Experiment 2, protocol v2 — CRC (encoder phikon, tiles <= 500)
Patients 358, deaths 74, censoring 79%; outer 5-fold x 3 seeds; inner 3-fold ensembles.

| Method | Harrell C (mean ± SD) | Uno C | IBS | Brier 12 m | Brier 24 m | Pooled C |
|---|---|---|---|---|---|---|
| Clinical-Cox | 0.681 ± 0.067 | 0.702 | 0.111 | 0.091 | 0.133 | 0.687 |
| Cox+zero-residual | 0.681 ± 0.067 | 0.702 | 0.113 | 0.094 | 0.132 | 0.687 |
| Omics-CoxNet | 0.496 ± 0.090 | 0.512 | 0.131 | 0.098 | 0.152 | 0.493 |
| Stacked-Cox | 0.696 ± 0.074 | 0.719 | 0.110 | 0.090 | 0.133 | 0.705 |
| Late-Fusion | 0.678 ± 0.077 | 0.699 | — | — | — | 0.689 |
| WSI-only | 0.631 ± 0.088 | 0.654 | 0.122 | 0.097 | 0.140 | 0.645 |
| Omics-only | 0.519 ± 0.087 | 0.523 | 0.123 | 0.098 | 0.142 | 0.518 |
| Concat-Fusion | 0.657 ± 0.084 | 0.672 | 0.120 | 0.096 | 0.139 | 0.675 |
| Attention-Fusion | 0.632 ± 0.097 | 0.646 | 0.121 | 0.096 | 0.140 | 0.653 |
| Gated-Fusion | 0.639 ± 0.088 | 0.660 | 0.120 | 0.096 | 0.139 | 0.663 |
| DCA | 0.630 ± 0.089 | 0.645 | 0.121 | 0.097 | 0.140 | 0.652 |
| PORPOISE-style | 0.595 ± 0.090 | 0.609 | 0.121 | 0.096 | 0.141 | 0.600 |
| MCAT-style | 0.618 ± 0.112 | 0.644 | 0.120 | 0.095 | 0.139 | 0.628 |
| Cox+WSI | 0.696 ± 0.075 | 0.722 | 0.114 | 0.093 | 0.133 | 0.710 |
| Cox+Omics | 0.643 ± 0.085 | 0.659 | 0.116 | 0.095 | 0.136 | 0.652 |
| Cox+Concat | 0.688 ± 0.074 | 0.704 | 0.116 | 0.093 | 0.135 | 0.701 |
| Cox+Attention | 0.670 ± 0.078 | 0.688 | 0.115 | 0.094 | 0.134 | 0.685 |
| Cox+Gated | 0.674 ± 0.074 | 0.700 | 0.116 | 0.094 | 0.135 | 0.689 |
| Cox+DCA | 0.674 ± 0.074 | 0.701 | 0.116 | 0.094 | 0.135 | 0.688 |

Primary contrast (pre-specified): Cox+DCA − Clinical-Cox, pooled out-of-fold C, patient-level bootstrap.

| Contrast | ΔC | 95% CI | Share of bootstrap ≤ 0 |
|---|---|---|---|
| Cox+DCA − Clinical-Cox | +0.002 | [-0.046, +0.053] | 0.469 |
| Cox+DCA − Cox+Attention | +0.004 | [-0.014, +0.022] | 0.361 |
| Cox+DCA − Cox+Concat | -0.012 | [-0.037, +0.013] | 0.834 |
| Cox+DCA − Cox+Gated | -0.001 | [-0.010, +0.010] | 0.560 |
| Cox+DCA − Stacked-Cox | -0.016 | [-0.052, +0.022] | 0.798 |
| Cox+DCA − Late-Fusion | +0.000 | [-0.038, +0.042] | 0.500 |
| Cox+WSI − Clinical-Cox | +0.024 | [-0.026, +0.081] | 0.201 |
| Cox+Omics − Clinical-Cox | -0.035 | [-0.072, +0.003] | 0.966 |
| DCA − Attention-Fusion | -0.002 | [-0.023, +0.021] | 0.573 |
| DCA − Concat-Fusion | -0.023 | [-0.053, +0.008] | 0.917 |
| DCA − Gated-Fusion | -0.011 | [-0.030, +0.009] | 0.865 |

Training diagnostics (median over all inner models):

| Model | Best epoch | Stop epoch | Best val NLL | Last train NLL | Grad norm | Seconds |
|---|---|---|---|---|---|---|
| Attention-Fusion | 10 | 30 | 0.632 | 0.420 | 0.75 | 0.6 |
| Concat-Fusion | 3 | 23 | 0.604 | 0.317 | 1.20 | 0.5 |
| Cox+Attention | 9 | 29 | 0.622 | 0.425 | 0.75 | 0.6 |
| Cox+Concat | 3 | 23 | 0.604 | 0.314 | 1.16 | 0.5 |
| Cox+DCA | 8 | 28 | 0.623 | 0.349 | 0.74 | 0.8 |
| Cox+Gated | 5 | 25 | 0.625 | 0.309 | 0.76 | 0.7 |
| Cox+Omics | 11 | 31 | 0.639 | 0.583 | 0.47 | 0.5 |
| Cox+WSI | 1 | 21 | 0.605 | 0.493 | 0.65 | 0.4 |
| Cox+zero-residual | 99 | 100 | 0.651 | 0.668 | 0.13 | 1.5 |
| DCA | 9 | 29 | 0.632 | 0.352 | 0.72 | 0.9 |
| Gated-Fusion | 8 | 28 | 0.623 | 0.303 | 0.75 | 0.7 |
| MCAT-style | 6 | 26 | 0.639 | 0.351 | 2.20 | 2.5 |
| PORPOISE-style | 7 | 27 | 0.628 | 0.295 | 2.31 | 2.3 |
| stage1-clin | 78 | 98 | 0.639 | 0.641 | 0.32 | 2.0 |
| stage1-omics | 19 | 39 | 0.671 | 0.651 | 1.48 | 0.8 |
| stage1-wsi | 8 | 28 | 0.648 | 0.379 | 2.05 | 1.9 |

Fusion mechanism (mean over folds; dNLL_remove_* = rise in test NLL when that modality is replaced by its training mean):

| Method | alpha_wsi | alpha_omics | alpha_clin | entropy | dNLL_remove_wsi | dNLL_remove_omics | dNLL_remove_clin | norm_wsi | norm_omics | norm_clin |
|---|---|---|---|---|---|---|---|---|---|---|
| Attention-Fusion | 0.353 | 0.327 | 0.320 | 0.871 | -0.001 | 0.004 | 0.025 | 7.390 | 4.003 | 1.914 |
| Cox+Attention | 0.366 | 0.342 | 0.292 | 0.859 | 0.010 | 0.001 | 0.003 | 7.390 | 4.003 | 1.914 |
| Cox+DCA | 0.364 | 0.311 | 0.325 | 0.852 | 0.008 | 0.004 | -0.002 | 7.390 | 4.003 | 1.914 |
| DCA | 0.359 | 0.287 | 0.354 | 0.857 | -0.012 | -0.010 | 0.007 | 7.390 | 4.003 | 1.914 |