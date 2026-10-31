# Experiment 2, protocol v2 — CRC (encoder resnet50, tiles <= 500)
Patients 358, deaths 74, censoring 79%; outer 5-fold x 3 seeds; inner 3-fold ensembles.

| Method | Harrell C (mean ± SD) | Uno C | IBS | Brier 12 m | Brier 24 m | Pooled C |
|---|---|---|---|---|---|---|
| Clinical-Cox | 0.681 ± 0.067 | 0.702 | 0.111 | 0.091 | 0.133 | 0.687 |
| Cox+zero-residual | 0.681 ± 0.067 | 0.702 | 0.113 | 0.094 | 0.132 | 0.687 |
| Omics-CoxNet | 0.496 ± 0.090 | 0.512 | 0.131 | 0.098 | 0.152 | 0.493 |
| Stacked-Cox | 0.675 ± 0.077 | 0.695 | 0.111 | 0.091 | 0.134 | 0.682 |
| Late-Fusion | 0.634 ± 0.080 | 0.663 | — | — | — | 0.651 |
| WSI-only | 0.555 ± 0.093 | 0.585 | 0.122 | 0.098 | 0.141 | 0.570 |
| Omics-only | 0.533 ± 0.107 | 0.535 | 0.122 | 0.098 | 0.142 | 0.535 |
| Concat-Fusion | 0.645 ± 0.085 | 0.669 | 0.119 | 0.096 | 0.138 | 0.662 |
| Attention-Fusion | 0.644 ± 0.095 | 0.669 | 0.117 | 0.095 | 0.135 | 0.660 |
| Gated-Fusion | 0.640 ± 0.095 | 0.657 | 0.118 | 0.096 | 0.137 | 0.657 |
| DCA | 0.637 ± 0.105 | 0.650 | 0.118 | 0.095 | 0.136 | 0.651 |
| PORPOISE-style | 0.531 ± 0.073 | 0.539 | 0.123 | 0.098 | 0.144 | 0.527 |
| MCAT-style | 0.582 ± 0.066 | 0.606 | 0.124 | 0.097 | 0.144 | 0.595 |
| Cox+WSI | 0.695 ± 0.082 | 0.719 | 0.114 | 0.095 | 0.133 | 0.709 |
| Cox+Omics | 0.647 ± 0.085 | 0.666 | 0.116 | 0.095 | 0.136 | 0.656 |
| Cox+Concat | 0.670 ± 0.076 | 0.690 | 0.116 | 0.094 | 0.136 | 0.682 |
| Cox+Attention | 0.678 ± 0.088 | 0.701 | 0.113 | 0.093 | 0.132 | 0.693 |
| Cox+Gated | 0.672 ± 0.089 | 0.693 | 0.116 | 0.095 | 0.135 | 0.686 |
| Cox+DCA | 0.674 ± 0.092 | 0.694 | 0.115 | 0.094 | 0.134 | 0.690 |

Primary contrast (pre-specified): Cox+DCA − Clinical-Cox, pooled out-of-fold C, patient-level bootstrap.

| Contrast | ΔC | 95% CI | Share of bootstrap ≤ 0 |
|---|---|---|---|
| Cox+DCA − Clinical-Cox | +0.003 | [-0.035, +0.046] | 0.467 |
| Cox+DCA − Cox+Attention | -0.004 | [-0.018, +0.011] | 0.696 |
| Cox+DCA − Cox+Concat | +0.007 | [-0.016, +0.033] | 0.282 |
| Cox+DCA − Cox+Gated | +0.004 | [-0.006, +0.013] | 0.222 |
| Cox+DCA − Stacked-Cox | +0.007 | [-0.028, +0.047] | 0.362 |
| Cox+DCA − Late-Fusion | +0.038 | [-0.001, +0.081] | 0.029 |
| Cox+WSI − Clinical-Cox | +0.023 | [-0.016, +0.070] | 0.140 |
| Cox+Omics − Clinical-Cox | -0.031 | [-0.069, +0.005] | 0.954 |
| DCA − Attention-Fusion | -0.010 | [-0.031, +0.009] | 0.842 |
| DCA − Concat-Fusion | -0.012 | [-0.047, +0.021] | 0.757 |
| DCA − Gated-Fusion | -0.006 | [-0.023, +0.009] | 0.773 |

Training diagnostics (median over all inner models):

| Model | Best epoch | Stop epoch | Best val NLL | Last train NLL | Grad norm | Seconds |
|---|---|---|---|---|---|---|
| Attention-Fusion | 11 | 31 | 0.615 | 0.524 | 0.75 | 0.7 |
| Concat-Fusion | 5 | 25 | 0.612 | 0.425 | 1.23 | 0.5 |
| Cox+Attention | 10 | 30 | 0.611 | 0.511 | 0.75 | 0.6 |
| Cox+Concat | 3 | 23 | 0.617 | 0.438 | 1.19 | 0.5 |
| Cox+DCA | 9 | 29 | 0.615 | 0.474 | 0.76 | 0.9 |
| Cox+Gated | 8 | 28 | 0.617 | 0.404 | 0.79 | 0.8 |
| Cox+Omics | 6 | 26 | 0.642 | 0.580 | 0.46 | 0.5 |
| Cox+WSI | 15 | 35 | 0.622 | 0.614 | 0.37 | 0.6 |
| Cox+zero-residual | 99 | 100 | 0.651 | 0.668 | 0.13 | 1.5 |
| DCA | 11 | 31 | 0.624 | 0.441 | 0.78 | 0.9 |
| Gated-Fusion | 9 | 29 | 0.630 | 0.393 | 0.83 | 0.8 |
| MCAT-style | 4 | 24 | 0.635 | 0.689 | 2.07 | 2.7 |
| PORPOISE-style | 7 | 27 | 0.643 | 0.343 | 2.56 | 2.9 |
| stage1-clin | 78 | 98 | 0.642 | 0.652 | 0.32 | 2.0 |
| stage1-omics | 22 | 42 | 0.671 | 0.650 | 1.47 | 0.8 |
| stage1-wsi | 36 | 56 | 0.654 | 0.656 | 1.02 | 4.6 |

Fusion mechanism (mean over folds; dNLL_remove_* = rise in test NLL when that modality is replaced by its training mean):

| Method | alpha_wsi | alpha_omics | alpha_clin | entropy | dNLL_remove_wsi | dNLL_remove_omics | dNLL_remove_clin | norm_wsi | norm_omics | norm_clin |
|---|---|---|---|---|---|---|---|---|---|---|
| Attention-Fusion | 0.289 | 0.337 | 0.373 | 0.832 | -0.002 | 0.002 | 0.033 | 4.244 | 3.937 | 1.880 |
| Cox+Attention | 0.295 | 0.372 | 0.333 | 0.845 | -0.000 | 0.004 | 0.006 | 4.244 | 3.937 | 1.880 |
| Cox+DCA | 0.285 | 0.380 | 0.335 | 0.840 | -0.001 | 0.002 | -0.002 | 4.244 | 3.937 | 1.880 |
| DCA | 0.300 | 0.343 | 0.357 | 0.835 | -0.005 | 0.001 | 0.026 | 4.244 | 3.937 | 1.880 |