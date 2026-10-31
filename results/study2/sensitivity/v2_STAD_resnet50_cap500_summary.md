# Experiment 2, protocol v2 — STAD (encoder resnet50, tiles <= 500)
Patients 278, deaths 108, censoring 61%; outer 5-fold x 3 seeds; inner 3-fold ensembles.

| Method | Harrell C (mean ± SD) | Uno C | IBS | Brier 12 m | Brier 24 m | Pooled C |
|---|---|---|---|---|---|---|
| Clinical-Cox | 0.646 ± 0.049 | 0.640 | 0.174 | 0.162 | 0.220 | 0.649 |
| Cox+zero-residual | 0.646 ± 0.049 | 0.640 | 0.187 | 0.178 | 0.229 | 0.649 |
| Omics-CoxNet | 0.547 ± 0.062 | 0.531 | 0.195 | 0.171 | 0.258 | 0.552 |
| Stacked-Cox | 0.642 ± 0.038 | 0.637 | 0.175 | 0.162 | 0.222 | 0.649 |
| Late-Fusion | 0.603 ± 0.038 | 0.592 | — | — | — | 0.616 |
| WSI-only | 0.505 ± 0.037 | 0.505 | 0.200 | 0.183 | 0.254 | 0.512 |
| Omics-only | 0.547 ± 0.063 | 0.536 | 0.201 | 0.183 | 0.258 | 0.550 |
| Concat-Fusion | 0.584 ± 0.055 | 0.581 | 0.192 | 0.177 | 0.240 | 0.594 |
| Attention-Fusion | 0.570 ± 0.045 | 0.561 | 0.193 | 0.178 | 0.241 | 0.587 |
| Gated-Fusion | 0.587 ± 0.049 | 0.579 | 0.194 | 0.179 | 0.242 | 0.610 |
| DCA | 0.580 ± 0.041 | 0.575 | 0.194 | 0.178 | 0.243 | 0.600 |
| PORPOISE-style | 0.585 ± 0.065 | 0.588 | 0.190 | 0.173 | 0.246 | 0.580 |
| MCAT-style | 0.566 ± 0.071 | 0.567 | 0.193 | 0.175 | 0.250 | 0.570 |
| Cox+WSI | 0.633 ± 0.037 | 0.626 | 0.188 | 0.180 | 0.229 | 0.640 |
| Cox+Omics | 0.647 ± 0.043 | 0.637 | 0.188 | 0.179 | 0.230 | 0.655 |
| Cox+Concat | 0.636 ± 0.045 | 0.628 | 0.191 | 0.179 | 0.233 | 0.644 |
| Cox+Attention | 0.639 ± 0.041 | 0.630 | 0.188 | 0.177 | 0.229 | 0.651 |
| Cox+Gated | 0.638 ± 0.049 | 0.630 | 0.189 | 0.178 | 0.231 | 0.652 |
| Cox+DCA | 0.636 ± 0.045 | 0.628 | 0.190 | 0.179 | 0.233 | 0.649 |

Primary contrast (pre-specified): Cox+DCA − Clinical-Cox, pooled out-of-fold C, patient-level bootstrap.

| Contrast | ΔC | 95% CI | Share of bootstrap ≤ 0 |
|---|---|---|---|
| Cox+DCA − Clinical-Cox | -0.000 | [-0.033, +0.032] | 0.503 |
| Cox+DCA − Cox+Attention | -0.002 | [-0.019, +0.015] | 0.615 |
| Cox+DCA − Cox+Concat | +0.005 | [-0.010, +0.021] | 0.265 |
| Cox+DCA − Cox+Gated | -0.003 | [-0.015, +0.009] | 0.738 |
| Cox+DCA − Stacked-Cox | -0.001 | [-0.029, +0.029] | 0.525 |
| Cox+DCA − Late-Fusion | +0.033 | [-0.004, +0.072] | 0.043 |
| Cox+WSI − Clinical-Cox | -0.008 | [-0.032, +0.016] | 0.768 |
| Cox+Omics − Clinical-Cox | +0.006 | [-0.028, +0.043] | 0.375 |
| DCA − Attention-Fusion | +0.013 | [-0.008, +0.037] | 0.122 |
| DCA − Concat-Fusion | +0.007 | [-0.018, +0.032] | 0.296 |
| DCA − Gated-Fusion | -0.010 | [-0.027, +0.006] | 0.895 |

Training diagnostics (median over all inner models):

| Model | Best epoch | Stop epoch | Best val NLL | Last train NLL | Grad norm | Seconds |
|---|---|---|---|---|---|---|
| Attention-Fusion | 13 | 33 | 1.105 | 0.895 | 1.03 | 0.6 |
| Concat-Fusion | 7 | 27 | 1.122 | 0.784 | 1.81 | 0.5 |
| Cox+Attention | 12 | 32 | 1.099 | 0.921 | 1.04 | 0.6 |
| Cox+Concat | 3 | 23 | 1.110 | 0.790 | 1.79 | 0.4 |
| Cox+DCA | 13 | 33 | 1.094 | 0.867 | 1.08 | 0.8 |
| Cox+Gated | 8 | 28 | 1.088 | 0.778 | 1.16 | 0.6 |
| Cox+Omics | 2 | 22 | 1.090 | 0.996 | 0.69 | 0.3 |
| Cox+WSI | 28 | 48 | 1.100 | 1.056 | 0.48 | 0.7 |
| Cox+zero-residual | 98 | 100 | 1.091 | 1.106 | 0.17 | 1.2 |
| DCA | 14 | 34 | 1.106 | 0.781 | 1.07 | 0.9 |
| Gated-Fusion | 12 | 32 | 1.118 | 0.710 | 1.19 | 0.7 |
| MCAT-style | 9 | 29 | 1.126 | 1.070 | 2.15 | 2.5 |
| PORPOISE-style | 9 | 29 | 1.136 | 0.639 | 2.76 | 2.3 |
| stage1-clin | 98 | 100 | 1.138 | 1.110 | 0.39 | 1.7 |
| stage1-omics | 13 | 33 | 1.154 | 1.015 | 2.00 | 0.5 |
| stage1-wsi | 19 | 39 | 1.156 | 1.092 | 1.47 | 2.7 |

Fusion mechanism (mean over folds; dNLL_remove_* = rise in test NLL when that modality is replaced by its training mean):

| Method | alpha_wsi | alpha_omics | alpha_clin | entropy | dNLL_remove_wsi | dNLL_remove_omics | dNLL_remove_clin | norm_wsi | norm_omics | norm_clin |
|---|---|---|---|---|---|---|---|---|---|---|
| Attention-Fusion | 0.335 | 0.345 | 0.320 | 0.883 | -0.015 | 0.002 | 0.007 | 2.998 | 3.945 | 1.717 |
| Cox+Attention | 0.347 | 0.356 | 0.298 | 0.866 | -0.018 | 0.003 | -0.004 | 2.998 | 3.945 | 1.717 |
| Cox+DCA | 0.336 | 0.346 | 0.318 | 0.888 | -0.035 | -0.020 | -0.026 | 2.998 | 3.945 | 1.717 |
| DCA | 0.350 | 0.322 | 0.328 | 0.895 | -0.030 | 0.020 | 0.008 | 2.998 | 3.945 | 1.717 |