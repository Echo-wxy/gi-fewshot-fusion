# Experiment 2, protocol v2 — STAD (encoder phikon, tiles <= 500)
Patients 278, deaths 108, censoring 61%; outer 5-fold x 3 seeds; inner 3-fold ensembles.

| Method | Harrell C (mean ± SD) | Uno C | IBS | Brier 12 m | Brier 24 m | Pooled C |
|---|---|---|---|---|---|---|
| Clinical-Cox | 0.646 ± 0.049 | 0.640 | 0.174 | 0.162 | 0.220 | 0.649 |
| Cox+zero-residual | 0.646 ± 0.049 | 0.640 | 0.187 | 0.178 | 0.229 | 0.649 |
| Omics-CoxNet | 0.547 ± 0.062 | 0.531 | 0.195 | 0.171 | 0.258 | 0.552 |
| Stacked-Cox | 0.645 ± 0.045 | 0.640 | 0.174 | 0.162 | 0.221 | 0.649 |
| Late-Fusion | 0.629 ± 0.035 | 0.620 | — | — | — | 0.646 |
| WSI-only | 0.538 ± 0.049 | 0.543 | 0.199 | 0.184 | 0.252 | 0.549 |
| Omics-only | 0.567 ± 0.050 | 0.557 | 0.200 | 0.183 | 0.254 | 0.572 |
| Concat-Fusion | 0.606 ± 0.043 | 0.600 | 0.194 | 0.178 | 0.241 | 0.620 |
| Attention-Fusion | 0.594 ± 0.049 | 0.587 | 0.192 | 0.180 | 0.238 | 0.613 |
| Gated-Fusion | 0.572 ± 0.055 | 0.564 | 0.196 | 0.180 | 0.244 | 0.589 |
| DCA | 0.578 ± 0.057 | 0.569 | 0.194 | 0.180 | 0.242 | 0.596 |
| PORPOISE-style | 0.580 ± 0.064 | 0.582 | 0.189 | 0.170 | 0.245 | 0.581 |
| MCAT-style | 0.541 ± 0.078 | 0.528 | 0.193 | 0.174 | 0.248 | 0.547 |
| Cox+WSI | 0.634 ± 0.038 | 0.631 | 0.187 | 0.180 | 0.226 | 0.643 |
| Cox+Omics | 0.656 ± 0.040 | 0.648 | 0.186 | 0.178 | 0.226 | 0.660 |
| Cox+Concat | 0.645 ± 0.041 | 0.638 | 0.189 | 0.180 | 0.229 | 0.648 |
| Cox+Attention | 0.651 ± 0.051 | 0.643 | 0.187 | 0.179 | 0.226 | 0.661 |
| Cox+Gated | 0.635 ± 0.062 | 0.629 | 0.189 | 0.178 | 0.231 | 0.644 |
| Cox+DCA | 0.631 ± 0.061 | 0.627 | 0.189 | 0.179 | 0.229 | 0.642 |

Primary contrast (pre-specified): Cox+DCA − Clinical-Cox, pooled out-of-fold C, patient-level bootstrap.

| Contrast | ΔC | 95% CI | Share of bootstrap ≤ 0 |
|---|---|---|---|
| Cox+DCA − Clinical-Cox | -0.007 | [-0.039, +0.027] | 0.668 |
| Cox+DCA − Cox+Attention | -0.019 | [-0.038, -0.001] | 0.978 |
| Cox+DCA − Cox+Concat | -0.006 | [-0.022, +0.009] | 0.772 |
| Cox+DCA − Cox+Gated | -0.002 | [-0.014, +0.009] | 0.635 |
| Cox+DCA − Stacked-Cox | -0.007 | [-0.037, +0.023] | 0.690 |
| Cox+DCA − Late-Fusion | -0.004 | [-0.038, +0.030] | 0.594 |
| Cox+WSI − Clinical-Cox | -0.006 | [-0.034, +0.021] | 0.663 |
| Cox+Omics − Clinical-Cox | +0.011 | [-0.022, +0.045] | 0.256 |
| DCA − Attention-Fusion | -0.017 | [-0.046, +0.011] | 0.880 |
| DCA − Concat-Fusion | -0.024 | [-0.055, +0.006] | 0.939 |
| DCA − Gated-Fusion | +0.008 | [-0.010, +0.026] | 0.196 |

Training diagnostics (median over all inner models):

| Model | Best epoch | Stop epoch | Best val NLL | Last train NLL | Grad norm | Seconds |
|---|---|---|---|---|---|---|
| Attention-Fusion | 13 | 33 | 1.123 | 0.876 | 1.03 | 0.6 |
| Concat-Fusion | 5 | 25 | 1.103 | 0.689 | 1.77 | 0.4 |
| Cox+Attention | 8 | 28 | 1.082 | 0.893 | 1.04 | 0.5 |
| Cox+Concat | 3 | 23 | 1.081 | 0.701 | 1.76 | 0.4 |
| Cox+DCA | 9 | 29 | 1.087 | 0.784 | 1.09 | 0.7 |
| Cox+Gated | 7 | 27 | 1.096 | 0.693 | 1.15 | 0.6 |
| Cox+Omics | 1 | 21 | 1.091 | 0.980 | 0.68 | 0.3 |
| Cox+WSI | 3 | 23 | 1.114 | 1.007 | 0.79 | 0.3 |
| Cox+zero-residual | 98 | 100 | 1.091 | 1.106 | 0.17 | 1.2 |
| DCA | 12 | 32 | 1.111 | 0.779 | 1.04 | 0.8 |
| Gated-Fusion | 10 | 30 | 1.122 | 0.690 | 1.13 | 0.6 |
| MCAT-style | 6 | 26 | 1.127 | 0.703 | 2.58 | 1.9 |
| PORPOISE-style | 8 | 28 | 1.135 | 0.510 | 2.83 | 1.9 |
| stage1-clin | 92 | 100 | 1.140 | 1.110 | 0.40 | 1.6 |
| stage1-omics | 10 | 30 | 1.149 | 1.021 | 2.00 | 0.5 |
| stage1-wsi | 3 | 23 | 1.176 | 0.782 | 2.90 | 1.3 |

Fusion mechanism (mean over folds; dNLL_remove_* = rise in test NLL when that modality is replaced by its training mean):

| Method | alpha_wsi | alpha_omics | alpha_clin | entropy | dNLL_remove_wsi | dNLL_remove_omics | dNLL_remove_clin | norm_wsi | norm_omics | norm_clin |
|---|---|---|---|---|---|---|---|---|---|---|
| Attention-Fusion | 0.339 | 0.314 | 0.347 | 0.922 | -0.013 | 0.011 | 0.004 | 5.262 | 4.006 | 1.707 |
| Cox+Attention | 0.339 | 0.325 | 0.336 | 0.914 | -0.013 | 0.011 | -0.008 | 5.262 | 4.006 | 1.707 |
| Cox+DCA | 0.344 | 0.317 | 0.339 | 0.911 | -0.016 | 0.004 | -0.007 | 5.262 | 4.006 | 1.707 |
| DCA | 0.339 | 0.301 | 0.360 | 0.896 | -0.018 | -0.002 | 0.006 | 5.262 | 4.006 | 1.707 |