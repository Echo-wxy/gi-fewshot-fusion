# Data sources

`src/download_data.py` downloads and verifies all inputs into `$GIFS_ROOT/data`.

| Dataset | Use | Source | License |
|---|---|---|---|
| HyperKvasir, labeled images (10,662 images, 23 classes) | Study 1 | https://datasets.simula.no/hyper-kvasir/ (Borgli et al., Scientific Data 2020) | CC BY 4.0 |
| Histological images for MSI vs. MSS classification in gastrointestinal cancer, FFPE samples (TCGA tumor tiles, 224 × 224 px, 0.5 µm/px, Macenko-normalized) | Study 2 | https://doi.org/10.5281/zenodo.2530835 (Kather et al., Nature Medicine 2019) | CC BY 4.0 |
| TCGA PanCancer Atlas omics and clinical data for STAD and COADREAD | Study 2 | cBioPortal datahub, commit e1b75cd25bf7436d16719b7d76117ad8f3f66e54 | Open access (TCGA) |
| MSigDB Hallmark gene sets | Study 2 | downloaded through gseapy | CC BY 4.0 |

Pretrained encoders: DINOv2 ViT-B/14 (timm, `vit_base_patch14_dinov2.lvd142m`), UNI (MahmoodLab/UNI, gated access), Phikon (owkin/phikon), ImageNet ResNet50 (torchvision).

## Cohort flow (Study 2)

| Step | STAD | COADREAD |
|---|---|---|
| Patients in the clinical file | 440 | 594 |
| With overall survival | 435 | 590 |
| Also with tumor tiles | 280 | 358 |
| Also with expression, copy-number and mutation data (analysis cohort) | 278 | 358 |

`results/study2/cohort_STAD.csv` and `cohort_CRC.csv` list the analyzed patients (TCGA barcodes, survival, age, sex, stage, tile counts).
