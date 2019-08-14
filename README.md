<center> <h2> iPath </h2> </center>

-------------------
Comprehensive pan-cancer analysis of the transcription of pre-defined gene sets at the individual level reveals novel biomarkers that links to clinical prognosis. **iPath** is an R package designated for the implementation of iPath algorithm, which is used to selected significiant pathways. These significant pathways demonstrate survival difference in survival for *TCGA* data (14 cancer types `Fig.a`). A schemetic overview of the glorithm is shown in (`Fig.b`). All the copyrights are explaned by Kenong Su <kenong.su@emory.edu> and Dr. Zhaohui (Steve) Qin <zhaohui.qin@emory.edu>.

![workflow](/assets/Fig.png)

### 1. Software Installation
```
library(devtools)
install_github("suke18/iPath")
R CMD INSTALL iPath_0.1.0.tar.gz # Alternatively, use this command line in the terminal.
```

### 2. Code Snippets
**(1). load expression data, gene set database (GSDB), and clinical data**
```r
library(iPath)
data("Y") # corresponding to BRCA cancer type
data("GSDB_example")
data("BRCA_cli")
```
**(2). iES score calculation per pathway per patient**
```r
iES_mat = iES_cal(Y, GSDB = GSDB_example)
```
**(3). iPath survival investigation**
```r
iPath_rslt = iES_surv(GSDB = GSDB_example, iES_mat = iES_mat,cli = BRCA_cli, qval=F)
```
**(4). Demonstration Figures**
```r
gs_str = "FARMER_BREAST_CANCER_APOCRINE_VS_LUMINAL"
water_fall(iES_mat = iES_mat, gs_str = gs_str) # waterfall plot for tumor and normal patients iES scores
density_fall(iES_mat = iES_mat, gs_str = gs_str) # densityfall plot for tumor and normal patients iES scores
iES_surv_one(GSDB_example, iES_mat, BRCA_cli, gs_str = gs_str) # one case of the survival analysis
```
