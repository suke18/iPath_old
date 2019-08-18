<center> <h3> iPath </h3> </center>

-------------------
Comprehensive pan-cancer analysis of the transcription of pre-defined gene sets at the individual level reveals novel biomarkers that links to clinical prognosis. **iPath** is an R package designated for the implementation of iPath algorithm, which is used to selected significiant pathways. These significant pathways demonstrate survival difference in survival for *TCGA* data (14 cancer types `Fig.a`). A schemetic overview of the glorithm is shown in (`Fig.b`). All the copyrights are explained by Kenong Su <kenong.su@emory.edu> and Dr. Zhaohui (Steve) Qin <zhaohui.qin@emory.edu>.

![workflow](/assets/Fig.png)

### 1. Software Installation
* Version 0.1.0 released
    + dependent bioconductor packages: Biobase (2.42.0), qvalue(2.14.1)
    + It can work on Windows, Mac and Linux platforms

```
library(devtools)
install_github("suke18/iPath")
R CMD INSTALL iPath_0.1.0.tar.gz # Alternatively, use this command line in the terminal.
```

### 2. Code Snippets
**(1). load expression data, gene set database (GSDB), and clinical data**
```r
library(iPath)
data("BRCA_exprs") # corresponding to BRCA cancer type
data("GSDB_example")
data("BRCA_cli")
```
**(2). iES score calculation per pathway per patient**
```r
iES_mat = iES_cal(BRCA_exprs, GSDB = GSDB_example)
```
**(3). iPath survival investigation**
```r
iPath_rslt = iES_surv(GSDB = GSDB_example, iES_mat = iES_mat,cli = BRCA_cli, qval=F)
```
**(4). Demonstration Figures**
```r
gs_str = "FARMER_BREAST_CANCER_APOCRINE_VS_LUMINAL"
water_fall(iES_mat = iES_mat, gs_str = gs_str)
density_fall(iES_mat = iES_mat, gs_str = gs_str)
iES_surv_one(GSDB_example, iES_mat, BRCA_cli, gs_str = gs_str)
```

### 2. Plots Illustration
**(1). waterfall plot** <br/>
For a specific pathway, iPath draws the waterfall plot which contains the iES scores for tumor and normal samples. The tSNE plot of iES score for the *C2* GSDB from *MSigDB* across these 14 cancer types is illustrated in (`Fig.c`).

![waterfall](/assets/Waterfall.jpeg)

**(2). densityfall plot**<br/>
For a specific pathway, iPath draws the density plots of the iES scores for tumor and normal samples. Because of the heterogeneity also wildly exists in normal samples from TCGA, iPath considers a mixture model fitting into normal samples.

![densityfall](/assets/densityfall.jpeg)

**(3). survival plot**<br/>
After classifying tumor samples into two groups: normal-like and perturbed, iPath performs the survival anlaysis based on the two groups.

![survivalone](/assets/Survival.jpeg)
