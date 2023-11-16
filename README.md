## Biostatistical analysis of melanoma patients’ transcriptomic data from open database TCGA

Authors: 
* Olga Bagrova
* Darya Shavronskaya
* Dmitry Zubkov
* Evgeniy Bakin
* Maksim Kuznetsov

### Introduction

Melanoma is a type of skin cancer developing as a result of melanocytes maltransformation. Skin cutaneous melanoma (SKCM) is highly malignant and predisposed to metastasis. Median survival time of SKCM stage IV patients is only 4–9 months, and the 5-year survival rate is less than 20% [^1][^2]. There is a need in biomarkers for early diagnostics, prognosis of survival time and prediction of the response to different therapeutic strategies.

**Goal**: to analyze transcriptomic data from SKCM samples to discover tendencies in gene expression and build a model for prediction of overall survival based on the level of expression.

**Tasks**:
* study the clinical and trancriptomic data from SKCM patients that were obtained within **The Cancer Genome Atlas (TCGA)** program[^3]
* perform unsupervised hierarchical clustering, interpret the results
* perform single sample Gene Set Enrichment Analysis (ssGSEA)
* by literature review and LASSO-regularized Cox regression select transcriptomic signatures for prediction of the patients' overall survival
* develop a model of survival prediction by Cox regression

### Data and pre-processing

Clinical data for patients and samples as well as normalized RNA-seq data were obtained from **TCGA** database [^3]. Due to the substantially different nature of metastatic and primary tumor samples we used **367 metastasis samples** only. 

**Pre-processing steps**:
* log-transformation (*log2(x + 1)*)
* selection of the 1500 most variable genes by MAD (median absolute deviation)
* median-centring
* MAD-scaling 

### Results

![Grand heatmap](/Dima_Zubkov/heatmap_metastasis.png)

Unsupervised hierarchical clustering (Euclidean distance, Ward.D2 algorithm) was applied to the preprocessed data. 1500 of the most variable genes were grouped into 4 clusters labeled from **G1** to **G4**. Gene Ontology (GO) enrichment analysis showed that clusters were enriched with GO terms related to immune regulation (**G1**), response to monoamines and catecholamines (**G2**) and nervous system development (**G3**). The rest of the genes fell into cluster **G4** not associated with any specific terms. **G1** cluster could indicate the difference in a tumor immune status (so called hot, cold or altered tumors) [^4]. Interestingly, **G2** cluster composition supports a hypothesis of neuroendocrine factors’ role in melanoma pathogenesis [^5]. **G3** cluster could reflect the degree of cell dedifferentiation since malignant cells have a tendency for that and skin and neural cells have common embryonic precursors [^6][^7].

![Simple heatmap](/Dima_Zubkov/Report_DZ_metastasis_files/figure-markdown_strict/boxplot-1.png)

Samples were grouped into 4 clusters as well (**S1-S4**). Immune-related cluster **G1** had the most variable median expression between the sample clusters. **S2** had the highest **G1** expression, **S4** had the lowest one, **G1** expression in **S1** and **S3** was slightly increased and decreased respectively. Survival analysis revealed the difference in overall and progression-free survival between patients corresponding to these samples: survival time was associated positively with the level **G1** expression (*HR = 0.64, 95% CI: 0.52—0.78* for overall survival).

![](/Dima_Zubkov/Report_DZ_metastasis_files/figure-markdown_strict/survival-1.png) ![](/Dima_Zubkov/Report_DZ_metastasis_files/figure-markdown_strict/survival-2.png)

However, it would be impractical to conduct the whole RNA-seq analysis for a patient-derived sample to predict their expected survival time. The preferable way would be to have a limited set of genes for the analysis. Therefore, we tried to find such genes.

We have reviewed the previous works on prognostic biomarkers for SKCM. 10 gene signatures were chosen according to two papers [^8][^9]: *SUCO*, *BTN3A1*, *LINC02908 (C9orf139)*, *MIR3667HG (C22orf34)*, *CCL4*, *CXCL10*, *CCL5*, *GZMB*, *C1QA*, *C1QB*. Next, we built our own LASSO-regularized Cox regression model on 1500 of the most variably expressed genes  to select the variables with the best predictive ability. Different preprocessing methods resulted in the several gene sets (complete analysis see [here](/Olga_Bagrova/3_GLMNET/Glmnet.md)).

![Upset plot](/Olga_Bagrova/3_GLMNET/Glmnet_files/figure-html/unnamed-chunk-86-1.png)

Each of the gene signatures found in the literature or detected by the modeling, were used for survival analysis.
All of the mRNAs/lncRNAs taken from the reviews were confirmed to have a strong prognostic potential. 

![HR literature](/Dasha_Shavronskaya/Shavronskaya-Cox-regression_files/figure-markdown_github/unnamed-chunk-12-1.png)

Of the model-derived genes, *KLRD1*, *GBP4*, *CCL8*, *GBP1P1*, *HSPA7* expression level positively correlated with survivability, while for *OCA2* and *WNK2* correlation was negative. 

![HR LASSO](/Dasha_Shavronskaya/Shavronskaya-Cox-regression_files/figure-markdown_github/unnamed-chunk-13-1.png)

*KLRD1*, *GBP4*, *CCL8*, *GBP1P1*, *OCA2* were described before as such prognostic biomarkers for SKCM [^10][^11][^12][^13]. *WNK2* showed the same tendencies in case of hepatocellular carcinoma [^14], while *HSPA7* had the opposite trend in glioblastoma [^15].

<img src="/Dasha_Shavronskaya/Shavronskaya-Cox-regression_files/figure-markdown_github/unnamed-chunk-24-1.png" width="425"/> <img src="/Dasha_Shavronskaya/Shavronskaya-Cox-regression_files/figure-markdown_github/unnamed-chunk-20-1.png" width="425"/>
<img src="/Dasha_Shavronskaya/Shavronskaya-Cox-regression_files/figure-markdown_github/unnamed-chunk-19-1.png" width="425"/> <img src="/Dasha_Shavronskaya/Shavronskaya-Cox-regression_files/figure-markdown_github/unnamed-chunk-23-1.png" width="425"/> 
<img src="/Dasha_Shavronskaya/Shavronskaya-Cox-regression_files/figure-markdown_github/unnamed-chunk-18-1.png" width="425"/> <img src="/Dasha_Shavronskaya/Shavronskaya-Cox-regression_files/figure-markdown_github/unnamed-chunk-22-1.png" width="425"/>
<img src="/Dasha_Shavronskaya/Shavronskaya-Cox-regression_files/figure-markdown_github/unnamed-chunk-21-1.png" width="425"/>

### Conclusion

* In transcriptomic data there were detected gene clusters related to immune regulation, response to monoamines and catecholamines and nervous system development
* Sample clusterisation by transcriptomics data resulted in sufficient difference in overall survival between sample clusters 
* As a result of LASSO-regularized Cox regression 7 gene signatures were detected: *KLRD1*, *GBP4*, *CCL8*, *GBP1P1*, *HSPA7*,  *OCA2* and *WNK2*
* A new prognostic model for SKCM patients' overall survival was developed based on both clinical and gene expression data
* HR for each of the gene signatures was estimated
* Further plans include the presentation of the prognostic model as a Shiny app

### Literature

[^1]:	Aubuchon, M. M. F., Bolt, L. J. J., Janssen-Heijnen, M. L. G., Verleisdonk-Bolhaar, S. T. H. P., Van Marion, A., & Van Berlo, C. L. H. (2017). Epidemiology, management and survival outcomes of primary cutaneous melanoma: A ten-year overview. Acta Chirurgica Belgica, 117(1), 29-35. doi: 10.1080/00015458.2016.1242214
[^2]:	Leonardi, G. C., Falzone, L., Salemi, R., Zanghì, A., Spandidos, D. A., McCubrey, J. A., . . . Libra, M. (2018). Cutaneous melanoma: From pathogenesis to therapy (Review) International Journal of Oncology (Vol. 52, pp. 1071-1080): Spandidos Publications.
[^3]:	http://www.cbioportal.org/study/summary?id=skcm_tcga_pan_can_atlas_2018
[^4]:	Galon, J., & Bruni, D. (2019). Approaches to treat immune hot, altered and cold tumours with combination immunotherapies Nature Reviews Drug Discovery (Vol. 18, pp. 197-218): Nature Publishing Group.
[^5]:	Scheau, C., Draghici, C., Ilie, M. A., Lupu, M., Solomon, I., Tampa, M., . . . Caruntu, C. (2021). Neuroendocrine Factors in Melanoma Pathogenesis. Cancers, 13, 2277-2277. doi: 10.3390/cancers
[^6]:	Belote, R. L., Le, D., Maynard, A., Lang, U. E., Sinclair, A., Lohman, B. K., . . . Judson-Torres, R. L. (2021). Human melanocyte development and melanoma dedifferentiation at single-cell resolution. Nature Cell Biology, 23(9), 1035-1047. doi: 10.1038/s41556-021-00740-8
[^7]:	Vandamme, N., & Berx, G. (2019). From neural crest cells to melanocytes: cellular plasticity during development and beyond Cellular and Molecular Life Sciences (Vol. 76, pp. 1919-1934): Birkhauser Verlag AG.
[^8]:	Tang, Y., Feng, H., Zhang, L., Qu, C., Li, J., Deng, X., . . . Peng, X. (2022). A novel prognostic model for cutaneous melanoma based on an immune-related gene signature and clinical variables. Scientific Reports, 12(1). doi: 10.1038/s41598-022-23475-4
[^9]:	Liang, Z., Pan, L., Shi, J., & Zhang, L. (2022). C1QA, C1QB, and GZMB are novel prognostic biomarkers of skin cutaneous melanoma relating tumor microenvironment. Scientific Reports, 12(1). doi: 10.1038/s41598-022-24353-9
[^10]:	Cursons, J., Souza-Fonseca-Guimaraes, F., Foroutan, M., Anderson, A., Hollande, F., Hediyeh-Zadeh, S., . . . Davis, M. J. (2019). A gene signature predicting natural killer cell infiltration and improved survival in melanoma patients. Cancer Immunology Research, 7(7), 1162-1174. doi: 10.1158/2326-6066.CIR-18-0500
[^11]:	Wang, Q., Wang, X., Liang, Q., Wang, S., Xiwen, L., Pan, F., . . . Li, D. (2018). Distinct prognostic value of mRNA expression of guanylate-binding protein genes in skin cutaneous melanoma. Oncology Letters, 15(5), 7914-7922. doi: 10.3892/ol.2018.8306
[^12]:	Azevedo, H., Pessoa, G. C., De Luna Vitorino, F. N., Nsengimana, J., Newton-Bishop, J., Reis, E. M., . . . Jasiulionis, M. G. (2020). Gene co-expression and histone modification signatures are associated with melanoma progression, epithelial-to-mesenchymal transition, and metastasis. Clinical Epigenetics, 12(1). doi: 10.1186/s13148-020-00910-9
[^13]:	Yang, P., Chen, W., Xu, H., Yang, J., Jiang, J., Jiang, Y., & Xu, G. (2021). Correlation of CCL8 expression with immune cell infiltration of skin cutaneous melanoma: potential as a prognostic indicator and therapeutic pathway. Cancer Cell International, 21(1). doi: 10.1186/s12935-021-02350-8
[^14]:	Ho, Y. J., Chang, J., Yeh, K. T., Gong, Z., Lin, Y. M., & Lu, J. W. (2020). Prognostic and clinical implications of WNK lysine deficient protein kinase 1 expression in patients with hepatocellular carcinoma. In Vivo, 34(5), 2631-2640. doi: 10.21873/invivo.12081
[^15]:	Zhao, R., Li, B., Zhang, S., He, Z., Pan, Z., Guo, Q., . . . Li, G. (2021). The N6-Methyladenosine-Modified Pseudogene HSPA7 Correlates With the Tumor Microenvironment and Predicts the Response to Immune Checkpoint Therapy in Glioblastoma. Frontiers in Immunology, 12. doi: 10.3389/fimmu.2021.653711
