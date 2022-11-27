#' PAM50 Molecular Intrinsic Subtypes Based Protemoics data set
#'
#' This data set contains 190 protein with abundance information from 859 breast
#' cancer (BRCA) patients from The Cancer Genome Atlas (TCGA). Normal-like BRCA
#' are excluded due to the small sample sizes. Proteomics data are obtained and normalized from
#' The Cancer Proteome Atlas (TCPA) which applied Reverse Phase Protein
#' Arrays (RPPAs) to quantify proteins or phosphoproteins. PAM50 subtypes are treated
#' as external covariates. Notably, two subtypes, Luminal A and Luminal B BRCA, are
#' combined in one category to achieve higher power. That leads to three groups
#' with 626 BRCA patients in Luminal A and B, 75 patients in Her2-enriched and
#' 158 in Basal-like.
"pam50"

#' Proteomics data with Stemness Indices and Age in Breast Cancer
#'
#' 189 protein abundance measurement across 616 breast cancer (BRCA) patients
#' from The Cancer Genome Atlas (TCGA) are included and treated as features.
#' External covariates include two independent stemness indices (SI) which are
#' obtained based on DNA methylation (mDNAsi) and (mRNAsi), and age of patients
#' The mDNAsi and mRNAsi provide the degrees of dedifferentiation on epigenetic and
#' gene expression level respectively, which range from 0 to 1 with lower values
#' implying tendency to normal-like cells.
"si"

#' Proteomics data sets for Gynecological and Breast Cancers
#'
#' This dataset contains proteomics data for breast cancer (BRCA) and a group of
#' gynecological cancers including cervical squamous cell carcinoma
#' and endocervical adenocarcinoma (CESC), ovarian serous cystadenocarcinoma (OV)
#' and uterine corpus endometrial carcinoma (UCEC). Abundance of 189 kinds of
#' protein were measured across 1,941 patients among which 892 patients had BRCA,
#' 173 had CESC, 436 had OV and 440 had UCEC.
"pan_gynae"

#' Spatial transcriptomics data set for breast cancer
#'
#' This data set was collected from biopsy of breast cancer at
#' a thickness of 16 \eqn{\mu}m whose locations can be classified into
#' three spatial regions as tumor, intermediate, and normal based on the
#' Hematoxylin and Eosin (H&E) staining image. The size of three regions are
#' 114, 67, and 69 respectively. The data includes measurement of 100 spatially expressed
#' genes with the lowest Benjamini-Hochberg (BH) adjusted p-value by applying SPARK method
#' at 250 spot locations.
"Spatial_transcriptomics"

