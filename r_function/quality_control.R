# quality_control.R ----
# Loads SNP presence/absence matrix and phenotypic data
# nicolas.godron@gmail.com

# INPUT: (NB: SNPs, pheno and covariates must contain the same number of individuals.)
# SNPs: Matrix of SNP absence/presence
# pheno: (1D) Matrix of phenotype
# MAF: Minor Allele Frequency filter parameter
# covariates: Covariates matrix

# OUPUT:
# Matrix of SNP absence/presence (QC-passed)
# (1D) Matrix of phenotype (QC-passed)

test_size <- function(snp_matrix, phenotypes, covariates) {
  if (! (nrow(snp_matrix) == length(phenotypes) & nrow(snp_matrix) == nrow(covariates)) ) {
    stop("Number of rows of the SNP matrix does not correspond to length of phenotypes.")
  }
  print("QC: SNPs, phenotypes and covariates matrices have corresponding sizes.")
}

remove_NA <- function(snp_matrix, phenotypes, covariates) {
  # This function:
  #     1. Gets rid of SNPs with NAs in their columns in snp_matrix
  #         and outputs names of individuals for which at least 1 SNP is missing.
  #     2. Gets rid of individuals with "NA" phenotypes (in snp_matrix and phenotypes).
  #     3. Gets rid of individuals with "NA" covariates.
  snp_matrix <- as.matrix(snp_matrix)
  phenotypes <- as.matrix(phenotypes)
  covariates <- as.matrix(covariates)
  
  n_col_NAs <- c(sum(colSums(is.na(snp_matrix) != 0)),
              sum(rowSums(is.na(phenotypes) != 0)),
              sum(rowSums(is.na(covariates) != 0)))
  
  snp_matrix <-
    snp_matrix[ , colSums(is.na(snp_matrix)) == 0]
  
  phenotypes <-
    phenotypes[complete.cases(phenotypes), ]
  
  covariates <-
    covariates[complete.cases(covariates) == 0, ]
  
  if (sum(n_col_NAs) != 0) {
    print(paste(n_col_NAs[1], "SNPs were removed."))
    print(paste(n_col_NAs[2], "individuals were removed because of missing phenotype."))
    print(paste(n_col_NAs[3], "individuals were removed because of missing covariates."))
  }
  
  n_col_not_bool <-
  
  ## To Do: QC of phenotypes here!
  # check unique(pheno) == 2 for pheno is binary

  out_NA <- list(snp_matrix, phenotypes, covariates)

  return(out_NA)
}

MAF_filter <- function(snp_matrix, minimum_allele_frequency) {
  # To Do!
}

QC <- function(SNPS, pheno, MAF, covars) {
  test_size(snp_matrix = SNPs, phenotypes = pheno, covariates = covars) # OK
  
  NA_free <- remove_NA(snp_matrix = SNPs, phenotypes = pheno) # OK for SNPs
  SNPs <- NA_free[[1]]
  pheno <- NA_free[[2]]
  
  MAF_OK <- MAF_filter(SNPs, MAF) # Filter SNPs with frequency lower than MAF
  SNPs <- MAF_OK
  QC_out <- NULL
  return(QC_out)
}

