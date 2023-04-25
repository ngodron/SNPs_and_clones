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

test_size <- function(snp_matrix, phenotypes) { # add COVAR
  if (nrow(snp_matrix) != length(phenotypes)) {
    stop("Number of rows of the SNP matrix does not correspond to length of phenotypes.")
  }
  print("QC: SNP matrix and phenotypes have corresponding sizes.")
}

remove_NA <- function(snp_matrix, phenotypes) {
  # This function:
  #     1. Gets rid of SNPs with NAs in their columns in snp_matrix
  #         and outputs names of individuals for which at least 1 SNP is missing
  #     2. Gets rid of individuals with "NA" phenotypes (in snp_matrix and phenotypes)
  #     
  na_snp <- which(is.na(snp_matrix), arr.ind = TRUE)
  n_removed <- nrow(na_snp)
  
  ind_with_na <- names(snp_matrix)[which(is.na(snp_matrix), arr.ind = TRUE)]
  
  na_pheno <- which(is.na(phenotypes))
  
  if (n_removed) {
    print("Removing SNPs with missing values:")
    # Can stats::complete.cases() be used column-wise?
    snp_matrix <- snp_matrix[, -na_snp[,2]]
    print(paste(n_removed, "SNPs were removed"))
  }
  
  ## To Do: QC of phenotypes here!
  # check unique(pheno) == 2 for pheno is binary

  out_NA <- list(snp_matrix, phenotypes)

  return(out_NA)
}

MAF_filter <- function(snp_matrix, minimum_allele_frequency) {
  # To Do!
}

QC <- function(SNPS, pheno, MAF, covariates) {
  test_size(snp_matrix = SNPs, phenotypes = pheno) # OK
  
  NA_free <- remove_NA(snp_matrix = SNPs, phenotypes = pheno) # OK for SNPs
  SNPs <- NA_free[[1]]
  pheno <- NA_free[[2]]
  
  MAF_OK <- MAF_filter(SNPs, MAF) # Filter SNPs with frequency lower than MAF
  SNPs <- MAF_OK
  QC_out <- NULL
  return(QC_out)
}

