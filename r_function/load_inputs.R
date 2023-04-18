# load_inputs.R ----
# Loads SNP presence/absence matrix and phenotypic data
# nicolas.godron@gmail.com

# INPUT:
# Directory path (char)
# SNP presence/absence tab-separated file path (char)
# Phenotype tab-separated file path (char)
# Index of phenotype column (int)

# OUPUT:
# Matrix of SNP absence/presence
# (1D) Matrix of phenotype

## Argparsing ----
### With argparse ----
# suppressPackageStartupMessages(library("argparse"))
# parser <- ArgumentParser()
# parser$add_argument("-d", "--dir", action = "store_true", default = FALSE,
#                     help = "Directory to be set as working directory in R.")
# parser$add_argument("-s", "--snp", action = "store_true", default = FALSE,
#                     help = "Path to SNP presence/absence tab-separated file (absolute, or relative to --dir).")
# parser$add_argument("-p", "--pheno", action = "store_true", default = FALSE,
#                     help = "Path to phenotype tab-separated file (absolute, or relative to --dir).")
# args <- parser$parse_args()

### Without argparse ----
arguments <- list()

# Directory in which are the SNP and phenotype files.
arguments$dir <- commandArgs(trailingOnly = TRUE)[1]

# Path to SNP presence/absence tab-separated file (absolute, or relative to --dir).
arguments$snp <- commandArgs(trailingOnly = TRUE)[2]

# Path to phenotype tab-separated file (absolute, or relative to --dir).
arguments$pheno <- commandArgs(trailingOnly = TRUE)[3]

# Index of phenotype column
arguments$pheno_index <- commandArgs(trailingOnly = TRUE)[4]

### For testing purposes ----
# arguments$dir <- "/home/nicolas/Kevin_These/Current_File/"
# arguments$snp <- "230331_matrix_snp_gene_613.csv"
# arguments$pheno <- "database_613.csv"
# arguments$pheno_index <- 3

### 
if (is.na(arguments$snp)) {
  stop("Path to SNP file was not provided (2nd trailing argument)")
}
    
if (is.na(arguments$pheno)) {
  stop("Path to phenotype file was not provided (3rd trailing argument)")
}

if (!is.na(arguments$dir)) {
  arguments$snp <- paste0(arguments$dir, arguments$snp)
  arguments$pheno <- paste0(arguments$dir, arguments$pheno)
}

## Inputs ----
# Test

loading <- function(SNP_path, pheno_path, pheno_index) {
  cat(SNP_path, "\n", pheno_path, "\n", pheno_index)
  
  SNP_matrix <- read.delim(file = SNP_path)
  SNP_matrix <- as.matrix(SNP_matrix)
  
  pheno_matrix <- read.delim(file = pheno_path)
  pheno_matrix <- pheno_matrix[pheno_index]
  pheno_matrix <- as.matrix(pheno_matrix)
  
  return(SNP_matrix, pheno_matrix)
}

matrices <- loading(arguments$snp, arguments$pheno, arguments$pheno_index)




