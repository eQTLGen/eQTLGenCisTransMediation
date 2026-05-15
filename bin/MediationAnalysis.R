#!/usr/bin/env Rscript

library(data.table)
library(optparse)
library(stringr)


# -----------------------------
# DEFINE ARGUMENTS
# -----------------------------
option_list <- list(
  make_option(
    c("--triplets"),
    type = "character",
    help = "File with variant - cis-eGene - trans-eGene triplets to analyse.",
    metavar = "FILE"
  ),
  make_option(
    c("--cov"),
    type = "character",
    help = "Covariates to correct models for",
    metavar = "STRING"
  ),
  make_option(
    c("--vcf_file"),
    type = "character",
    help = "Genotype vcf file to parse genotype data",
    metavar = "FILE"
  ),
  make_option(
    c("--exp_file"),
    type = "character",
    help = "Preprocessed expression matrix where gene is extracted",
    metavar = "FILE"
  ),
  make_option(
    c("--cov_file"),
    type = "character",
    help = "Covariate file where covariates are extracted",
    metavar = "FILE"
  )
)

parser <- OptionParser(
  option_list = option_list,
  description = "Run cis-trans mediation Sobel test for given trans-eQTL"
)

args <- parse_args(parser)

# Functions
# Parse genotypes

ParseGenotypes <- function(vcf, variants) {
  # variants = character vector of variant IDs

  # 1. sample IDs
  sample_ids <- fread(
    cmd = paste("bcftools query -l", vcf),
    header = FALSE
  )$V1

  # 2. build bcftools filter for multiple variants
  # example: ID="rs1" || ID="rs2" || ID="rs3"
  variant_filter <- paste0(
    'ID="', variants, '"',
    collapse = " || "
  )

  cmd <- paste0(
    "bcftools view -i '", variant_filter, "' ", vcf,
    " | bcftools query -f \"%CHROM\\t%POS\\t%ID[\\t%DS]\\n\""
  )

  dt <- fread(cmd = cmd, header = FALSE)

  # 3. assign names
  setnames(
    dt,
    c("CHROM", "POS", "ID", sample_ids)
  )

  # 4. long format
  dt <- melt(
    dt,
    id.vars = c("CHROM", "POS", "ID"),
    variable.name = "sample",
    value.name = "DS"
  )

  chr <- unique(dt$CHROM)
  pos <- unique(dt$POS)
  ID <- unique(dt$ID)

  dt <- dcast(
    dt,
    sample ~ ID,
    value.var = "DS"
  )

  return(list(chr = chr, pos = pos, ID = ID, geno = dt))

}

# ZtoP
ZtoP <- function(Z, largeZ = FALSE, log10P = TRUE) {
  if (!is.numeric(Z)) {
    message("Some of the Z-scores are not numbers! Please check why!")
    message("Converting the non-numeric vector to numeric vector.")
    Z <- as.numeric(Z)
  }

  if (largeZ == TRUE) {
    P <- log(2) + pnorm(abs(Z), lower.tail = FALSE, log.p = TRUE)

    if (largeZ == TRUE & log10P == TRUE) {
      P <- -(P * log10(exp(1)))
    }
  } else {
    P <- 2 * pnorm(abs(Z), lower.tail = FALSE)

    if (min(P) == 0) {
      P[P == 0] <- .Machine$double.xmin
      message("Some Z-score indicates very significant effect and P-value is truncated on 2.22e-308. If relevant, consider using largeZ = TRUE argument and logarithmed P-values instead.")
    }
  }

  return(P)
}

# Sobel test

# Sobel test
SobelTest <- function(trans_gene, cis_gene, variant, covariates, data){

    variant <- paste0("`", variant, "`")

    # Model for mediator
    cis_formula <- as.formula(paste(cis_gene,  "~", variant, " + ", paste(covariates, collapse = " + ")))

    # Trans-eQTL model
    trans_formula <- as.formula(paste(trans_gene, "~", variant, " + ", paste(covariates, collapse = " + ")))

    # Model for mediated trans-effect
    mediated_formula <- as.formula(paste(trans_gene, "~", variant, " + ", cis_gene, "+", paste(covariates, collapse = " + ")))

    # Fit models 
    cis_model <- lm(cis_formula, data = data) 
    trans_model <- lm(trans_formula, data = data) 
    mediated_model <- lm(mediated_formula, data = data)

    if (!str_detect(rownames(coef(summary(cis_model)))[2], "\\`")){variant <- str_remove_all(variant, "\\`")}

    cis_beta  <- coef(summary(cis_model))[variant, "Estimate"]
    cis_se <- coef(summary(cis_model))[variant, "Std. Error"]

    mediated_beta  <- coef(summary(mediated_model))[cis_gene, "Estimate"]
    mediated_se <- coef(summary(mediated_model))[cis_gene, "Std. Error"]

    # Calcuate Sobel's T-statistic
    sobel_t <- (cis_beta * mediated_beta) / sqrt(cis_beta^2 * mediated_se^2 + mediated_beta^2 * cis_se^2)
    sobel_aroian_t <- (cis_beta * mediated_beta) / sqrt(cis_beta^2 * mediated_se^2 + mediated_beta^2 * cis_se^2 + cis_se^2 * mediated_se^2)
    sobel_goodman_t <- (cis_beta * mediated_beta) / sqrt(cis_beta^2 * mediated_se^2 + mediated_beta^2 * cis_se^2 - cis_se^2 * mediated_se^2)

    # Calculate Sobel's P-values
    if(!is.na(sobel_t)){sobel_p <- ZtoP(sobel_t)} else {sobel_p <- NA}
    if(!is.na(sobel_aroian_t)){sobel_aroian_p <- ZtoP(sobel_aroian_t)} else {sobel_aroian_p <- NA}
    if(!is.na(sobel_goodman_t)){sobel_goodman_p <-ZtoP(sobel_goodman_t)} else {sobel_goodman_p <- NA}

    # Calculate the % of mediated
    trans_beta  <- coef(summary(trans_model))[variant, "Estimate"]
    trans_se  <- coef(summary(trans_model))[variant, "Std. Error"]
    trans_beta_adjusted  <- coef(summary(mediated_model))[variant, "Estimate"]
    trans_se_adjusted  <- coef(summary(mediated_model))[variant, "Std. Error"]

    prop_mediated <- (trans_beta - trans_beta_adjusted) / trans_beta

    R_cis_trans <- cor(data[[cis_gene]], data[[trans_gene]])

    results <- data.table(
        variant = variant, 
        cis_gene = cis_gene, 
        trans_gene = trans_gene,
        cor_cis_trans = R_cis_trans,
        n = nrow(data),
        cis_beta = cis_beta,
        cis_se = cis_se,
        trans_beta = trans_beta,
        trans_se = trans_se,
        trans_beta_adjusted = trans_beta_adjusted,
        trans_se_adjusted = trans_se_adjusted,
        prop_mediated = prop_mediated,
        sobel_t = sobel_t,
        sobel_p = sobel_p,
        sobel_aroian_t = sobel_aroian_t,
        sobel_aroian_p = sobel_aroian_p,
        sobel_goodman_t = sobel_goodman_t,
        sobel_goodman_p = sobel_goodman_p
        )

results$variant <- str_remove_all(results$variant, "\\`")

return(results)

}

############
# Analysis #
############

# Results template
SobelResults <- data.table(
  variant = NA, 
  cis_gene = NA, 
  trans_gene = NA, 
  cor_cis_trans = NA,
  n = NA, 
  cis_beta = NA, 
  cis_se = NA, 
  trans_beta = NA, 
  trans_se = NA, 
  trans_beta_adjusted = NA, 
  trans_se_adjusted = NA, 
  prop_mediated = NA, 
  sobel_t = NA, 
  sobel_p = NA, 
  sobel_aroian_t = NA, 
  sobel_aroian_p = NA, 
  sobel_goodman_t = NA, 
  sobel_goodman_p = NA
)[-1]

# Triplet table
triplets_to_test <- fread(args$triplets)

# Output file name
output_file_name <- paste0(args$variant, "_", args$trans_gene, "_SobelResults.txt")

# read in expression data
exp <- fread(args$exp_file)
colnames(exp)[1] <- "SampleID"
exp$SampleID <- as.character(exp$SampleID)

# Extract cis and trans genes
exp_to_analysis <- exp[, colnames(exp) %in% c("SampleID", triplets_to_test$cis_gene, triplets_to_test$trans_gene), with = FALSE]

# Remove those genes which are not present in this dataset
triplets_to_test <- triplets_to_test[cis_gene %in% colnames(exp_to_analysis)]
triplets_to_test <- triplets_to_test[trans_gene %in% colnames(exp_to_analysis)]

# read in covariate data
cov <- fread(args$cov_file)
colnames(cov)[1] <- "SampleID"
cov$SampleID <- as.character(cov$SampleID)

# split the covariates to use
covariates <- unlist(str_split(args$cov, ","))
cov <- cov[, colnames(cov) %in% c("SampleID", covariates), with = FALSE]

exp_to_analysis <- merge(exp_to_analysis, cov, by = "SampleID")

variants_to_test <- unique(triplets_to_test$variant)

# read in genotypes
geno <- ParseGenotypes(variants = variants_to_test, vcf = args$vcf_file)

geno_mat <- geno$geno
colnames(geno_mat)[1] <- "SampleID"
geno_mat$SampleID <- as.character(geno_mat$SampleID)

# Keep only those which are present in genotype data
variants_to_test <- variants_to_test[variants_to_test %in% colnames(geno_mat)]
triplets_to_test <- triplets_to_test[variant %in% variants_to_test]

nr_variants <- length(variants_to_test)
snp_iterator <- 0
gene_iterator <- seq(from = 1, to = 10000, by = 10)

input_data <- merge(geno_mat, exp_to_analysis, by = "SampleID")

if (nrow(triplets_to_test) < 1 | ncol(geno_mat) < 2 | ncol(exp_to_analysis) < length(covariates) + 2){

message("No data to analyse!")

fwrite(SobelResults, paste0(args$triplets, "_SobelResults.txt"))


} else {

for (i in 1:nrow(triplets_to_test)) {

# Sobel test
SobelResultsTemp <- SobelTest(
  trans_gene = triplets_to_test$trans_gene[i], 
  cis_gene = triplets_to_test$cis_gene[i], 
  variant = triplets_to_test$variant[i], 
  covariates = covariates, 
  data = input_data
)
SobelResults <- rbind(SobelResults, SobelResultsTemp)

if (i %in% gene_iterator){message(paste("Analysed genes:", i))}

}

fwrite(SobelResults, paste0(args$triplets, "_SobelResults.txt"))

}
