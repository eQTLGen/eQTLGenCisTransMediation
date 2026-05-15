#!/usr/bin/env Rscript

library(data.table)
library(optparse)
library(dplyr)

# -----------------------------
# DEFINE ARGUMENTS
# -----------------------------
option_list <- list(
  make_option(
    c("--finemapped"),
    type = "character",
    help = "eQTLGen file with fine-mapped cis- and trans-eQTLs.",
    metavar = "FILE"
  )
)

parser <- OptionParser(
  option_list = option_list,
  description = "Prepare files with variant - cis-eGene - trans-eGene triplets to analyse."
)

args <- parse_args(parser)

# Functions
library(data.table)

make_batches <- function(res, target_rows_per_batch = 1000) {
  # res must contain:
  # chromosome, variant, trans_gene, cis_gene

  setorder(res, chromosome, variant)

  # count rows per variant
  vars <- res[
    , .(n_rows = .N),
    by = .(chromosome, variant)
  ]

  setorder(vars, chromosome, variant)

  # assign batches within chromosome
  vars[, batch := {
    
    b <- integer(.N)
    current_batch <- 1
    current_size  <- 0

    for (i in seq_len(.N)) {

      if (current_size > 0 &&
          current_size + n_rows[i] > target_rows_per_batch) {
        
        current_batch <- current_batch + 1
        current_size  <- 0
      }

      b[i] <- current_batch
      current_size <- current_size + n_rows[i]
    }

    b

  }, by = chromosome]

  # add batch labels back
  res <- merge(
    res,
    vars[, .(chromosome, variant, batch)],
    by = c("chromosome", "variant"),
    all.x = TRUE
  )

  return(res)
}


write_batches <- function(res, outdir = ".", prefix_cols = TRUE) {
  # writes:
  # chromosome_batch1.txt
  # chromosome_batch2.txt

  dir.create(outdir, showWarnings = FALSE)

  for (chr in unique(res$chromosome)) {
    for (b in unique(res[chromosome == chr, batch])) {

      out <- res[
        chromosome == chr & batch == b,
        .(variant, trans_gene, cis_gene)
      ]

      fname <- paste0(chr, "_batch", b, ".txt")

      fwrite(
        out,
        file.path(outdir, fname),
        sep = "\t"
      )
    }
  }
}

dt <- fread(args$finemapped)

dt[, chromosome := as.character(chromosome)]
dt[, seqid      := as.character(seqid)]

cis   <- dt[type == "cis"]
trans <- dt[type == "trans"]

# cis gene positions
cis <- unique(cis[, c(5, 7, 8, 9), with = FALSE])
trans <- unique(trans[, c(12, 5, 2, 3, 4), with = FALSE])

colnames(cis) <- c("cis_gene", "chromosome", "start", "end")
colnames(trans) <- c("variant", "trans_gene", "chromosome", "start", "end")

setkey(cis, chromosome, start, end)
setkey(trans, chromosome, start, end)

res <- foverlaps(trans, cis, type = "any", nomatch = 0L)

res <- res[, c(1, 5, 6, 2), with = FALSE]

# Remove edge cases where cis and trans are for same gene in same (large) locus
res <- res[cis_gene != trans_gene]

setorder(res, chromosome)

batches <- make_batches(res)

write_batches(batches)
