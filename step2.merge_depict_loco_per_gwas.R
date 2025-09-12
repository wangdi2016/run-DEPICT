#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(tools) })

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: merge_depict_loco_per_gwas.R <GWAS_OUT_DIR> <OUT_PREFIX>\n",
       "Example: merge_depict_loco_per_gwas.R /.../run30gwas/gwas1 /.../run30gwas/gwas1/gwas1")
}

gwas_dir <- args[1]      # e.g. /.../run30gwas/gwas1
out_prefix <- args[2]    # e.g. /.../run30gwas/gwas1/gwas1

# helper
read_all <- function(pattern) {
  files <- list.files(gwas_dir, pattern=pattern, recursive=TRUE, full.names=TRUE)
  if (length(files)==0) return(NULL)
  rbindlist(lapply(files, fread), use.names=TRUE, fill=TRUE, idcol="SourceFile")
}

# --- Gene prioritization ---
gp <- read_all("geneprioritization\\.txt$")
if (!is.null(gp)) {
  # add CHR from path (…/chr5/…)
  gp[, CHR := sub("^.*chr([0-9XYM]+).*$", "\\1", SourceFile)]
  # normalize some common column names
  rn <- names(gp)
  ren <- function(old,new) if (old %in% rn && !(new %in% rn)) setnames(gp, old, new, skip_absent=TRUE)
  ren("GeneName","Gene"); ren("gene","Gene"); ren("Symbol","Gene")
  ren("LocusName","Locus"); ren("LocusID","Locus")
  ren("p","P"); ren("q","FDR"); ren("Fdr","FDR")

  fwrite(gp, paste0(out_prefix, ".LOCO_geneprioritization_all.tsv"), sep="\t", quote=FALSE, na="NA")

  if ("Locus" %in% names(gp)) {
    # rank within locus using FDR -> P -> Score
    if ("FDR" %in% names(gp)) gp[, rankkey := rank(FDR, ties.method="average"), by=Locus]
    else if ("P" %in% names(gp)) gp[, rankkey := rank(P, ties.method="average"), by=Locus]
    else if ("Score" %in% names(gp)) gp[, rankkey := rank(-Score, ties.method="average"), by=Locus]
    if ("rankkey" %in% names(gp)) {
      setorder(gp, Locus, rankkey)
      gp_top1 <- gp[, .SD[1], by=Locus]
      fwrite(gp_top1, paste0(out_prefix, ".LOCO_geneprioritization_top1_by_locus.tsv"),
             sep="\t", quote=FALSE, na="NA")
    }
  }
}

# --- Gene set enrichment ---
gse <- read_all("genesetenrichment\\.txt$")
if (!is.null(gse)) {
  gse[, CHR := sub("^.*chr([0-9XYM]+).*$", "\\1", SourceFile)]
  fwrite(gse, paste0(out_prefix, ".LOCO_genesetenrichment_all.tsv"), sep="\t", quote=FALSE, na="NA")
}

# --- Tissue enrichment ---
te <- read_all("tissueenrichment\\.txt$")
if (!is.null(te)) {
  te[, CHR := sub("^.*chr([0-9XYM]+).*$", "\\1", SourceFile)]
  fwrite(te, paste0(out_prefix, ".LOCO_tissueenrichment_all.tsv"), sep="\t", quote=FALSE, na="NA")
}

