#' R script for testing scATAC analysis w. Seurat and Signac
#' Author: Line Wulff
#' Date (created): 24-05-16
#' # Based on https://stuartlab.org/signac/articles/BMHAPBS8wk_vignette

#' #### ---- Initiate libraries ---- ####
library(ggplot2)
library(stringr)
library(ggrastr)
library(viridis)
library(scales)
library(Signac)
library(Seurat)
library(biovizBase)
library(EnsDb.Mmusculus.v79)

#### ---- variables used throughout script ---- ####
projdir <- getwd()
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
proj_data_dir <- "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/"

#### ---- Creating object from downloaded data ---- ####
## # extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
## # change to UCSC style since the data was mapped to mm10
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
# GCRm38 should be same as mm10, mm10 is the UCSC format - but same in terms of sequences, changing name, since ref gen. was mm10
genome(annotations) <- "mm10"

## Sample 1
samp_dir <- paste(proj_data_dir, "BM-HA107-PBS-8wk", sep = "")
counts <- Read10X_h5(filename = paste(samp_dir,"filtered_peak_bc_matrix.h5",sep = "/"))
motifs <- Read10X_h5(filename = paste(samp_dir,"filtered_tf_bc_matrix.h5",sep = "/"))

metadata <- read.csv(
  file = paste(samp_dir,"singlecell.csv",sep="/"),
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = paste(samp_dir,'fragments.tsv.gz',sep="/"),
  min.cells = 10,
  min.features = 200
)

BMHAPBS8wk <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  motifs = motifs,
  meta.data = metadata
)

Annotation(BMHAPBS8wk) <- annotations

## Sample 2
samp_dir <- paste(proj_data_dir, "BM-PBS-PBS-8wk", sep = "")
counts <- Read10X_h5(filename = paste(samp_dir,"filtered_peak_bc_matrix.h5",sep = "/"))
motifs <- Read10X_h5(filename = paste(samp_dir,"filtered_tf_bc_matrix.h5",sep = "/"))

metadata <- read.csv(
  file = paste(samp_dir,"singlecell.csv",sep="/"),
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = paste(samp_dir,'fragments.tsv.gz',sep="/"),
  min.cells = 10,
  min.features = 200
)

BMPBSPBS8wk <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  motifs = motifs,
  meta.data = metadata
)

Annotation(BMPBSPBS8wk) <- annotations

## after last sample remove the unnecesseary objects to save env. space
rm(chrom_assay, counts, metadata, annotations)

## Looking at the files at a glance
BMHAPBS8wk[['peaks']]
granges(BMHAPBS8wk)
head(BMHAPBS8wk@assays$peaks@data)
BMHAPBS8wk@assays$peaks@annotation


