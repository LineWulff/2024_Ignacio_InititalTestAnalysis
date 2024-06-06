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

#### ---- variables used throughout script ---- ####
projdir <- getwd()
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
proj_data_dir <- "/Volumes/Promise RAID/Line/projects/24_TI_IgnacioWulff/samples/BM-HA107-PBS-8wk"

#### ---- Creating object from downloaded data ---- ####
counts <- Read10X_h5(filename = paste(proj_data_dir,"filtered_peak_bc_matrix.h5",sep = "/"))

metadata <- read.csv(
  file = paste(proj_data_dir,"singlecell.csv",sep="/"),
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = paste(proj_data_dir,'fragments.tsv.gz',sep="/"),
  min.cells = 10,
  min.features = 200
)

BMHAPBS8wk <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
