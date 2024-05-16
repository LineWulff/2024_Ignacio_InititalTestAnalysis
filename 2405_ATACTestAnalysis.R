#' R script for testing scATAC analysis w. Seurat and Signac
#' Author: Line Wulff
#' Date (created): 24-05-16
#' # Based on https://stuartlab.org/signac/articles/pbmc_vignette

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


