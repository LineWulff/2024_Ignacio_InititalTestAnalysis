#' R script for testing scATAC analysis w. Seurat and Signac
#' Author: Line Wulff
#' Date (created): 24-05-16
#' # Based on https://stuartlab.org/signac/articles/BMHAPBS8wk_vignette

#### ---- Initiate libraries ---- ####
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
samp <- "BM-HA107-PBS-8wk"
samp_dir <- paste(proj_data_dir, samp, sep = "")
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


BMHAPBS8wk$blacklist_fraction <- FractionCountsInRegion(
  object = BMHAPBS8wk, 
  assay = 'peaks',
  regions = blacklist_mm10
)

Annotation(BMHAPBS8wk) <- annotations
BMHAPBS8wk@meta.data$orig.ident <- samp
BMHAPBS8wk@meta.data$tissue <- unlist(str_split(BMHAPBS8wk@meta.data$orig.ident,"-"))[seq(1,nrow(BMHAPBS8wk@meta.data)*4,4)]
BMHAPBS8wk@meta.data$colonization <- unlist(str_split(BMHAPBS8wk@meta.data$orig.ident,"-"))[seq(2,nrow(BMHAPBS8wk@meta.data)*4,4)]
BMHAPBS8wk@meta.data$stimulation <- unlist(str_split(BMHAPBS8wk@meta.data$orig.ident,"-"))[seq(3,nrow(BMHAPBS8wk@meta.data)*4,4)]
BMHAPBS8wk@meta.data$timepoint <- unlist(str_split(BMHAPBS8wk@meta.data$orig.ident,"-"))[seq(4,nrow(BMHAPBS8wk@meta.data)*4,4)]


## Sample 2
samp <- "BM-PBS-PBS-8wk"
samp_dir <- paste(proj_data_dir, samp, sep = "")
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
BMPBSPBS8wk@meta.data$orig.ident <- samp
BMPBSPBS8wk@meta.data$tissue <- unlist(str_split(BMPBSPBS8wk@meta.data$orig.ident,"-"))[seq(1,nrow(BMPBSPBS8wk@meta.data)*4,4)]
BMPBSPBS8wk@meta.data$colonization <- unlist(str_split(BMPBSPBS8wk@meta.data$orig.ident,"-"))[seq(2,nrow(BMPBSPBS8wk@meta.data)*4,4)]
BMPBSPBS8wk@meta.data$stimulation <- unlist(str_split(BMPBSPBS8wk@meta.data$orig.ident,"-"))[seq(3,nrow(BMPBSPBS8wk@meta.data)*4,4)]
BMPBSPBS8wk@meta.data$timepoint <- unlist(str_split(BMPBSPBS8wk@meta.data$orig.ident,"-"))[seq(4,nrow(BMPBSPBS8wk@meta.data)*4,4)]

## after last sample remove the unnecesseary objects to save env. space
rm(chrom_assay, counts, metadata, annotations)

## Looking at the files at a glance
BMHAPBS8wk[['peaks']]
granges(BMHAPBS8wk)
head(BMHAPBS8wk@assays$peaks@data)
BMHAPBS8wk@assays$peaks@annotation
head(BMHAPBS8wk@meta.data)



#### ---- Sample QC ---- ####
## Five measures:
## Nucleosome banding pattern - calculate, saved as nucleosome_signal
## Transcription starting site - calculate, saved TSS.Enrichment
## Total number of fragments in peaks - precalc. from CR,
## Fraction of fragments in peaks - precalc. from CR, 
## Ratio reads in generic blacklist regions - precalc from CR, 

# compute nucleosome signal score per cell
# ratio of reads per cell from  mononucleosome (147-294 bp) to nucleosome free (<147 bp)
BMHAPBS8wk <- NucleosomeSignal(object = BMHAPBS8wk)

# compute TSS enrichment score per cell
BMHAPBS8wk <- TSSEnrichment(object = BMHAPBS8wk, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
BMHAPBS8wk$pct_reads_in_peaks <- BMHAPBS8wk$peak_region_fragments / BMHAPBS8wk$passed_filters * 100
BMHAPBS8wk$blacklist_ratio <- BMHAPBS8wk$blacklist_region_fragments / BMHAPBS8wk$peak_region_fragments

VlnPlot(
  object = BMHAPBS8wk,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 3
)

## Nucleosome banding patterns
# grouping cells based on their mononuc/nfr ration, here 2:1 - 2x mononucleosome bound to nfr 
BMHAPBS8wk$nucleosome_group <- ifelse(BMHAPBS8wk$nucleosome_signal > 2, 'NS > 2', 'NS < 2') 
# plotted with fragment histo 
length(BMHAPBS8wk$nucleosome_group[BMHAPBS8wk$nucleosome_group=='NS > 2'])
FragmentHistogram(BMHAPBS8wk, group.by = 'nucleosome_group', region = "chr1-1-20000000")+geom_vline(xintercept = 147)

## TSS enrivhment
BMHAPBS8wk$high.tss <- ifelse(BMHAPBS8wk$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(BMHAPBS8wk, group.by = 'high.tss') + NoLegend()
TSSPlot(BMHAPBS8wk, group.by = 'nucleosome_group') + NoLegend()


#### ---- Subset and save thresholds ---- ####
BMHAPBS8wk <- subset(
  x = BMHAPBS8wk,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)
BMHAPBS8wk