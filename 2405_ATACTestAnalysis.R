#' R script for producing outputs in from Fabien's rewilded vs. SPF mice
#' Author: Line Wulff
#' Date (created): 23-08-18

#' #### ---- Initiate libraries ---- ####
library(ggplot2)
library(stringr)
library(ggrastr)
library(viridis)
library(scales)
library(Signac)

#### ---- variables used throughout script ---- ####
projdir <- getwd()
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
