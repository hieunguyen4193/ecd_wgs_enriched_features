gc()
rm(list = ls())

new.pkgs <- c("tidyverse")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }
}

library(dplyr)
library(tidyverse)
library(ggplot2)
library(vroom)
library(limma)
library(GenomicRanges)
library(comprehenr)

options(scipen = 999)
path.to.save.output <- "/home/hieunguyen/src/ecd_wgs_enriched_features/enrichFeature_pipeline/bed_files/bin600/"
path.to.storage <- "/media/HNSD01/storage/resources"

# read in nucleosome map from healthy samples
nucleosome.map <- read.csv(file.path(path.to.storage, "rpr_map_EXP0779.sorted.bed"), header = FALSE, sep = "\t")
colnames(nucleosome.map) <- c("chrom", "start", "end", "region", "V5", "strand", "mid.point.start", "mid.point.end")

nucleosome.map <- nucleosome.map %>% 
  rowwise() %>%
  mutate(new.start = mid.point.start - 300) %>%
  mutate(new.end = mid.point.end + 300) %>%
  mutate(new.len = new.end - new.start)

nucleosome.map$new.start <- as.character(nucleosome.map$new.start)
nucleosome.map$new.end <- as.character(nucleosome.map$new.end)

write.table(nucleosome.map[, c("chrom", "new.start", "new.end", "new.len")], 
            file.path(path.to.save.output, "nucleosome_map_bin600.bed"),
            row.names = FALSE, 
            col.names = FALSE, 
            sep = "\t", 
            quote = FALSE)

##### EOF