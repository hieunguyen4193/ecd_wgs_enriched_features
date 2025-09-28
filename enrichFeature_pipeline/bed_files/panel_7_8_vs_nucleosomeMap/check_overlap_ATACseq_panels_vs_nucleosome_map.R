gc()
rm(list = ls())

library(dplyr)
library(tidyverse)
library(ggplot2)
library(vroom)
library(limma)
library(GenomicRanges)
library(comprehenr)

inputdir <- "/media/hieunguyen/HNSD01/storage/TCGA_ATAC-seq"
outdir <- "/media/hieunguyen/HNSD01/outdir"

path.to.main.src <- "/media/hieunguyen/HNSD01/src/ATAC-seq-TCGA-data-analysis"
dir.create(file.path(path.to.main.src, "TCGA_atac_seq_beds"), showWarnings = FALSE, recursive = TRUE)  

PROJECT <- "ATAC-seq-TCGA-data-analysis"

path.to.storage <- "/media/hieunguyen/HNSD01/storage/TCGA_ATAC-seq"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.main.src <- "/media/hieunguyen/HNSD01/src/ATAC-seq-TCGA-data-analysis"
path.to.save.output <- file.path(path.to.main.src, "panel_7_8_vs_nucleosomeMap")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

all.beds <- Sys.glob(file.path(path.to.main.src, "atacseq_bed", "*.bed"))
names(all.beds) <- to_vec(
  for (item in all.beds){
    str_replace(basename(item), ".bed", "")
  }
)

# read in nucleosome map from healthy samples
nucleosome.map <- read.csv(file.path(path.to.main.src, "rpr_map_EXP0779.sorted.bed"), header = FALSE, sep = "\t")
colnames(nucleosome.map) <- c("chrom", "start", "end", "region", "V5", "strand", "mid.point.start", "mid.point.end")
nucleosome.map.grange <- makeGRangesFromDataFrame(df = nucleosome.map, 
                                                  seqnames.field = "chrom", 
                                                  start.field = "start",
                                                  end.field = "end", 
                                                  keep.extra.columns = TRUE)

for (bed.name in names(all.beds)){
  print(sprintf("Working on bed file %s", bed.name))
  atacseq.beddf <- read.csv(all.beds[[bed.name]], sep = "\t", header = TRUE)
  colnames(atacseq.beddf) <- c("chrom", "start", "end")
  atacseq.beddf$bed.name <- bed.name
  
  atacseq.grange <- makeGRangesFromDataFrame(df = atacseq.beddf, 
                                             seqnames.field = "chrom", 
                                             start.field = "start",
                                             end.field = "end", 
                                             keep.extra.columns = TRUE)
  
  overlap.res <- findOverlapPairs(atacseq.grange, nucleosome.map.grange)
  
  firstdf <- first(overlap.res) %>% data.frame()
  colnames(firstdf) <- to_vec(
    for (item in colnames(firstdf)){
      sprintf("ATACseq_%s", item)
    }
  )
  
  seconddf <- second(overlap.res) %>% data.frame()
  colnames(seconddf) <- to_vec(
    for (item in colnames(seconddf)){
      sprintf("NuMAP_%s", item)
    }
  )
  
  overlapdf <- cbind(firstdf, seconddf) %>%
    rowwise() %>%
    mutate(chrom = ATACseq_seqnames) %>%
    mutate(new.start = min(ATACseq_start, NuMAP_start)) %>%
    mutate(new.end = max(ATACseq_end, NuMAP_end)) %>%
    subset(select = c(chrom, new.start, new.end)) %>%
    mutate(region = sprintf("%s:%s-%s", chrom, new.start, new.end)) %>%
    mutate(region.len = new.end - new.start)
  
  overlapdf <- overlapdf[!duplicated(overlapdf$region),]
  colnames(overlapdf) <- c("chrom", "start", "end", "region", "region.len")
  
  write.table(overlapdf, file.path(path.to.save.output, sprintf("%s_overlap_NucleosomeMap.bed", bed.name)), 
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}


