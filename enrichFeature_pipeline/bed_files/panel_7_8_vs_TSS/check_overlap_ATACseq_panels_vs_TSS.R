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
path.to.save.output <- file.path(path.to.main.src, "panel_7_8_vs_TSS")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

all.beds <- Sys.glob(file.path(path.to.main.src, "atacseq_bed", "*.bed"))
names(all.beds) <- to_vec(
  for (item in all.beds){
    str_replace(basename(item), ".bed", "")
  }
)

# read in TSS map from healthy samples
tss.dir <- "/media/hieunguyen/HNSD01/src/ATAC-seq-TCGA-data-analysis/tss_beds"

for (tss.source in c("biomart", "UCSC")){
  up.flanking.size <- 1000
  down.flanking.size <- 1000
  
  tss.input.file <- sprintf(file.path(tss.dir, sprintf("%s/promoter_regions_up_%s_down_%s.bed", tss.source, up.flanking.size, down.flanking.size)))
  tssdf <- read.csv(tss.input.file, sep = "\t", header =  FALSE)
  colnames(tssdf) <- c("chrom", "start", "end", "region")
  
  tss.grange <- makeGRangesFromDataFrame( df = tssdf, 
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
    
    overlap.res <- findOverlapPairs(atacseq.grange, tss.grange)
    
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
    
    write.table(overlapdf, file.path(path.to.save.output, sprintf("%s_overlap_TSS_up_%s_down_%s_%s.bed", bed.name, up.flanking.size, down.flanking.size, tss.source)), 
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}


