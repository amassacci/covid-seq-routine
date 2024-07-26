args <- commandArgs(trailingOnly = T)
print(args)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

run_name <- gsub("run_","",args[1])
outfolder <- paste0("/data/sequence/covid-19/covid-seq/run_",run_name,"/")

datalist = list()
vocs <- read.csv("/data/sequence/covid-19/covid-seq/voc.bed", sep="\t", header = FALSE, col.names = c("Chrom","Start","End","Variant","Lineage"))
coverage.files <- list.files(pattern = "_coverage_vocs.tsv", path = paste0(outfolder,"./QC"), full.names = T)
for (file in coverage.files) {
  data <- read.csv(file, sep="\t", header = FALSE, col.names = c("Chrom","Start","End","Coverage"))
  data <- setDT(data)[, .ind:= cumsum(c(TRUE,Start[-1]!=End[-.N]))] %>%
    group_by(.ind) %>%
    summarize(Coverage = round(mean(Coverage)),
              Start = min(Start),
              End = max(End))
  join <- merge(vocs,data,by=c("Start","End")) %>%
    select(Lineage,Variant,Start,End,Coverage)
  join <- mutate(join, "Sample" = gsub("_coverage_vocs.tsv", "", basename(file), fixed = TRUE))
  datalist[[basename(file)]] <- join
}

big_data = do.call(rbind, datalist)
rownames(big_data) <- NULL

wider <- select(big_data, Sample,Variant,Coverage) %>%
  pivot_wider(names_from = Variant, values_from = Coverage)
write.table(wider, paste0(outfolder,"all_samples_",run_name,"_QC_vocs.tsv"), sep = "\t", quote = FALSE, row.names = F, col.names = T)
