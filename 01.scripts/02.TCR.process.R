
library(Seurat)
library(tidyverse)
library(reticulate)

#---------- load data

sample.dirs <- list.dirs("/home/pauling/projects/08.lung.pd1/01.data/01.vdj")
sample.dirs <- sample.dirs[stringr::str_detect(sample.dirs,"tumor")]

all.vdj <- list()

for (i in 1:length(sample.dirs)) {
  print(sample.dirs[i])
  all.vdj[[i]] <- readr::read_csv(file.path(sample.dirs[i],"filtered_contig_annotations.csv"))
}

sample.meta <- tibble(
  samples = sample.dirs
) %>%
  tidyr::separate(samples, c(paste0("Col",".",1:8),"sample"), sep = "\\/") %>%
  dplyr::select(sample) %>%
  as.data.frame() %>%
  tidyr::separate(sample, c("patient","tumour","num"), sep = "\\_", remove = F) %>%
  dplyr::mutate(sample = stringr::str_remove(sample,".vdj"))

for (i in 1:length(all.vdj)) {
  print(i)
  all.vdj[[i]] <- process.config(all.vdj[[i]])
  all.vdj[[i]] <- process.tcr(all.vdj[[i]])
  all.vdj[[i]] <- all.vdj[[i]] %>% 
    dplyr::mutate(patient = sample.meta$patient[i]) %>% 
    dplyr::mutate(cellid = paste0(sample.meta$sample[i], ".", all.vdj[[i]]$CellName)) %>%
    dplyr::filter(!is.na(`CDR3(Alpha1)`)) %>%
    dplyr::filter(!is.na(`CDR3(Beta1)`))
}

tcr.data <- Reduce(rbind, all.vdj)
tcr.data %>% readr::write_rds("/home/pauling/projects/08.lung.pd1/03.all.tcr.rds.gz", compress = "gz")