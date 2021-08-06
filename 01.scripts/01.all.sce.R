
library(Seurat)
library(tidyverse)
library(reticulate)

#---------- load data

sample.dirs <- list.dirs("/home/pauling/projects/08.lung.pd1/01.data")
sample.dirs <- sample.dirs[stringr::str_detect(sample.dirs,"tumor")]

load.data <- function(path){ 
  gene_path <- "/home/pauling/projects/01_data/03_gene/coding_gene.rds.gz"
  cd_gene <- readr::read_rds(gene_path)
  
  matr <- Read10X(path)
  ig.name <- rownames(matr)
  over_gene <- c(intersect(rownames(matr), cd_gene$gene_name))
  matr <- matr[over_gene,]
  total_umi <- colSums(as.matrix(matr))
  gene_count <- colSums(as.matrix(matr) > 0, na.rm = T)
  matr <- matr[, c(total_umi < 25000 & total_umi > 500 & gene_count > 500)]

  return(matr)
}
all.expr <- list()

for (i in 1:length(sample.dirs)) {
  print(sample.dirs[i])
  all.expr[[i]] <- load.data(path = sample.dirs[i])
}

sample.meta <- tibble(
  samples = sample.dirs
) %>%
  tidyr::separate(samples, c(paste0("Col",".",1:7),"sample"), sep = "\\/") %>%
  dplyr::select(sample) %>%
  as.data.frame() %>%
  tidyr::separate(sample, c("patient","tumour","num"), sep = "\\_", remove = F)

for (i in 1:length(all.expr)) {
  colnames(all.expr[[i]]) <- paste0(sample.meta$sample[i],".",colnames(all.expr[[i]]))
}

all.expr <- Reduce(cBind, all.expr)

#---------- meta data

meta <- tibble(cellid = colnames(all.expr)) %>%
  tidyr::separate(cellid, c("patient"), sep = "\\_", remove = F)

meta <- as.data.frame(meta)
rownames(meta) <- meta$cellid

#---------- Clustering

sce <- CreateSeuratObject(counts = all.expr, meta.data = meta)


sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes, do.center = F, do.scale = F)
sce %>% readr::write_rds("/home/pauling/projects/08.lung.pd1/02.all.rds.gz", compress = "gz")

#ent.res <- ROGUE::SE_fun(sce@assays$RNA@counts)

#sce <- RunPCA(sce, features = ent.res$Gene[1:1500])
#sce <- RunUMAP(sce, dims = 1:10)

#sce <- bbknn.batch(sce, r = 1)
#DimPlot(sce, label = T)
