
library(Seurat)
library(tidyverse)
library(reticulate)
library(ks)
library(RColorBrewer)

#---------- load data -----------#

sce <- readr::read_rds("/raid1/pauling/projects/08.lung.pd1/02.all.rds.gz")
all.tcr  <- readr::read_rds("/raid1/pauling/projects/08.lung.pd1/03.all.tcr.rds.gz")

sce.filt <- sce[,colnames(sce) %in% all.tcr$cellid]
all.tcr <- all.tcr %>% dplyr::filter(cellid %in% colnames(sce.filt))

sce.filt <- NormalizeData(sce.filt, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(sce.filt)
sce.filt <- ScaleData(sce.filt, features = all.genes, do.center = F, do.scale = F)

#---------- Tex clones ----------#

all.tcr <- all.tcr %>%
  dplyr::mutate(clone.id = paste(`Identifier(Alpha1)`,`Identifier(Beta1)`, sep = ":")) %>%
  as.tibble()

use.clones <- unique(all.tcr$clone.id)

use.clones.cells <- all.tcr %>%
  dplyr::select(patient, clone.id, cellid) %>%
  dplyr::group_by(patient, clone.id) %>%
  tidyr::nest() %>%
  dplyr::mutate(num = purrr::map_dbl(data, nrow)) %>%
  dplyr::filter(num >= 1) %>%
  dplyr::mutate(cloneid = paste0("clone.",1:nrow(.)))

mcells <- unique(use.clones.cells$cloneid)

cal.genes <- c("CXCL13","CD8A")
meta.expr <- Matrix::Matrix(data = 0, nrow = length(mcells), ncol = length(cal.genes))

for (i in 1:length(mcells)) {
  tmp.cell <- use.clones.cells$data[[i]]$cellid
  print(i)
  if (length(tmp.cell) == 1) {
    meta.expr[i,] <- sce.filt@assays$RNA@scale.data[cal.genes,tmp.cell]
  }else if(length(tmp.cell) > 1){
    meta.expr[i,] <- Matrix::rowMeans(sce.filt@assays$RNA@scale.data[cal.genes,tmp.cell])
  }
}

pr.matr <- as.matrix(meta.expr)
rownames(pr.matr) <- mcells
colnames(pr.matr) <- cal.genes

tibble(
  CD8A = pr.matr[,"CD8A"],
  CXCL13 = pr.matr[,"CXCL13"],
  size = use.clones.cells$num,
  patient = use.clones.cells$patient,
  clone.id = use.clones.cells$clone.id
) -> pda

pda <- as.data.frame(pda)
rownames(pda) <- rownames(pr.matr)

pda %>% readr::write_rds("/home/pauling/projects/08.lung.pd1/04.clone.tex.info.rds.gz", compress = "gz")

#------------- Tex cell frequency -------------#

tex.clones.pro <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD8A > 1 & CXCL13 > 0, "Tex", "Other"))

response = c('non-MPR','MPR','MPR','non-MPR','non-MPR','MPR','non-MPR','MPR','non-MPR','non-MPR','non-MPR','non-MPR','non-MPR','MPR','MPR')

used.color <- c("#FA8072", "#6CA6CD")

tex.clones.pro %>%
  dplyr::group_by(patient, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(response = response) -> pda

pda %>%
  ggplot(aes(response, freq)) +
  geom_boxplot(fill = c(used.color), outlier.colour = "black", outlier.shape = 1, width = 0.6) +
  #geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Frequency of Tex/tumour-reactive\nclones in T cells"
  ) -> p_boxplot

dat <- ggplot_build(p_boxplot)$data[[1]]
dat$xsp <- 1/2 * (dat$xmax + dat$xmin) - 1/4 * (dat$xmax - dat$xmin)
dat$xep <- 1/2 * (dat$xmax + dat$xmin) + 1/4 * (dat$xmax - dat$xmin)
p_boxplot +
  geom_segment(data = dat,  aes(x = xmin, xend = xmax, y = middle, yend = middle), colour = "black", size = 0.6) +
  geom_segment(data = dat,  aes(x = xsp, xend = xep, y = ymin, yend = ymin), colour = "black", size = 0.6) +
  geom_segment(data = dat,  aes(x = xsp, xend = xep, y = ymax, yend = ymax), colour = "black", size = 0.6)

wilcox.test(pda$freq[pda$response == "MPR"], pda$freq[pda$response != "MPR"])

#-------------- Waterfall plot -----------------#

Tex.prop <- readr::read_rds("/raid1/pauling/projects/08.lung.pd1/04.clone.tex.info.rds.gz")

tex.clones.pro <- Tex.prop %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD8A > 1 & CXCL13 > 0, "Tex", "Other"))

response = c('non-MPR','MPR','MPR','non-MPR','non-MPR','MPR','non-MPR','MPR','non-MPR','non-MPR','non-MPR','non-MPR','non-MPR','MPR','MPR')

tex.clones.pro %>%
  dplyr::group_by(patient, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(response = response) %>% 
  dplyr::arrange(desc(freq)) -> pda

used.color <- ifelse(pda$response == "MPR", used.color[1],used.color[2])

barplot(pda$freq-0.4,
        col=used.color,
        border=used.color,
        space=0.5,
        ylim=c(-0.3,0.2),
        cex.axis=1)
