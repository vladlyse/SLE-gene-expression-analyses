# Workflow partially forked from @kpatel427
# Data: 
#   GSE118254 from NCBI GEO,
#   SLE-associated gene list from Bentham et al., Nature Genetics 2015


# prepare the environment
require(dplyr)
require(tidyverse)
require(GEOquery)
require(ggplot2)
require(GGally)
require(M3C)
Sys.setenv('VROOM_CONNECTION_SIZE' = 131072 * 10)
setwd('/home/v/Documents/Single-cell RNA-seq/Projects/')


# prepare the metadata
ge.meta <- getGEO('GSE118254', GSEMatrix = T)
cols.sel <- c('title', 'cell subtype:ch1', 'sorting markers:ch1', 'subject status:ch1')
ge.meta.phe <- pData(phenoData(ge.meta[[1]]))
ge.meta.sel <- 
  select(ge.meta.phe, all_of(cols.sel)) %>% 
  rename(cell_subtype = 'cell subtype:ch1', 
         markers = 'sorting markers:ch1', 
         status = 'subject status:ch1') %>%
  mutate(title = gsub('RNA-seq', '', title)) %>%
  mutate(title = gsub('[[]]', '', title)) %>%
  mutate(title = gsub(' ', '', title)) %>%
  mutate(status = gsub('Systemic lupus erythematosus ', '', status)) %>%
  mutate(status = gsub('[()]', '', status))
  

# prepare the data
ge <- read.csv('data.GSE118254.csv')
ge <- na.omit(ge)
ge <- ge[!(names(ge) %in% c('ENTREZID', 'length'))] %>% rename(gene = SYMBOL)
ge <- ge[!duplicated(ge[c('gene')]), ]
names(ge) <- gsub(names(ge), pattern = '.fpkm', replacement = '')
ge.re <- ge %>% gather(key = 'samples', value = 'fpkm', -gene, na.rm = T)
ge.re <- left_join(ge.re, ge.meta.sel, by = c('samples' = 'title'))


# import genes associated with SLE and filter the gene expression data by them
genes <- scan('data.genes_bentham2015.origsort.txt', what='', sep='\n')
# ge.filter <- ge.re
ge.filter <- filter(ge.re, gene %in% genes)
ge.filter.H <- filter(ge.filter, status == 'Healthy')
ge.filter.L <- filter(ge.filter, status == 'SLE')


# TODO make scatterplot matrix with correlation tests for the imported genes


# make heatmaps
p_heatmap.samples.H <- ggplot(ge.filter.H, aes(x = samples, y = gene, fill = fpkm)) + 
  geom_tile() + scale_fill_gradient(low = 'white', high = '#7c2b83') + 
  labs(title = 'Heatmap of SLE-associated genes in healthy samples') +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
ggsave('plots.heatmap.samples.healthy.png', p_heatmap.samples.H, width = 10, height = 10)

p_heatmap.samples.L <- ggplot(ge.filter.L, aes(x = samples, y = gene, fill = fpkm)) + 
  geom_tile() + scale_fill_gradient(low = 'white', high = '#7c2b83') + 
  labs(title = 'Heatmap of SLE-associated genes in SLE samples') + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
ggsave('plots.heatmap.samples.SLE.png', p_heatmap.samples.L, width = 10, height = 10)

p_heatmap.cells.H <- ggplot(ge.filter.H, aes(x = cell_subtype, y = gene, fill = fpkm)) + 
  geom_tile() + scale_fill_gradient(low = 'white', high = '#7c2b83') + 
  labs(title = 'Heatmap of SLE-associated genes in B cells of healthy samples') + 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5))
ggsave('plots.heatmap.cells.healthy.png', p_heatmap.cells.H, width = 10, height = 10)

p_heatmap.cells.L <- ggplot(ge.filter.L, aes(x = cell_subtype, y = gene, fill = fpkm)) + 
  geom_tile() + scale_fill_gradient(low = 'white', high = '#7c2b83') + 
  labs(title = 'Heatmap of SLE-associated genes in B cells of SLE samples') + 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5))
ggsave('plots.heatmap.cells.SLE.png', p_heatmap.cells.L, width = 10, height = 10)



# visualise the whole gene sets as cell type clusters
ge.index <- ge
rownames(ge.index) <- ge.index$gene
ge.index$gene <- NULL
ge.index <- ge.index[ , order(names(ge.index))]

ge.filter.index <- filter(ge, gene %in% genes)
rownames(ge.filter.index) <- ge.filter.index$gene
ge.filter.index$gene <- NULL
ge.filter.index <- ge.filter.index[ , order(names(ge.filter.index))]


p_pca <- pca(ge.index, labels=ge.meta.sel$cell_subtype, legendtitle = "B cell subtype", legendtextsize = 10, axistextsize = 10, dotsize=2) + 
  labs(title = 'Principal Component Analysis of all available genes in B cells')
ggsave('plots.dim.pca.png', p_pca, width = 10, height = 5)
p_pca.filter <- pca(ge.filter.index, labels=ge.meta.sel$cell_subtype, legendtitle = "B cell subtype", legendtextsize = 10, axistextsize = 10, dotsize=2) + 
  labs(title = 'Principal Component Analysis of SLE-associated genes in B cells')
ggsave('plots.dim.pca.filter.png', p_pca.filter, width = 10, height = 5)

p_umap <- umap(ge.index, labels=as.factor(ge.meta.sel$cell_subtype), legendtitle = "B cell subtype", controlscale=TRUE, scale=3) + 
  labs(title = 'Uniform Manifold Approximation and Projection of all available genes in B cells')
ggsave('plots.dim.umap.png', p_umap, width = 10, height = 5)
p_umap.filter <- umap(ge.filter.index, labels=as.factor(ge.meta.sel$cell_subtype), legendtitle = "B cell subtype", controlscale=TRUE, scale=3) + 
  labs(title = 'Uniform Manifold Approximation and Projection of SLE-associated genes in B cells')
ggsave('plots.dim.umap.filter.png', p_umap.filter, width = 10, height = 5)


# TODO convert FPKM to gene counts
