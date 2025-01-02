---
description: 'last update: 9 sept 24'
---

# differential gene expression analysis

## setting env

```r
output="--"
input.dir = "--"
```

## loading libraries

```r
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(ggpubr)
library(rstatix)
library(ggprism)
library(ggfortify)
library(RColorBrewer)
```

## genome annotation file (GENCODE v44)

* obtaining gene id and symbol from .gtf file

```r
# R.utils::gunzip("F://annotation/gencode.v44.primary_assembly.annotation.gtf.gz")
x <- rtracklayer::import("F://annotation/gencode.v44.primary_assembly.annotation.gtf")
exons.list.per.gene <- split(x[x$type=="exon",], x[x$type=="exon",]$gene_id) 
exonic.gene.sizes <- data.frame(size=sum(GenomicRanges::width(GenomicRanges::reduce(exons.list.per.gene))))

gene <- as.data.frame(exons.list.per.gene@unlistData@elementMetadata@listData[c("gene_id", "gene_name", "gene_type")])
gene$seqlvl <- as.numeric(exons.list.per.gene@unlistData@seqnames)
gene <- unique(gene) 
gene$index <- seq(1:nrow(gene))
gene <- cbind(gene, exonic.gene.sizes[gene$gene_id,"size"]) 
colnames(gene)[6] = "size"

#62754 -> 47496
gene <- gene %>% 
  separate_wider_delim(cols = gene_id, delim = ".", names = c("ensembl_id", "version"), cols_remove = FALSE) %>%
  #gene[gene$seqlvl %in% seq(1,25),] %>%
  filter(!grepl("pseudogene", gene_type),
         !grepl("_PAR_Y", gene_id),
         !gene_type %in% c("artifact")) %>%
  arrange(gene_id)
# write.csv(gene, paste0(output, "filtered_gene_rnaseq.csv"), row.names = F)
```

## preparing count table from STAR outputs

column 1: gene ID

column 2: counts for unstranded RNA-seq

column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)

column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

* selecting the 2nd row

```r
rawfile <- list.files(path = "data/", pattern = ".out.tab", full.names = T)
filelist <- list()
filelist[["gene_id"]] <- read.csv(rawfile[1], sep = "\t")[4:62757, 1]
names  = gsub("data/", "", gsub("ReadsPerGene.out.tab", "", gsub("-","_", rawfile)))

for (file in 1:length(names)){
  x <- read.csv(rawfile[file], sep = "\t")
  filelist[[names[file]]] = x[4:62757,2]
}
raw_count = do.call("cbind", filelist) 

write.csv(raw_count, paste0(output, "data/count_matrix.csv"), row.names = F)
```

## setting up deseq

### getting the count matrix

```r
folder = "filtered"
input.file = paste0(input.dir, "all_featurecounts.tsv")

## importing count matrix + data cleaning
df <- read.csv(input.file, sep="\t")
rownames(df) <- df$Geneid
count_matrix <- df[,7:ncol(df)]

colnames(count_matrix) <- c("Direct_A1", "Direct_A2", "Direct_A3", "Mono_A1", "Mono_A2", "Mono_A3", "Mono_D1", "Mono_D2", "Mono_D3", "Trans_A1", "Trans_A2", "Trans_A3")
names=colnames(count_matrix)
```

```r
gene <- read.csv("../rna_seq/filtered_gene_rnaseq.csv")
rownames(gene) <- gene$gene_id

count_matrix <- read.csv("data/count_matrix(filtered).csv") %>%
  dplyr::rename(gene_id = X) %>%
  filter(gene_id %in% gene$gene_id) %>%
  arrange(gene_id)
rownames(count_matrix) <- count_matrix$gene_id
count_matrix <- count_matrix[,2:ncol(count_matrix)]
```

### filtering the count matrix

```r
## removing low variance counts and lowly expressed counts
mean_nolym <- apply(count_matrix,1,mean)
var_nolym <- apply(count_matrix,1,var)
cv_nolym <- abs(var_nolym/mean_nolym)

meanexp <- edgeR::filterByExpr(count_matrix, min.count = 1, min.prop = 2/ncol(count_matrix), min.total.count=1)
filt_genes <- subset(cv_nolym, subset = cv_nolym > 0 & meanexp)

count_matrix <- count_matrix[names(filt_genes),]
gene <- gene[names(filt_genes),] # left 15011 genes

# loading sample info matrix
sampleinfo <- data.frame(ID=names, treatment=factor(c(rep("Direct_A",3), rep("Mono_A", 3), rep("Mono_D", 3), rep("Trans_A", 3)), levels=c("Mono_D", "Mono_A", "Direct_A", "Trans_A")))
rownames(sampleinfo) <- sampleinfo$ID


# tpm normalisation
tpm_matrix <- counts_to_tpm(count_matrix, gene$size)
tpm_matrix <- cbind(as.matrix(gene[,'gene_name']), tpm_matrix)
colnames(tpm_matrix)[1] = "Gene_Symbol"
write.csv(tpm_matrix, file = "data/all_tpm_counts.csv")

```

### PCA

* to visualise the data: identify patterns, clusters, outliers
* to identify if batch effects are confounding the differences

```r
p <- pca_plot(exprTable = count_matrix, annot=sampleinfo, Batch="treatment", title = "PCA before Analysis")
print(p)

ggsave(paste0("plot/",folder,"PCA_before_analysis.tiff"), compression = "lzw", width = 15, height = 12, dpi = 300, units = "cm")
```

<figure><img src="../.gitbook/assets/image (18).png" alt="" width="563"><figcaption></figcaption></figure>

## params

* "Direct\_A1", "Direct\_A2", "Direct\_A3", "Mono\_A1", "Mono\_A2", "Mono\_A3", "Mono\_D1", "Mono\_D2", "Mono\_D3", "Trans\_A1", "Trans\_A2", "Trans\_A3"

```{r}
folder = "filtered/well/" # "well/" # "treatment/"
name = "all"

x <- count_matrix[,c(1:6)]
y <- sampleinfo[c(1:6),]
```

## running DESeq

```{r}
des <- DESeqDataSetFromMatrix(countData = x,
                              colData = y,
                              design = ~ treatment)

nrow(des) # no of genes = 16874

## fit statistical model
des <- DESeq(des)

## saving normalized (median of ratios method of normalization)
normalized_counts <- counts(des, normalized=TRUE)
# write.csv(normalized_counts, paste0("result/",folder,"norm_counts.csv"), row.names = F)
```

### VarianceStabilizingTransformation + Visualisation

#### Between-sample distance matrix

```{r}
vsd <- vst(des)

# calculate between-sample distance matrix
sampleDistMatrix <- as.matrix(dist(t(assay(vsd))))

tiff(paste0("plot/",folder,"hm_vst_",name,".tiff"), compression = "lzw", width = 15, height = 15, res = 300, units = "cm")
pheatmap(sampleDistMatrix)
dev.off()

pheatmap(sampleDistMatrix)
```

<figure><img src="../.gitbook/assets/image (20).png" alt="" width="563"><figcaption></figcaption></figure>

#### PCA after VST

```{r}
plotPCA(object = vsd, intgroup = c("treatment")) + 
  geom_text(aes(label = name), position = position_nudge(y = 0.5)) + 
  theme_bw() + ggtitle("After VarianceStabilizingTransformation")

ggsave(paste0("plot/",folder,"PCA_vst_",name,".tiff"), compression = "lzw", width = 15, height = 12, dpi = 300, units = "cm")
```

<figure><img src="../.gitbook/assets/image (21).png" alt="" width="563"><figcaption></figcaption></figure>

### Differential Gene Analysis

* Table for differentially expressed genes with |log2FoldChange| > 1; padj < 0.0001
* 1243

```r
print(resultsNames(des))
comparison="treatment_Direct_A_vs_Mono_A" 
```

```r
de <- results(object = des, name=comparison, pAdjustMethod = "fdr", alpha = 0.05) # default is 0.1 for padj
de_shrink <- lfcShrink(dds = des,
                       coef=comparison,
                       type="apeglm")

de_shrink <- data.frame(de_shrink@listData)

de_shrink$GeneID <- rownames(de_shrink)
de_shrink <- left_join(de_shrink, gene, by=c("GeneID"="gene_id")) %>%
  rename(Symbol = gene_name) %>%
  arrange(padj)

write.csv(de_shrink, paste0("result/",folder,"deseq2_",comparison,".csv"), row.names = F)

# top significant DEGs
topdeg <- de_shrink %>%
  filter(log2FoldChange > 0.58 | log2FoldChange < -0.58) %>%
  filter(padj < 0.05) %>%
  arrange(padj)

write.csv(de_shrink, paste0("result/",folder,"deseq2_",comparison,"_sig_0.58_0.05.csv"), row.names = F)

```

## plots

### volcano plot

* Table for differentially expressed genes with |log2FoldChange| > 1; padj < 0.01

```r
folder= "filtered/treatment/"
title="Trans AraC vs Mono DMSO"
comparison="treatment_Trans_A_vs_Mono_D" 
de_shrink <- read.csv(paste0("result/",folder,"deseq2_",comparison,".csv"))

p1 <- volc.plot(de_shrink, paste("Differential Gene Expression:",title,"(p=0.0001)"), fc=2, sig=0.0001, labelCol = "Symbol")
ggsave(paste0("plot/",folder,"volc_plot_",comparison,".tiff"),plot=p1, width = 20, height = 12, units = "cm", dpi=400)

show(p1)
```

<figure><img src="../.gitbook/assets/image (22).png" alt="" width="563"><figcaption></figcaption></figure>

### gsea

```r
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
set.seed(10)

# loading gene set
hall_set <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, ensembl_gene)
go_set <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, ensembl_gene) %>%
  filter(grepl("GOBP", gs_name))
```

```r
# obtain gene_list
original_gene_list <- de_shrink$log2FoldChange
names(original_gene_list) <- de_shrink$ensembl_id

gene_list<-na.omit(original_gene_list) # omit any NA values 
gene_list = sort(gene_list, decreasing = TRUE) # sort the list in decreasing order

# running gsea
em2 <- GSEA(gene_list, TERM2GENE = hall_set, seed = TRUE, verbose = F, minGSSize=20, eps=0)

# save
write.csv(em2@result, paste0("result/",folder,"HM_",comparison,".csv"), row.names = F)

# plotting
em2@result %>%
  mutate(Description=gsub("_"," ", gsub("HALLMARK_","", Description))) %>%
  dot.plot(., title="MSigDB Hallmark Gene Set", subtitle = title,
           savefile=paste0("plot/",folder,"HM_NES_plot_",comparison,".tiff"))
```

<figure><img src="../.gitbook/assets/image (23).png" alt="" width="563"><figcaption></figcaption></figure>

```{r}
# running gsea
em2 <- GSEA(gene_list, TERM2GENE = go_set, seed = TRUE, verbose = F, minGSSize=20, eps=0)

# save
write.csv(em2@result, paste0("result/",folder,"GOBP_",comparison,".csv"), row.names = F)

em2@result %>%
  mutate(Description=gsub("_"," ", gsub("GOBP_","", Description))) %>%
  dot.plot(., title="MSigDB GO BP Gene Set", subtitle = title,
           savefile=paste0("plot/",folder,"GOBP_NES_plot_",comparison,".tiff"))

```

<figure><img src="../.gitbook/assets/image (24).png" alt="" width="563"><figcaption></figcaption></figure>

## obtaining the LFC of core enrichment gene (DEGs) of selected pathways

```r
folder= "filtered/well/"
# title="Direct AraC vs Mono AraC"
comparison="treatment_Direct_A_vs_Mono_A" 
de_shrink <- read.csv(paste0("result/",folder,"deseq2_",comparison,".csv"))
result <- read.csv(paste0("result/",folder,"GOBP_",comparison,".csv"), row.names = 1)

selected = c("GOBP_CELL_KILLING", "GOBP_POSITIVE_REGULATION_OF_CELL_KILLING", "GOBP_REGULATION_OF_CELL_KILLING")
selected = result[selected,c("Description", "core_enrichment")]

for (i in 1:nrow(selected)) {
  gene <- strsplit(selected[i,2],split="/")[[1]]
  df = de_shrink[de_shrink$ensembl_id %in% gene,]
  write.csv(df, paste0("result/",folder,selected[i,1], "_",comparison,".csv"), row.names = F)
}
```

```
```

