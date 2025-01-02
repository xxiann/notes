# bulk RNAseq in R

## loading libraries

```r
require(tidyverse)
```

## GTF files to CSV function

to obtain gene list containing gene\_id, gene\_name, gene\_biotype, seqlvl, size as reference

data format of GENCODE: [https://www.gencodegenes.org/pages/data\_format.html](https://www.gencodegenes.org/pages/data_format.html)

```r
## coded for GENCODE gtf files
gtf_to_csv <- function(gtffile, savefile){
  x <- rtracklayer::import(gtffile)
  exons.list.per.gene <- split(x[x$type=="exon",], x[x$type=="exon",]$gene_id) 
  exonic.gene.sizes <- data.frame(size=sum(GenomicRanges::width(GenomicRanges::reduce(exons.list.per.gene))))

  gene <- as.data.frame(exons.list.per.gene@unlistData@elementMetadata@listData[c("gene_id", "gene_name", "gene_type")])
  gene$seqlvl <- as.numeric(exons.list.per.gene@unlistData@seqnames)
  gene <- unique(gene) 
  gene$index <- seq(1:nrow(gene))
  gene <- cbind(gene, exonic.gene.sizes[gene$gene_id,"size"]) 
  colnames(gene)[6] = "size"

  #removing pseudogenes, artifacts and pseudoautosomal regions
  gene <- gene %>% 
    separate_wider_delim(cols = gene_id, delim = ".", names = c("ensembl_id", "version"), cols_remove = FALSE) %>%
    #gene[gene$seqlvl %in% seq(1,25),] %>% # limiting to primary chromosomes
    filter(!grepl("pseudogene", gene_type),
          !grepl("_PAR_Y", gene_id),
          !gene_type %in% c("artifact")) %>%
    arrange(gene_id)
  
  write.csv(gene, savefile, row.names = F)
  return(gene)
}
```

## gene count normalisation

reference: [DGE\_workshop-hbctraining](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)

normalization is essential for exploratory data analysis, visualization of data, and whenever you are exploring or comparing counts between or within samples such as in differential expression analysis



## counts to TPM

```r
counts_to_tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
```

