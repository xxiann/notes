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

\*\*\*\*\*\*\* TBC with details of normalisation methods

## counts to TPM

```r
counts_to_tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
```



## gene symbols, gene ids

* known issues of excel with symbols
  * Gene symbols like "DEC1" or "MAR1" are automatically converted to dates (e.g., "1-Dec" or "1-Mar"), causing data corruption especially for files stored in excel workbook
  * gene duplication like MARC1 and MARCH1
  * attempt in correcting such errors - [https://kuanrongchan-date-to-gene-converter-date-gene-tool-oodr7h.streamlit.app/](https://kuanrongchan-date-to-gene-converter-date-gene-tool-oodr7h.streamlit.app/)
* standardisation of gene symbols - old versions
  * HGNC - HUGO gene name nomenclature committee&#x20;
  * gene name file - [https://www.genenames.org/download/](https://www.genenames.org/download/)
  * multi symbol checker - [https://www.genenames.org/tools/multi-symbol-checker/](https://www.genenames.org/tools/multi-symbol-checker/)

```r
## finding indexes for gene symbol; used when no gene_id can be obtained
# Alias symbol - the input has been added by an HGNC editor as aa alias of the approved symbol.
# Previous symbol - the input was previously an approved gene symbol, but the gene has since been updated with the approved symbol shown.
gene_hgnc <- read.csv("../hgnc_complete_set.txt", sep="\t")
hgnc.index <- function(df.gene, hgnc.gene){
  require(data.table)
  df <- data.table(index=character(), og=character())
  
  find_name <- toupper(df.gene)
  temp2 <- gsub("\\|$", "", paste(hgnc.gene$alias_symbol, hgnc.gene$prev_symbol, sep="|"))
  temp2 <- strsplit(toupper(temp2),"\\|") # alias & prev symbol
  temp3 <- toupper(hgnc.gene$symbol) # standardised symbol
  
  pb <- txtProgressBar(min = 0, max = length(find_name), style = 3, width = 50, char = "=")
  start_time <- Sys.time()
  for (i in seq_along(find_name)){
    index <- match(find_name[i], temp3) # checks the standardised symbol
    if (is.na(index)) { 
      index <- which(sapply(temp2, function(x) find_name[i] %in% x)) # checks if in alias+prvs symbol
    }
    df <- rbind(df, data.table(index=paste(index, collapse=","), og=find_name[i]))
    setTxtProgressBar(pb, i)
  }
  end_time <- Sys.time()
  close(pb)
  
  time_taken <- end_time - start_time
  print(paste("Time taken: ", time_taken))
  
  # Check for rows with multiple detected outputs and remove duplicates for those that alr hv exact match
  df[, index := sapply(index, function(x) {
    if (grepl(",", x)) {
      indices <- unlist(strsplit(x, ","))
      unique_indices <- unique(indices)
      # Remove duplicates seen in other rows
      unique_indices <- unique_indices[!unique_indices %in% unlist(df$index[!grepl(",", df$index)], ",")]
      paste(unique_indices, collapse = ",")
    } else {
      x
    }
  })]
  
  return(df)
}
```



### biomaRt - ensembl database

[https://asia.ensembl.org/info/data/biomart/biomart\_r\_package.html](https://asia.ensembl.org/info/data/biomart/biomart_r_package.html)

```r
library(biomaRt)
ensembl75 = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=75)

gene_description <- getBM(attributes=c('ensembl_gene_id','external_gene_id','description'), 
    filters = 'ensembl_gene_id', values = gene$gene_id, mart =ensembl75) 
```



