# gene set enrichment analysis in R

## loading libraries

```r
require(clusterProfiler)
require(enrichR)
require(msigdbr)
require(tidyverse)
```

## loading gene sets

```r
## loading gene sets from msigdbr
dbs <- c(GO = "GO_Biological_Process_2023", HALL = "MSigDB_Hallmark_2020", KEGG = "KEGG_2021_Human", REACTOME = "Reactome_2022")

#' @param genes vector of genes considered - remove gs with <= 5 genes
load_gs <- function(category = NULL, subcategory = NULL, genes=NULL) {
  require(msigdbr)
  gs <- msigdbr(category = category, subcategory = subcategory) %>% 
    dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_description, gs_exact_source)
  
  if (!is.null(genes)) {
    gs <- gs %>%
      group_by(gs_name) %>%
      mutate(gene.in = ensembl_gene %in% genes,
             filter.out = sum(gene.in)) %>%
      filter(filter.out > 5) %>% # more than 5 genes are present in the gene set
      select(-c(gene.in, filter.out))
  }
  
  return(gs)
}
```

## running GSEA

### 1. Using clusterProfiler

documentation: [https://yulab-smu.top/biomedical-knowledge-mining-book/](https://yulab-smu.top/biomedical-knowledge-mining-book/)

tutorial: [https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/](https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/)

```r
## runs clusterProfiler::GSEA
#' runs GSEA on selected gene set
run_msigdb <- function(tmp, ensbl.col="gene", rank.col="avg_log2FC", 
                       selected_set,
                       title="", subtitle=NULL,
                       basename="output/", save.op=NULL, ...){
  
  original_gene_list <- tmp[[rank.col]] 
  names(original_gene_list) <- tmp[[ensbl.col]]
  
  ## GO:BP
  gene_list<-na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  gse <- GSEA(gene_list, TERM2GENE = selected_set[,c("gs_name", "ensembl_gene")], seed = TRUE, verbose = F, minGSSize=10, eps=0)
  
  if (nrow(gse@result)!=0){
    gse@result$Description <- str_to_title(str_replace_all(gse@result$Description, c("^[^_]*_"="", "_"=" ")))
    
    if (!is.null(basename)){
      p <- gsea.dotplot(gse@result, n=10, title=title, savefile=paste0(basename, "_msigdb.tiff"), subtitle=subtitle, ...)
    }
    
    if (!is.null(save.op)){
      x <- as.data.frame(unique(subset(selected_set, select=-c(ensembl_gene, gene_symbol))) )
      gse@result <- left_join(gse@result, x, by= c("ID"="gs_name"))
      
      write.csv(gse@result, file=paste0(save.op, "_msigdb.csv"), row.names = F)}
    
  }
  return(gse)
}
```
