# gene set enrichment analysis

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

