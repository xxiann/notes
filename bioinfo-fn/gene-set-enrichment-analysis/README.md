---
layout:
  title:
    visible: true
  description:
    visible: true
  tableOfContents:
    visible: true
  outline:
    visible: true
  pagination:
    visible: true
---

# gene set enrichment analysis in R

## loading libraries

```r
require(clusterProfiler)
require(enrichR)
require(msigdbr)
require(tidyverse)
```

## loading gene sets

### [GSEA msigdbr](https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp)

Molecular Signatures Database is a collection of annotated gene sets. It contains 8 major collections:&#x20;

* H: hallmark gene sets&#x20;
* C1: positional gene sets&#x20;
* C2: curated gene sets&#x20;
* C3: motif gene sets&#x20;
* C4: computational gene sets&#x20;
* C5: GO gene sets&#x20;
* C6: oncogenic signatures&#x20;
* C7: immunologic signatures

```r
## loading gene sets from msigdbr
dbs <- c(GO = "GO_Biological_Process_2023", HALL = "MSigDB_Hallmark_2020", KEGG = "KEGG_2021_Human", REACTOME = "Reactome_2022")

#' @param genes vector of genes considered - remove gs with <= 5 genes
load_gs <- function(category = NULL, subcategory = NULL, genes=NULL) {
  require(msigdbr)
  gs <- msigdbr(category = category, subcategory = subcategory) %>% 
    dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_description, gs_exact_source)
  
  if (!is.null(genes)) {
    # more than 5 genes are present in the gene set
    # subset genes in the gene set to what is available
    gs <- gs %>%
      group_by(gs_name) %>%
      mutate(gene.in = ensembl_gene %in% genes,
             filter.out = sum(gene.in)) %>%
      filter(filter.out > 5, gene.in) %>% 
      select(-c(gene.in, filter.out))
  }
  
  return(gs)
}
```

other colnames of the db from msigdbr

```
[1] "gs_cat"             "gs_subcat"          "gs_name"            "gene_symbol"        "entrez_gene"
[6] "ensembl_gene"       "human_gene_symbol"  "human_entrez_gene"  "human_ensembl_gene" "gs_id"
[11] "gs_pmid"            "gs_geoid"           "gs_exact_source"    "gs_url"             "gs_description"
```

## running GSEA

Youtube video

* [Gene Set Enrichment Analysis (GSEA) – simply explained!](https://www.youtube.com/watch?v=egO7Lt92gDY)
* [How to interpret GSEA results and plot - simple explanation of ES, NES, leading edge and more!](https://www.youtube.com/watch?v=Yi4d7JIlAsM)

tutorial video

* [Gene Set Enrichment Analysis (GSEA) with fgsea - easy R tutorial](https://youtu.be/itaZcIXZAwg?si=8DYJEy1ypirYLDMa)
* [https://biostatsquid.com/fgsea-tutorial-gsea/](https://biostatsquid.com/fgsea-tutorial-gsea/)

### creating ranked gene list

* ranking of genes determines the top/down as the most 'significant’ genes - ie highest and lowest ranked genes
* 3 ways
  1. p value
  2. log2FC
  3. signed fold change \* -log10pvalue
* highest ranked genes can be associated with the genes that are associated with the specific condition/treatment of the samples and those that are the bottom would be associated with the other condition/treatment

<figure><img src="../../.gitbook/assets/image (2).png" alt="" width="563"><figcaption></figcaption></figure>

* Filtering for significant pathways using FDR/p-adjusted&#x20;
  1. Sort by adjusted p-value < 0.25 (Recommended)
  2.  Filter data with adjusted pvalue and sort by absolute NES value.

      [https://www.biostars.org/p/375584/](https://www.biostars.org/p/375584/)
* NES - is corrected for size of gene sets

### 1. using clusterProfiler

documentation: [https://yulab-smu.top/biomedical-knowledge-mining-book/](https://yulab-smu.top/biomedical-knowledge-mining-book/)

tutorial: [https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/](https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/)

```r
## runs clusterProfiler::GSEA
#' runs GSEA on selected gene set
run_msigdb <- function(tmp, ensbl.col="gene", 
                       rank.col="avg_log2FC", 
                       selected_set_list,
                       pvalueCutoff = 0.25,
                       title="",
                       basename="plot/rnaseq/gsea/msigdb/",
                       save.op=NULL){
  require(clusterProfiler)
  
  ## ranked gene list
  original_gene_list <- tmp[[rank.col]] 
  names(original_gene_list) <- tmp[[ensbl.col]]
  gene_list<-na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  ## running gsea
  list_gse = list()
  for (i in names(selected_set_list)) {
    selected_set = selected_set_list[[i]]
    gse <- GSEA(gene_list, TERM2GENE = selected_set[,c("gs_name", "ensembl_gene")], seed = TRUE, verbose = F, minGSSize=10,  pvalueCutoff = pvalueCutoff)
    
    ##saving
    if (nrow(gse@result)!=0){ #str_to_title
      gse@result$Description <- str_replace_all(gse@result$Description, c("^[^_]*_"="", "_"=" "))
      
      ## plotting
      if (!is.null(basename)){
        gsea.dotplot(gse@result, n=10, title=title, savefile=paste0(basename, "_",i,".tiff"))
      }
      
      ## saving in csv
      if (!is.null(save.op)){
        x <- as.data.frame(unique(subset(selected_set, select=-c(ensembl_gene, gene_symbol))) )
        gse@result <- left_join(gse@result, x, by= c("ID"="gs_name"))
        write.csv(gse@result, file=paste0(save.op, "_",i,".csv"), row.names = F)
        }
    }
    
    list_gse[[i]] = gse
  }
  
  return(list_gse)
}
```



## Over-representation analysis

over-representation analysis ≠ gene expression analysis

* To perform the over-representation analysis, we need a list of background genes and a **list of significant genes**. For our **background dataset** we will use all genes tested for differential expression (all genes in our results table). For our significant gene list we will use genes with p-adjusted values less than 0.05 (we could include a fold change threshold too if we have many DE genes).

### 1. using clusterProfiler

```r
#' @description
#' uses clusterProfiler::enricher - which is a general fn
#' outputs gene ratio and p adj (affected by bgRatio) - ranked by pvalue
#' @param selected_set_list has to be a named list
run_msigdb_enrich <- function(over, under, title = c("", ""), 
                      basename="plot/rnaseq/gsea/msigdb/", selected_set_list,
                       color=c("firebrick","seagreen"), save.op=NULL){
  require(clusterProfiler)
  
  list_over = list()
  list_under = list()
  for (i in names(selected_set_list)) {
    selected_set = selected_set_list[[i]]
    
    ## over-representation
    res <- enricher(over, TERM2GENE=selected_set[,c("gs_name", "ensembl_gene")]) 
    
    if (nrow(res@result)!=0){ #str_to_title
      res@result$Description <- str_replace_all(res@result$Description, c("^[^_]*_"="", "_"=" "))
      
      if (!is.null(basename)){
        enricher.barplot(res, fill=color[2], title=title[2], savefile = paste0(basename,"_",i,"_over.tiff"))
      }
    
      if (!is.null(save.op)){
        x <- as.data.frame(unique(subset(selected_set, select=-c(ensembl_gene, gene_symbol))) )
        res@result <- left_join(res@result, x, by= c("ID"="gs_name"))
        
        write.csv(res@result, file=paste0(save.op, "_",i,"_over.csv"), row.names = F)}
    }
    
    list_over[[i]] = res
  
    ## under-representation
    res <- enricher(under, TERM2GENE=selected_set[,c("gs_name", "ensembl_gene")]) 
    
    if (nrow(res@result)!=0){ #str_to_title
      res@result$Description <- str_replace_all(res@result$Description, c("^[^_]*_"="", "_"=" "))
      
      if (!is.null(basename)){
        enricher.barplot(res, fill=color[2], title=title[2], savefile = paste0(basename,"_",i,"_under.tiff"))
      }
      
      if (!is.null(save.op)){
        res@result <- left_join(res@result, x, by= c("ID"="gs_name"))
        write.csv(res@result, file=paste0(save.op, "_",i,"under.csv"), row.names = F)}
    }
    
    list_under[[i]] = res
  }
  return(list(list_under, list_over))
}
```

### 2. using enrichR

<details>

<summary>output for <code>enrichR::enrichr</code></summary>

* **Odds Ratio (OR)**:
  * The odds ratio is a statistical measure used in enrichment analysis to assess the association between a gene set (e.g., a pathway or functional category) and a specific condition (e.g., differentially expressed genes).
  * Specifically, it quantifies how much more likely the gene set is to be enriched (overrepresented) among the differentially expressed genes compared to what would be expected by chance.
  * An OR greater than 1 indicates enrichment, while an OR less than 1 suggests depletion.
  *   It’s calculated as: where:

      OR=odds of chanceodds of enrichment=b⋅ca⋅d

      * (a) = Number of genes in the gene set that are differentially expressed.
      * (b) = Number of genes not in the gene set that are differentially expressed.
      * (c) = Number of genes in the gene set that are not differentially expressed.
      * (d) = Number of genes not in the gene set that are not differentially expressed.

- **Combined Score**:
  * The combined score represents the overall significance of enrichment for a gene set.
  * It combines multiple statistical metrics (such as p-values, Z-scores, or Fisher’s exact test results) into a single score.
  * Higher combined scores indicate stronger evidence of enrichment.
  * It helps prioritize gene sets based on their biological relevance.
  * The specific formula for calculating the combined score depends on the enrichment method used (e.g., Fisher’s exact test, hypergeometric test, etc.).

</details>

```r
#' @description
#' uses enrichR::enrichr - which accesses online databases
#' outputs combined score = c = log(p) * z; which is used to rank significance
#' uses gene symbol
#' @param selected_set_list has to be a named vector
run_enrich <- function(over, under, title = c("", ""), 
                       basename="plot/rnaseq/gsea/msigdb/",
                       selected_set_list,
                       color=c("firebrick","seagreen"), 
                       save.op=NULL){
  require(enrichR)
  
  list_over = list()
  list_under = list()
  for (i in names(selected_set_list)) {
    selected_set = selected_set_list[[i]]
    ## over-representation
    enrich_results <- enrichr(genes = over, databases = selected_set)[[1]] %>%
      separate_wider_delim(Term, " (", names=c("Term", "id"), too_few="align_start", too_many="merge", cols_remove=T) 
    
    if (!is.null(basename)){
      enrich_results %>%
        filter(Adjusted.P.value<0.05) %>%
        enrichR.barplot(., fill=color[1], title=title[1], savefile = paste0(basename,"_",i,"_over.tiff"))
    }
    
    if (!is.null(save.op)){
      if (nrow(enrich_results)!=0){
        write.csv(enrich_results, file=paste0(save.op, "_",i,"_over.csv"), row.names = F)}
    }
    
    list_over[[i]] = enrich_results
    
    ## under-representation
    enrich_results <- enrichr(genes = under, databases = selected_set)[[1]] %>%
      separate_wider_delim(Term, " (", names=c("Term", "id"), too_few="align_start", too_many="merge", cols_remove=T) 
    
    
    if (!is.null(basename)){
      enrich_results %>%
        filter(Adjusted.P.value<0.05) %>%
        enrichR.barplot(., fill=color[2], title=title[2], savefile = paste0(basename,"_",i,"_under.tiff"))
    }
    
    if (!is.null(save.op)){
      if (nrow(enrich_results)!=0){
        write.csv(enrich_results, file=paste0(save.op, "_",i,"_under.csv"), row.names = F)}
    }
    
    list_under[[i]] = enrich_results
  }
  return(list(list_under, list_over))
}
```





