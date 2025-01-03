# batch correction

references: [https://www.biostars.org/p/266507/](https://www.biostars.org/p/266507/)

Batch effect refers technical variation and differences caused when genomic data are produced in batches due to logistical or practical restrictions and can have unfavorable impact on downstream biological analysis.

* batch effect is usually included in DESeq as a covariate instead
  * [visualising in without batch variation in PCA after VST](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-after-vst-are-there-still-batches-in-the-pca-plot)
* edgeR accepts batch corrected dataframe

## visualising using pca

Reference: \
[https://liulab-dfci.github.io/RIMA/Preprocessing.html#batch-effect-removal](https://liulab-dfci.github.io/RIMA/Preprocessing.html#batch-effect-removal) [https://github.com/liulab-dfci/RIMA\_pipeline/blob/master/src/preprocess/pca.R](https://github.com/liulab-dfci/RIMA_pipeline/blob/master/src/preprocess/pca.R)

```r
pca_plot <- function(exprTable, annot,title, Batch, label.hjust = 1.1, label.vjust = 0) {
  batch_n <- length(unique(as.character(annot[colnames(exprTable),Batch])))
  
  print(paste("there are ", batch_n, " batches in your data"))
  
  df <- cbind.data.frame(t(exprTable),batch = as.character(annot[colnames(exprTable),Batch]))
  pca_plot <- autoplot(prcomp(t(exprTable)), data = df, col = 'batch', size = 5, label = TRUE, label.label = colnames(exprTable), label.hjust = label.hjust, label.vjust = label.vjust)+ 
    labs(title=title)+ 
    scale_color_manual(values = brewer.pal(name = "Set1", n = 9)[1:batch_n])+
    scale_x_continuous(expand = c(0.1,0.1)) +
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=12,face = "bold",hjust=1),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold")) +
    guides(color = guide_legend(override.aes = aes(fill = NA, label = "")))
  
  return(pca_plot)
}
```

## sva - [combat-seq](https://github.com/zhangyuqing/ComBat-seq)

* reference
  * [https://bookdown.org/jean\_souza/PreProcSEQ/batch-effect-removal.html#removing-the-batch-effect-of-normalized-counts](https://bookdown.org/jean_souza/PreProcSEQ/batch-effect-removal.html#removing-the-batch-effect-of-normalized-counts)
  * [https://rnabio.org/module-03-expression/0003/06/02/Batch-Correction/](https://rnabio.org/module-03-expression/0003/06/02/Batch-Correction/)

<details>

<summary><code>Combat_seq</code> takes in raw count matrix from genomic studies</summary>

* `counts`. This is your matrix of gene expression read counts (raw counts). Each row is a gene, each column is a sample, and each cell has an integer count for the number of RNA-seq counts observed for that gene/sample combination. In R we want this data to be passed into ComBat-Seq in matrix format (use `as.matrix()` if neccessary).

- `batch`. This is a vector describing the batches you are concerned about. For example, if you have eight samples that you created RNA-seq data for, but for the first four you used library kit (A) and for the last four samples you used library kit (B). In this situation you would define your `batch` vector as: `c(1,1,1,1,2,2,2,2)`.

* `group = NULL`. This is a vector describing your biological condition of interest. For example, if your experiment involved _pairs_ of drug treated and untreated cells, and you did 4 biological replicates. You would define your `group` vector as: c(1,2,1,2,1,2,1,2).

- `covar_mod = NULL`. If you have multiple biological conditions of interest, you can define these with `covar_mod` (covariates) instead of `group`. For example, lets assume we have the same experiment as described above, except that we did four replicates (treated vs untreated pairs), but we alternated use of male and female cells for each of the replicates. You then would define a covariate matrix to supply to `covar_mod` as follows:

* `full_mod = TRUE`. If TRUE include the condition of interest in model. Generally we believe this should be set to the default TRUE. We have yet to find a cohesive explanation for a situation where one would want this to be FALSE.

- `shrink = FALSE`. Whether to apply shrinkage on parameter estimation.

* `shrink.disp = FALSE`. Whether to apply shrinkage on dispersion.

- `gene.subset.n = NULL`. Number of genes to use in empirical Bayes estimation, only useful when shrink = TRUE.

</details>

```r
suppressMessages(library(sva))
library(SummarizedExperiment)

## for ngs RNAseq data
# include condition (group variable)
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group, full_mod=TRUE)

# do not include condition
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=NULL, full_mod=FALSE)
```

`ComBat` - \[array data only] the input data are assumed to be cleaned and normalized before batch effect removal.

```r
## for array data only
batchremoval.combat <- function(mat, model, meta, batch, output=""){
  meta$Batch = meta[,batch]
  mat <- ComBat(dat = mat, 
              batch = meta$Batch, 
              mod = model)
              
  mat <- cbind(GeneID=rownames(mat),mat)
  write.csv(mat,paste0(output,"_tpm_matrix_cb.csv"), row.names = FALSE)
}
```



## limma - removeBatchEffects

[https://github.com/liulab-dfci/RIMA\_pipeline/blob/master/src/preprocess/batch\_removal.R](https://github.com/liulab-dfci/RIMA_pipeline/blob/master/src/preprocess/batch_removal.R)

`removeBatchEffect` takes in a numeric matrix of **log-expression values** and returns it with with batch and covariate effects removed

```r
suppressMessages(library(limma))

batchremoval.limma <- function(expr.dat, meta, design, Condition, output=""){
  samples <- subset(meta, !is.na(meta[,Condition])) # removes NA
  print(samples)
  
  expr.dat <- log2(expr.dat + 1)
  expr.dat <- expr.dat[,rownames(samples)]
  
  ## filtering out genes with low variance among samples  
  CVFILTER <- 0
  mean_nolym <- apply(expr.dat,1,mean)
  var_nolym <- apply(expr.dat,1,var)
  cv_nolym <- abs(var_nolym/mean_nolym)
  filt_genes <- subset(cv_nolym, cv_nolym > CVFILTER)
  
  ## Select those genes that pass variance filtering
  exprZero <- expr.dat
  expr.dat <- expr.dat[rownames(expr.dat) %in% names(filt_genes),]
  exprZero <- subset(exprZero, !(rownames(exprZero) %in% names(filt_genes)))
  
  
  ## assume no batch
  print('No batches used !')
  expr.limma <- expr.dat
  expr.limma = rbind(expr.limma,exprZero)
  
  expr.limma <- cbind(GeneID=rownames(expr.limma),expr.limma)
  write.csv(expr.limma, paste0(output,"_log2tpm_matrix.csv"), row.names = FALSE)
  
  ## assume batch
  samples$Batch <- samples[,design] 
  
  print('Running limma for batch removal')
  expr.limma = tryCatch(
  removeBatchEffect(as.matrix(expr.dat),
                      samples$Batch),
  error = function(e){
  print(e)
  })
  expr.limma = rbind(expr.limma,exprZero)
  expr.limma <- cbind(GeneID=rownames(expr.limma),expr.limma)
  write.csv(expr.limma, paste0(output,"_log2tpm_matrix_limma.csv"), row.names = FALSE)
}
```

