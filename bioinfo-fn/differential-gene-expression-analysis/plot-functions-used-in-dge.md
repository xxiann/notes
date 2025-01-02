# plot functions used in dge

```r
###### pca plot #########
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

###### volcano plot #########
volc.plot <- function(data, title, fc=0.58, sig=5e-2, labelCol="Symbol") {
  data <- as.data.frame(data)
  data$labelCol <- data[,labelCol]
  sigDEG <- data %>%
    dplyr::mutate(Significance = ifelse((log2FoldChange > fc & padj < sig), "Up-regulated", ifelse((log2FoldChange < -fc & padj < sig), "Down-regulated", "Not Significant")))
  
  sigDEG <- sigDEG %>%
    dplyr::mutate(label = ifelse(Significance %in% c("Up-regulated","Down-regulated"), labelCol, NA))
  
  sigDEG$Significance <- factor(sigDEG$Significance, levels = c("Up-regulated","Down-regulated", "Not Significant")) 
  
  max <- max(abs(sigDEG$log2FoldChange))
  
  vplot <- ggplot(sigDEG) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = Significance, alpha=0.3)) +
    geom_vline(xintercept = c(-fc,fc), linetype = "dashed") +
    geom_hline(yintercept = -log10(sig), linetype = "dashed") +
    ggrepel::geom_label_repel(aes(x = log2FoldChange, y = -log10(padj), label = label),box.padding = 0.35, point.padding = 0.5, segment.color = "grey50", min.segment.length = 0.1, max.overlaps = 20, na.rm = TRUE) +
    scale_colour_manual(values = c("red", "blue", "grey")) +
    ggtitle(title) + 
    xlab("Log2 Fold Change") + 
    ylab("-Log10 Adjusted p-value") +
    xlim(-max, max) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10), 
          axis.title=element_text(size=12),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5, size=14, face="bold"))
  return(vplot)
}

###### GSEA dot plot #########
dot.plot <- function(df, n=10, title, savefile, subtitle=NULL, ...){
  if (nrow(df)==0){return("Empty df")}
  # by NES
  result <- df %>%
    filter(p.adjust < 0.05) %>%
    mutate(ordering = abs(NES)) %>%
    group_by(sign(NES)) %>%
    slice_max(order_by = ordering, n=n, with_ties = F) %>%
    mutate(Description = stringr::str_trunc(Description, 60))
  
  plot <- ggplot(result, 
                 aes(x = NES, y = reorder(Description,NES))) + 
    geom_point(aes(size = setSize, color = p.adjust)) +
    theme_bw(base_size = 12) +
    scale_colour_gradient(limits=c(0, 0.05), low="#ff6a00", high= "#ffffff") +
    scale_y_discrete(labels = label_wrap_gen(30)) +
    ylab(NULL) +
    ggtitle(title, subtitle)+
    xlim(-5, 5) +labs(x = "Normalized Enrichment Score (NES)", face="bold")+
    theme(axis.text.y=element_text(size=12, colour = "black"), 
          axis.title=element_text(size=12),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
    geom_vline(xintercept=0, linetype="dotted")+
    geom_segment(aes(x=0, xend = NES, y = Description, yend= Description), linetype=2) 
  
  if (!is.null(savefile)){
    ggsave(savefile, ...)
  }
  
  return (plot)  
}
```
