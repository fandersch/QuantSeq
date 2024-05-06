if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("DESeq2")
# BiocManager::install("edgeR")
# BiocManager::install("apeglm")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("annotate")
# BiocManager::install("EnhancedVolcano")
# 
# install.packages("PerformanceAnalytics")
# install.packages("naniar")
# install.packages("ggrepel")
# install.packages("RColorBrewer")
# install.packages("annotate")
# install.packages("pheatmap")
# install.packages("cluster")
# options(connectionObserver = NULL)

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(PerformanceAnalytics)
library(tidyverse)
library(DESeq2)
library(apeglm)
library(matrixStats)
library(ggrepel)
library(annotate)
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(cluster)
library(EnhancedVolcano)

args = commandArgs(trailingOnly=TRUE)


# output_dir <- "DEA_mm"
# input_file <- "input/input_mm.txt"
# count_file <- "featureCounts_umi_mm/counts.tsv"
# species <- "mm"
# filter_out_samples <- ""

output_dir <- args[1]
input_file <- args[2]
count_file <- args[3]
species <- args[4]

sample_ids_filtered<-""
if(!is.na(args[5])){
  sample_ids_filtered <- str_split(args[5], pattern = ",") %>% unlist %>% as.numeric
}

apply_lfc_shrink <- F
if(!is.na(args[6])){
  if( args[6] == "lfcShrinkage"){
    apply_lfc_shrink <- T
  }
}

dir.create(output_dir)

if(species == "hs"){
  #human input
  cts_human <- as.matrix(read.csv(count_file, sep="\t",row.names="entrez_id"))
  cts_human<-cts_human[,7:ncol(cts_human)]
  class(cts_human) <- "numeric"
  
  gene_symbols_human <- data.frame(getSYMBOL(rownames(cts_human[,]), data='org.Hs.eg.db'), rownames(cts_human[,]))
  colnames(gene_symbols_human)<- c("gene_symbol", "entrez_id")
}

if(species == "mm"){
  #mouse input
  cts_mouse <- as.matrix(read.csv(count_file, sep="\t",row.names="entrez_id"))
  cts_mouse<-cts_mouse[,7:ncol(cts_mouse)]
  class(cts_mouse) <- "numeric"
  
  gene_symbols_mouse <- data.frame(getSYMBOL(rownames(cts_mouse[,]), data='org.Mm.eg.db'), rownames(cts_mouse[,]))
  colnames(gene_symbols_mouse)<- c("gene_symbol", "entrez_id")
}


#get sampleIDs
sample_info <- read_tsv(input_file) %>%
  rowwise %>%
  distinct

sample_info$`condition` <- gsub('(_|/| )', '.', sample_info$`condition`)
sample_info$`group` <- gsub('(_|/| )', '.', sample_info$`group`)
#get sample groups
sample_groups <- sample_info$group %>% unique

#filter out bad samples
sample_info <- sample_info %>%
  filter(!sample_id %in% sample_ids_filtered)

###################
# PCA all samples #
###################
output_dir_group<- paste0(output_dir, "/PCA_all/")
#specify output directory
dir.create(output_dir_group, recursive = T)

#make ordered metadata matrix
coldata <- sample_info %>%
  dplyr::select(sample_id, condition) %>%
  distinct %>%
  arrange(sample_id) %>%
  column_to_rownames(var = "sample_id")

coldata$condition <- factor(coldata$condition)
cols <- rownames(coldata)

if(species == "mm"){
  cts <- cts_mouse
  gene_symbols <- gene_symbols_mouse
}
if(species == "hs"){
  cts <- cts_human
  gene_symbols <- gene_symbols_human
}

colnames(cts) <- colnames(cts) %>% str_remove(pattern="^X")

#make dds from cts, coldata
dds <- DESeqDataSetFromMatrix(countData = cts[,as.character(c(cols))],
                              colData = coldata,
                              design = ~ condition)

#perform diff expression analysis
dds <- DESeq(dds)

# Input is a matrix of log transformed values
vst <- varianceStabilizingTransformation(dds, blind=T)
vst_mat <- assay(vst)
vst_mat <- vst_mat %>% cbind(rowVars(vst_mat))
vst_mat_ordered <- vst_mat[order(vst_mat[,ncol(vst_mat)], decreasing = T),]

pca <- prcomp(t(vst_mat_ordered[1:500,1:ncol(vst_mat)-1]))
sdev<-pca$sdev
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(coldata, pca$x) %>% 
  cbind(sdev) %>% 
  rownames_to_column(var="sample_id")

write_tsv(df, file=paste0(output_dir_group, "matrix_PCA.txt"), col_names = T)

PoV = df$sdev ^ 2 / sum(pca$sdev ^ 2) * 100

#plot pca with labels
ggplot(df, aes(x= PC1, y= PC2, colour=condition, label=paste0(sample_id,"_",condition))) +
  xlab(paste("PC1 (", round(PoV[1],digits=2), " % variance)",sep="")) +
  ylab(paste("PC2 (", round(PoV[2],digits=2), " % variance)",sep="")) +
  geom_point() + 
  geom_text_repel(aes(label=paste0(sample_id,"_",condition)), size=2, max.overlaps=100, force = 2) + 
  theme(aspect.ratio=3/4)

ggsave(paste0(output_dir_group, "PCA_sample_id_PC1_PC2.png"), height=10,width=10)


#############
# group DEA #
#############

df_LFC_all_group <- NULL
df_LFC_all <- NULL
# dev.off()
for(z in 1:length(sample_groups)){
  print(sample_groups[z])
  output_dir_group<- paste0(output_dir, "/", sample_groups[z], "/")
  #specify output directory
  dir.create(output_dir_group, recursive = T)
  
  #get group species 
  species <- sample_info %>%
    filter(group==sample_groups[z]) %>%
    .$species %>%
    unique
  
  #get sample conditions per group
  sample_conditions <- sample_info %>%
    filter(group==sample_groups[z]) %>%
    dplyr::select(condition, control) %>%
    distinct
  
  #get tretment conditions
  sample_group_treatments <- sample_conditions %>%
    filter(control == 0) %>%
    .$condition
  
  #get control
  sample_group_control <- sample_conditions %>%
    filter(control == 1) %>%
    .$condition
  
  #make ordered metadata matrix
  coldata <- sample_info %>%
    filter(group==sample_groups[z]) %>% 
    arrange(sample_id) %>%
    column_to_rownames(var = "sample_id") %>%
    dplyr::select(condition)
  
  coldata$condition <- factor(coldata$condition)
  cols <- rownames(coldata)
  
  if(species == "mouse"){
    cts <- cts_mouse
    gene_symbols <- gene_symbols_mouse
  }
  if(species == "human"){
    cts <- cts_human
    gene_symbols <- gene_symbols_human
  }
  
  colnames(cts) <- colnames(cts) %>% str_remove(pattern="^X")
  
  #make dds from cts, coldata
  dds <- DESeqDataSetFromMatrix(countData = cts[,as.character(c(cols))],
                                colData = coldata,
                                design = ~ condition)
  
  #set ref condition
  if(length(sample_group_control)>0){
    dds$condition <- relevel(dds$condition, ref = sample_group_control)
  }
 
  
  #perform diff expression analysis
  dds <- DESeq(dds)
  
  #  Export PCA plot (default deseq PCA on 500 most variable genes).
  # png(paste0(output_dir_group,sample_groups[z], "_PCA.png"),width = 300, height = 300, units='mm', res = 100)
  # plotPCA(varianceStabilizingTransformation(dds, blind=T), intgroup = "condition", ntop = 4000)
  # plotPCA(varianceStabilizingTransformation(dds, blind=T), intgroup = "condition", ntop = 300)
  # ggsave(paste0(output_dir_group,sample_groups[z], "_PCA.png"))
  # dev.off()
  
  # Input is a matrix of log transformed values
  vst <- varianceStabilizingTransformation(dds, blind=T)
  vst_mat <- assay(vst)
  vst_mat <- vst_mat %>% cbind(rowVars(vst_mat))
  vst_mat_ordered <- vst_mat[order(vst_mat[,ncol(vst_mat)], decreasing = T),]
  
  pca <- prcomp(t(vst_mat_ordered[1:500,1:ncol(vst_mat)-1]))
  sdev<-pca$sdev
  # Create data frame with metadata and PC3 and PC4 values for input to ggplot
  df <- cbind(coldata, pca$x) %>% 
    cbind(sdev) %>% 
    rownames_to_column(var="sample_id")
  
  write_tsv(df, file=paste0(output_dir_group,sample_groups[z], "matrix_PCA.txt"), col_names = T)
  
  PoV = df$sdev ^ 2 / sum(pca$sdev ^ 2) * 100
  
  
  #plot pca with labels
  ggplot(df, aes(x= PC1, y= PC2, colour=condition, label=sample_id)) +
    xlab(paste("PC1 (", round(PoV[1],digits=2), " % variance)",sep="")) +
    ylab(paste("PC2 (", round(PoV[2],digits=2), " % variance)",sep="")) +
    geom_point() + 
    geom_label_repel(aes(label=sample_id), size=5) + 
    theme(aspect.ratio=3/4)

  ggsave(paste0(output_dir_group,sample_groups[z], "_PCA_sample_id_PC1_PC2.png"), height=10,width=10)

  resultsNames(dds)
  norm_counts <- counts(dds, normalized=T) %>%
    round(2) %>%
    as.data.frame() %>%
    rownames_to_column("entrez_id") %>%
    left_join(gene_symbols) %>%
    dplyr::select(entrez_id, gene_symbol, everything()) %>%
    arrange(gene_symbol)
    
  write_tsv(norm_counts, file=paste0(output_dir_group,sample_groups[z], "_normalized_counts.tsv"))
  
  #################################################
  #hierachical clustering 10% most variable genes #
  #################################################
  #sample distance heatmap
  vst <- varianceStabilizingTransformation(dds, blind=T)
  vst_mat <- assay(vst)
  sampleDists <- dist(t(vst_mat), method = "euclidean")
  sampleDistMatrix <- as.matrix(sampleDists)
  
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  pdf(file = paste0(output_dir_group,sample_groups[z],"all_genes_sample_distance.pdf"), width = 10, height=10)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  dev.off()
  
  # Calculate the variance for each gene
  variances <- apply(vst_mat, 1, var)
  mean <- apply(vst_mat, 1, mean)
  
  # Determine the upper quartile variance/mean cutoff value
  upper_var <- quantile(variances, 0.25)
  upper_mean <- quantile(mean, 0.25)
  
  # Filter the data choosing only genes that are variable
  df_by_var_mean <- data.frame(vst_mat) %>%
    dplyr::filter(variances > upper_var, mean > upper_mean)
  
  colors <- colorRampPalette(brewer.pal(9, "RdBu"))(255)
  
  pdf(file = paste0(output_dir_group,sample_groups[z],"all_genes_clustering_heatmap.pdf"),height=10,width=10)
  pheatmap(
    df_by_var_mean,
    cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
    cluster_cols = TRUE, # Cluster the columns of the heatmap (samples),
    show_rownames = FALSE, # There are too many genes to clearly show the labels
    col=colors,
    scale = "row" # Scale values in the direction of genes (rows)
  )
  dev.off()
  
  if(length(sample_group_control)>0){
    for(i in 2:length(resultsNames(dds))){
      print(resultsNames(dds)[i])
      sample_group_treatment <- resultsNames(dds)[i] %>% str_replace("condition_", "") %>% str_replace(paste0("_vs_", sample_group_control), "")
      #results
      # res <- results(dds, name=resultsNames(dds)[i], lfcThreshold = 0.585, alpha=0.05) 
      res <- results(dds, name=resultsNames(dds)[i]) 
      
      #LFC shrinkage
      if(apply_lfc_shrink){
        res <- lfcShrink(dds, res=res, coef=resultsNames(dds)[i], type="apeglm")
      }
      
      #MA-plot
      png(paste0(output_dir_group,resultsNames(dds)[i], "_MA.png"),width = 400, height = 300, units='mm', res = 100)
      DESeq2::plotMA(res, ylim=c(-7,7))
      dev.off()
    
      #Plot counts
      #plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
      
      #volcano plot
      # create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
      # this can be achieved with nested ifelse statements
      keyvals <- ifelse(
        res$log2FoldChange < -1 & -log10(res$pvalue) > 5, 'blue',
        ifelse(res$log2FoldChange > 1 & -log10(res$pvalue) > 5, 'red',
               'black'))
      keyvals[is.na(keyvals)] <- 'black'
      names(keyvals)[keyvals == 'red'] <- 'up'
      names(keyvals)[keyvals == 'black'] <- 'neutral'
      names(keyvals)[keyvals == 'blue'] <- 'down'
      
      labels <- data.frame(entrez_id=rownames(res)) %>%
        left_join(gene_symbols)
      
      EnhancedVolcano(res,
                      lab = labels$gene_symbol,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      colCustom=keyvals,
                      labSize = 2.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.2,
                      max.overlaps=25,
                      arrowheads=F
      )
      
      ggsave(paste0(output_dir_group,resultsNames(dds)[i], "_volcano.pdf"),height=10,width=10)

      col_treatment <- sample_info %>%
        filter(group==sample_groups[z]) %>%
        filter(condition == sample_group_treatments[i-1]) %>%
        .$sample_id %>% 
        as.character
      
      col_control <- sample_info %>%
        filter(group==sample_groups[z]) %>%
        filter(condition == sample_group_control) %>%
        .$sample_id %>% 
        as.character
      
      df <- norm_counts[,c("entrez_id", "gene_symbol", col_treatment, col_control)] %>% as.data.frame
      df$treatment_mean <- rowMeans(norm_counts[,c(col_treatment)]) %>% unlist
      df$control_mean <- rowMeans(norm_counts[,c(col_control)]) %>% unlist
  
      df <- df %>%
        full_join(
          as.data.frame(res) %>% 
          rownames_to_column(var = "entrez_id"),
          by= "entrez_id") %>% 
        mutate_at(vars(c("baseMean", "log2FoldChange", "lfcSE", "treatment_mean", "control_mean", all_of(col_treatment), all_of(col_control))), round, 2) %>%
        dplyr::select(gene_symbol, entrez_id, everything()) %>%
        arrange(padj)
        
      #Exporting results to CSV files
      write_tsv(df, file=paste0(output_dir_group,resultsNames(dds)[i], ".tsv"))
      
      #correlation counts replicates
      png(paste0(output_dir_group,resultsNames(dds)[i], "_corr_replicates.png"),width = 250, height = 250, units='mm', res = 200)
      
      norm_counts_analytics<- norm_counts[,c(col_treatment, col_control)] %>%
        as.data.frame() %>%
        mutate(rowMean = rowMeans(norm_counts[,c(col_treatment, col_control)])) %>%
        filter(rowMean>0) %>%
        dplyr::select(-rowMean)
      
      chart.Correlation(norm_counts_analytics, histogram=TRUE, pch=19)
      dev.off()
      
      #pepare LFC matrix 
      df_LFC <- df %>%
        dplyr::select(gene_symbol, entrez_id, log2FoldChange, pvalue, padj)
      
      names(df_LFC)[names(df_LFC) == "log2FoldChange"] <- paste0(resultsNames(dds)[i], "_LFC")
      names(df_LFC)[names(df_LFC) == "pvalue"] <- paste0(resultsNames(dds)[i], "_pvalue")
      names(df_LFC)[names(df_LFC) == "padj"] <- paste0(resultsNames(dds)[i], "_padj")
      
      if(is.null(df_LFC_all_group)){
        df_LFC_all_group <- df_LFC
      }else{
        df_LFC_all_group <- df_LFC_all_group %>%
          full_join(df_LFC) 
      }
      
      if(is.null(df_LFC_all)){
        df_LFC_all <- df_LFC
      }else{
        df_LFC_all <- df_LFC_all %>%
          full_join(df_LFC) 
      }
    }
    write_tsv(df_LFC_all_group %>% arrange(gene_symbol), file=paste0(output_dir_group,"LFCs_Pvalue_Padj_", sample_groups[z],".tsv"), na = "")
    df_LFC_all_group<-NULL
  }
}

write_tsv(df_LFC_all %>% arrange(gene_symbol), file=paste0(output_dir,"/LFCs_Pvalue_Padj_all.tsv"), na = "")
