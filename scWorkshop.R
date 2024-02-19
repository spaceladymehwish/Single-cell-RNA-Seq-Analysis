#install a group of packages collectively 
install.packages(c("ggplot2", "BiocManager", "sctransform",
          "devtools", "cowplot", "matrixStats",
          "ggbeeswarm", "ggnewscale", "msigdbr", "bit64", "scater",
          "AnnotationDbi",
          "SingleR", "clusterProfiler", "celldex",
          "dittoSeq", "DelayedArray",
          "DelayedMatrixStats",
          "limma", "SingleCellExperiment",
          "SummarizedExperiment",
          "slingshot", "batchelor",
          "clustree", "edgeR"))

#or install a group of packages collectively

# Load the packages
pkgs <- c("ggplot2", "BiocManager", "sctransform",
          "devtools", "cowplot", "matrixStats",
          "ggbeeswarm", "ggnewscale", "msigdbr", "bit64", "scater",
          "AnnotationDbi",
          "SingleR", "clusterProfiler", "celldex",
          "dittoSeq", "DelayedArray",
          "DelayedMatrixStats",
          "limma", "SingleCellExperiment",
          "SummarizedExperiment",
          "slingshot", "batchelor",
          "clustree", "edgeR")
### if packages need to be installed from github:
### devtools::install_github("namespace/repo")
## install Bioconductor --------------------------------------------------------
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

## install and check package loading -------------------------------------------
for (pkg in basename(pkgs)) {
  BiocManager::install(pkg, ask = FALSE, update = FALSE)
  
  if (! library(pkg, character.only = TRUE, logical.return = TRUE)) {
    write(paste0("Installation of package ",
                 pkg,
                 " exited with non-zero exit status"),
          stdout())
    quit(status = 1, save = "no")
  }
}


pkgs

for (pkg in pkgs) {
  library(pkg, character.only = TRUE)
}

##Set Working Directory and Read Data Set

setwd("D:/Projects/Data Analysis on ScRNA-Seq/scRNA-Seq-Analysis/Sample dataset")
sampleinfo <- read.csv("sampleinfo.csv")
datadirs <- file.path(".", sampleinfo$SampleName)
names(datadirs) <- gsub("_", "-", sampleinfo$SampleName)
datadirs

sparse_matrix <- Seurat::Read10X(data.dir = datadirs)

seu <- Seurat::CreateSeuratObject(counts = sparse_matrix,
                                  project = "Primary",
                                  min.cells = 3,
                                  min.features = 100)

metadata_table <- seu@meta.data
write.csv(metadata_table, "Meta_Data_Table.csv", row.names=TRUE)
Cells <- table(seu@active.ident)
write.csv(Cells, "Cells.csv", row.names=TRUE)

##Preprocess and Quality Control
  

mypalette <- c("#FF0000", "#00796B", "#27CED7", "#FF00CC", "#00FF00","purple","pink","orange","yellow")
mypalette2 <- c("#27CED7", "#FF00CC", "#00FF00","#FF2A00", "#00468B","purple","pink","orange","yellow")
mypalette3 <- c("#27CED7", "#FF00CC", "#00FF00","#FF2A00", "#00468B","purple","pink","orange","yellow","#FF0000", "#00796B","brown")


hist(seu$nCount_RNA, col = mypalette, main = paste0("Histogram of RNA counts per cell"))
hist(seu$nFeature_RNA, col = mypalette, main = paste0("Histogram of Gene counts per cell"))


##Feature Analysis
  
 
Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", col = mypalette2)
Seurat::VlnPlot(seu, features = "nCount_RNA", cols = mypalette2)
Seurat::VlnPlot(seu, features = "nFeature_RNA", cols = mypalette)

Seurat::VlnPlot(seu, features = c("nCount_RNA",
                                  "nFeature_RNA"))

seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^MT-",
                                    col.name = "percent.mito")

seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^RP[SL]",
                                    col.name = "percent.ribo")

seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^HB[^(P)]",
                                    col.name = "percent.globin")
head(seu@meta.data)



Seurat::VlnPlot(seu, features = "percent.mito", cols = mypalette2, pt.size = 0)

Seurat::VlnPlot(seu, features = "percent.ribo", cols = mypalette2, pt.size = 0)


Seurat::VlnPlot(seu, features = "percent.globin", cols = mypalette2, pt.size = 0)

Seurat::VlnPlot(seu, features = c("percent.mito",
                                  "percent.ribo",
                                  "percent.globin"))


Seurat::FeatureScatter(seu,
                       feature1 = "percent.mito",
                       feature2 = "percent.ribo", cols = mypalette2)




##########################Feature Selection and Visualization
  

library(ggplot2)
library(Matrix)
library(Seurat)

most_expressed_boxplot <- function(object, ngenes = 20){
  
  cts <- Seurat::GetAssayData(object, assay = "RNA", slot = "counts")
  
  cts <- t(cts)/colSums(cts)*100
  medians <- apply(cts, 2, median)
  
  most_expressed <- order(medians, decreasing = T)[ngenes:1]
  most_exp_matrix <- as.matrix((cts[,most_expressed]))
  
  most_exp_df <- stack(as.data.frame(most_exp_matrix))
  colnames(most_exp_df) <- c("perc_total", "gene")
  
  boxplot <- ggplot(most_exp_df, aes(x=gene, y=perc_total)) +
    geom_boxplot() +
    coord_flip()
  return(boxplot)
}
most_expressed_boxplot(seu, 20)


#########################Data Filtering###########
#nfeature_RNA: Number of genes detected per cell greater than 500 and less than 2500 and percent.mito: Percentage of mitochondrial genes less than 15%.  

seu <- subset(seu, subset = nFeature_RNA > 500 &
                nFeature_RNA < 2500 &
                percent.mito < 15)

Seurat::VlnPlot(seu, features = c("nFeature_RNA",
                                  "percent.mito"))


########################Data Normalization#########
  
Seurat::GetAssayData(seu)[1:30,1:30] 
Assay_Data_Before_Normalization <- Seurat::GetAssayData(seu)[1:30,1:30]
write.csv(Assay_Data_Before_Normalization, "Assay_Data_Before_Normalization.csv", row.names=TRUE)


#########################Identify Variable Features#######
  

seu <- Seurat::NormalizeData(seu,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)
seu <- Seurat::FindVariableFeatures(seu,
                                    selection.method = "vst",
                                    nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(seu), 10)
top10

vf_plot <- Seurat::VariableFeaturePlot(seu)
Seurat::LabelPoints(plot = vf_plot,
                    points = top10, repel = TRUE)


############################Data Scaling########
  

seu <- Seurat::ScaleData(seu)
seu <- Seurat::ScaleData(seu,
                         features = rownames(seu))

plots <- Seurat::VlnPlot(seu, features = c("percent.mito",
                                   "percent.ribo",
                                   "percent.globin"),
                 pt.size =0,
                 combine = FALSE, cols = mypalette2)

for(i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + geom_boxplot() + theme(legend.position = 'none')
}

CombinePlots(plots)

###############Principal Component Analysis (PCA)#######
  
seu <- Seurat::RunPCA(seu)
Seurat::DimPlot(seu, reduction = "pca")
Seurat::DimPlot(seu, reduction = "pca", cols = mypalette2, pt.size = 1)
Seurat::DimHeatmap(seu, dims = 1:12, cells = 500, balanced = TRUE)
Seurat::ElbowPlot(seu, ndims = 40)


################## Dimensionality Reduction -- UMAP#####
#Run UMAP to visualize the data in lower-dimensional space.


seu <- Seurat::RunUMAP(seu, dims = 1:25)

Seurat::DimPlot(seu, reduction = "umap", cols = mypalette2)
seu <- Seurat::RunUMAP(seu, dims = 1:25, n.neighbors = 30) 
Seurat::DimPlot(seu, reduction = "umap", cols = mypalette2)
seu <- Seurat::RunUMAP(seu, dims = 1:5)      
Seurat::DimPlot(seu, reduction = "umap", cols = mypalette2)
seu <- Seurat::RunUMAP(seu, dims = 1:50)     
Seurat::DimPlot(seu, reduction = "umap", cols = mypalette2)

seu <- Seurat::RunUMAP(seu, dims = 1:30)   
Seurat::DimPlot(seu, reduction = "umap", cols = mypalette2)
seu <- Seurat::RunTSNE(seu, dims = 1:30)  
Seurat::DimPlot(seu, reduction = "tsne", group.by = 'orig.ident', cols = mypalette2)

###############Integration#################  
 


library(harmony) 
seu <- seu %>% RunHarmony('orig.ident', plot_convergence = F)


library(ggplot2)
library(tidyverse)
seu.embed <- Embeddings(seu, "harmony")


############Clustering Analysis############

library(ggraph)
library(clustree)
seu_clusters_UMAP <- seu %>%
  RunUMAP(reduction = "harmony", dims = 1:25) %>%
  FindNeighbors(reduction = "harmony", dims = 1:25) %>%
  FindClusters(resolution = seq(0.1, 0.8, by=0.1) )

seu_clusters_TSNE <- seu %>%
  RunTSNE(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 0.8, by=0.1) ) 

clustree::clustree(seu_clusters_TSNE@meta.data[,grep("RNA_snn_res", colnames(seu_clusters_TSNE@meta.data))],
                   prefix = "RNA_snn_res.")


#########################Visualization of Clusters################
  
 

Clusters_tSNE <- DimPlot(seu_clusters_TSNE, reduction = 'tsne', group.by = 'RNA_snn_res.0.3', raster = FALSE, cols = mypalette3 ) 
Clusters_tSNE
library(celldex)
library(SingleR)
seu_clusters_TSNE <- Seurat::SetIdent(seu_clusters_TSNE, value = seu_clusters_TSNE$RNA_snn_res.0.3)

tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")
Seurat::VlnPlot(seu_clusters_TSNE,
                features = tcell_genes,
                ncol = 2, pt.size = 0, cols = mypalette3)
seu_clusters_TSNE <- Seurat::AddModuleScore(seu_clusters_TSNE,
                                            features = list(tcell_genes),
                                            name = "tcell_genes")
Seurat::VlnPlot(seu_clusters_TSNE,
                "tcell_genes1",
                pt.size = 0, cols = mypalette3)


########################Cell Cycle Scoring#################

mypalette4 <- c("#27CED7", "#FF0000", "#FF9900")

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
seu_clusters_TSNE <- Seurat::CellCycleScoring(seu_clusters_TSNE,
                                              s.features = s.genes,
                                              g2m.features = g2m.genes)
Seurat::DimPlot(seu_clusters_TSNE, reduction = 'tsne', group.by = "Phase",
                cols = mypalette4, label = TRUE)


ref <- celldex::HumanPrimaryCellAtlasData()
class(ref)
table(ref$label.main)


###########################SingleR Analysis############################

seu_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu_clusters_TSNE, slot = "data"),
                                ref = ref,
                                labels = ref$label.main)
head(seu_SingleR)

SingleR::plotScoreHeatmap(seu_SingleR)

SingleR::plotDeltaDistribution(seu_SingleR)

singleR_labels <- seu_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- "none"

seu_clusters_TSNE$SingleR_annot <- singleR_labels


##Visualization of Results**

mypalette0 <- c("#00796B", "#27CED7", "#FF00CC", "#00FF00", "#00FFFF", "#FF0000", "#00468B", "#FDAF91", "#5050FF", "#350E20", "#999999", "#7B4173", "#FF9900", "#358000", "#0000CC", "#99CCFF", "#FFCCCC", "#004C00", "#CCFFFF", "#CC99FF", "#9900CC", "#996600", "#666600", "#CCFF00", "#FFCC00", "#000000", "#FF420E", "#79CC3D", "#7E0021", "#FFF7F3", "#6699FF", "#CCCC99")
dittoSeq::dittoDimPlot(seu_clusters_TSNE, reduction = 'tsne', "SingleR_annot", size = 0.7)
dittoSeq::dittoBarPlot(seu_clusters_TSNE, var = "SingleR_annot", group.by = "orig.ident")
dittoSeq::dittoBarPlot(seu_clusters_TSNE, 
                       var = "SingleR_annot", 
                       group.by = "RNA_snn_res.0.3")
Seurat::VlnPlot(seu_clusters_TSNE,
                features = "percent.ribo",
                pt.size = 0, cols = mypalette3)


#########Differential gene expression###########
library(edgeR)
library(limma)
de_genes <- Seurat::FindAllMarkers(seu_clusters_TSNE, min.pct = 0.25,
                                   only.pos = TRUE)

de_genes <- subset(de_genes, de_genes$p_val_adj<0.05)
write.csv(de_genes, "de_genes_FindAllMarkers.csv", row.names = F, quote = F)
write.table(de_genes, "de_genes_FindAllMarkers.txt", row.names = FALSE, quote
            = FALSE)
library(dplyr)
top_specific_markers <- de_genes %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)
top_specific_markers <- de_genes %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)
dittoSeq::dittoDotPlot(seu_clusters_TSNE, vars =
                         unique(top_specific_markers$gene),
                       group.by = "SingleR_annot")


