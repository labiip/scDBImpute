rm(list = ls())
gc()

X = read.csv("klein_expression.csv",row.names = 1)
label = read.csv('Klein_cell_label.csv',header = T,row.names = 1)

class.label<-as.matrix(label)
lpsdata=as.matrix(X)
yan1<-lpsdata
label<-class.label
#label = label+1
## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE---------------
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png')

## ---- message=FALSE, warning=FALSE-----------------------------------------
library(SingleCellExperiment)
library(SC3)
library(scater)
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(yan1),
    logcounts = log2(as.matrix(yan1) + 1)
    #logcounts = yan1
  ), 
  colData = label
)

# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

# define spike-ins
#isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)
is.spike <- grepl("^ERCC", rownames(sce))
sce <- splitAltExps(sce, ifelse(is.spike, "ERCC", "gene"))

## --------------------------------------------------------------------------
#dir.create(tempdir())
sce <- sc3(sce,gene_filter = FALSE, ks = 4, biology = TRUE)

sc3_n_clusters <- paste0("sc3_", as.character(4),"_clusters")
sc3_n_log2_outlier_score <- paste0("sc3_", as.character(4), "_log2_outlier_score")

sc3_n_clusters
sc3_n_log2_outlier_score

sce <- runPCA(sce)
plotReducedDim(sce, 
               "PCA", 
               colour_by = sc3_n_clusters, 
               size_by = sc3_n_log2_outlier_score)

plotReducedDim(sce, 
               "PCA", 
               colour_by = celltype)

sce <- runTSNE(sce, perplexity=10)
plotTSNE(sce, colour_by = celltype)

sce <- runUMAP(sce)
plotUMAP(sce, colour_by = celltype)

sc3_plot_consensus(
  sce, k = 4, 
  show_pdata = c(
    "cell_type1", 
    "log10_total_features",
    sc3_n_clusters, 
    sc3_n_log2_outlier_score
  )
)
#聚类对角性质的度量，对角线全部为红色 则为最好
#sc3_plot_silhouette(sce, k = size_row)
#绘制相关聚类的基因表达量图谱
sc3_plot_expression(
  sce, k = 4, 
  show_pdata = c(
    "cell_type1", 
    "log10_total_features",
    sc3_n_clusters, 
    sc3_n_log2_outlier_score
  )
)
#计算簇内稳定性
sc3_plot_cluster_stability(sce, k = 4)
#并绘制了最低p值的50个基因的基因表达谱
sc3_plot_de_genes(sce, k = 4)

col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
rezult_labels <- colData(sce)[[sc3_n_clusters, exact = FALSE]]
rezult_labels
#label = label[,2]
library("mclust")
if (require("mclust")) {
  adjustedRandIndex(label, rezult_labels)
}
#sc3_export_results_xls(sce)
library("aricode")
x<-as.vector(as.matrix(label))
y<-as.vector(as.matrix(rezult_labels))
NMI(x, y)
