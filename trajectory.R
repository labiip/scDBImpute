
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("monocle")
rm(list = ls())
gc()
library(monocle)
readsCount = read.csv('Petropoulos_expressoin.csv', header = T, row.names = 1)            
data <- as(as.matrix(readsCount), "dgCMatrix")
time = read.csv('petropoulos_time.csv',header = F)
time = as.data.frame(time[,2])
colnames(time)[1] = 'time'

time = time[,1]
time = as.character(time)
time = as.factor(time)

sce <- CreateSeuratObject(counts = data, project = "klein", min.cells = 3, min.features = 200)
sce

library(monocle)
#导入注释好的seurat对象（已注释）  
##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息) 
expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix') 
##提取表型信息到p_data(phenotype_data)里面
pbmc@meta.data$orig.ident <- time
p_data <- pbmc@meta.data
##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加 
##提取基因信息 如生物类型、gc含量等 
f_data <- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc)) 
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)  
#构建CDS对象 
pd <- new('AnnotatedDataFrame', data = p_data)  
fd <- new('AnnotatedDataFrame', data = f_data) 
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。 
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())



cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2,reduction_method ='ICA')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds,color_by = "State")
plot_cell_trajectory(cds,color_by = "cell_type2")

