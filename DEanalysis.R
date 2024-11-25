library("edgeR")
library("statmod")

chu = read.csv('chu_expression.csv',header = T,row.names = 1)
celltype_label = read.csv('D:/MATLAB/bin/chu/chu_celltype_label.csv',header = F)
celltype = as.data.frame(celltype_label[,2])
colnames(celltype)[1] = 'cell type'

c1 = which(celltype == 'H1 ESC')
c2 = which(celltype == 'DEC')
H1 = chu_sc[,c1]
DEC = chu_sc[,c2]
countData = cbind(H1, DEC)

countData = round(countData)
group <- c(rep("H1",212),rep("DEC",138))

#数据预处理
#（1）构建 DGEList 对象
dgelist <- DGEList(counts = countData, group = group)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
table(keep)

dgelist$dgelist$lib.size <- colSums(dgelist$counts)

#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

# 将归一化后的数据赋值给 dge 变量
dge = dgelist_norm

# 创建设计矩阵，用于指定差异分析模型
design <- model.matrix(~0 + factor(group))
rownames(design) <- colnames(dge)
colnames(design) <- levels(factor(group))

#（1）估算基因表达值的离散度
# 估计数据的离散度 —— common离散度、trended离散度、tagwise离散度
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

#dge <- estimateDisp(dge, design, robust = TRUE)

# 在估计的模型基础上进行 广义线性模型 (GLM) 拟合
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, contrast = c(-1, 1))
nrDEG <- topTags(lrt, n = nrow(dge))

# 将差异表达基因结果转换为数据框形式
DEG_edgeR <- as.data.frame(nrDEG)
res1 = DEG_edgeR[1:200,]


library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db) ##加载小鼠
library(org.Hs.eg.db) ##加载人类

genes = rownames(res1)
out = chu_sc[genes,]
out = chu_bulk[genes,]
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(out,entrezID=entrezIDs)
rt=out[is.na(out[,"entrezID"])==F,]                                 #去除基因id为NA的基因
gene=rt$entrezID
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)



go.res <- data.frame(kk) # 将GO结果转为数据框，方便后续分析（不转也可以，看个人习惯）
goBP <- subset(go.res,subset = (ONTOLOGY == "BP"))[1:10,]
goCC <- subset(go.res,subset = (ONTOLOGY == "CC"))[1:10,]
goMF <- subset(go.res,subset = (ONTOLOGY == "MF"))[1:10,]

go.df <- rbind(goBP,goCC,goMF)
go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))

library(stringr)
goBP$logP = -log10(goBP$p.adjust)
goBP=goBP[order(goBP$logP,decreasing = T),]
goBP$Description <- factor(goBP$Description,levels = rev(goBP$Description))
go_bar <- ggplot(data = goBP, # 绘图使用的数据
                 aes(x = Description, y = logP,fill = p.adjust))+ # 横轴坐标及颜色分类填充
  theme_classic()+
  scale_fill_viridis_c()+
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = NULL,y = "-log10(P.adj)",title = "Biological process")+ # 设置坐标轴标题及标题
  theme(text=element_text(size=16,  family="serif"))+
  theme(axis.title = element_text(size = 20), # 坐标轴标题大小
        axis.text.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold",angle = 90,vjust = 0.75,hjust = 0.85), # 坐标轴标签大小
        plot.title = element_text(size = 24,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 18), # 图例标题大小
        legend.text = element_text(size = 18), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
go_bar


goCC$logP = -log10(goCC$p.adjust)
goCC=goCC[order(goCC$logP,decreasing = T),]
goCC$Description <- factor(goCC$Description,levels = rev(goCC$Description))
go_bar <- ggplot(data = goCC, # 绘图使用的数据
                 aes(x = Description, y = logP,fill = p.adjust))+ # 横轴坐标及颜色分类填充
  theme_classic()+
  scale_fill_viridis_c()+
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  #coord_flip()+
  theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 30))+ # 设置term名称过长时换行
  labs(x = NULL,y = "-log10(P.adj)",title = "Cellular component")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 20), # 坐标轴标题大小
        axis.text.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold",angle = 90,vjust = 0.75,hjust = 0.85), # 坐标轴标签大小
        plot.title = element_text(size = 24,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 18), # 图例标题大小
        legend.text = element_text(size = 18), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
go_bar


goMF$logP = -log10(goMF$p.adjust)
goMF=goMF[order(goMF$logP,decreasing = T),]
goMF$Description <- factor(goMF$Description,levels = rev(goMF$Description))
go_bar <- ggplot(data = goMF, # 绘图使用的数据
                 aes(x = Description, y = logP,fill = p.adjust))+ # 横轴坐标及颜色分类填充
  theme_classic()+
  scale_fill_viridis_c()+
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  #coord_flip()+
  theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 30))+ # 设置term名称过长时换行
  labs(x = NULL,y = "-log10(P.adj)",title = "Molecular function")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 20), # 坐标轴标题大小
        axis.text.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold",angle = 90,vjust = 0.75,hjust = 0.85), # 坐标轴标签大小
        plot.title = element_text(size = 24,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 18), # 图例标题大小
        legend.text = element_text(size = 18), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
go_bar


# 差异基因热图
library(pheatmap)
countData = as.data.frame(countData)
countData = log1p(countData)
deg_opt <- DEG_edgeR %>% filter(DEG_edgeR$change != "stable")
deg_opt <- DEG_edgeR[1:200,] %>% filter(DEG_edgeR[1:200,]$change != "stable")
exp_brca_heatmap <- countData %>% filter(rownames(countData) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(exp_brca_heatmap) 

p1 <- pheatmap(exp_brca_heatmap, show_colnames = F, show_rownames = F,
               #scale = "row",
               cluster_cols = F,
               border = F,
               annotation_col = annotation_col,breaks = seq(0, 10, length.out = 100)) 
p1



