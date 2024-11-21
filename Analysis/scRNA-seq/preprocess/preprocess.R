library(Seurat)
library(readr)
library(ggplot2)
library(stringr)

run_Seurat <- function(dat) {
	dat <- NormalizeData(object = dat, normalization.method = "LogNormalize", 
		scale.factor = 10000)
	dat <- FindVariableFeatures(object = dat, 
		selection.method = "vst", nfeatures = 2000)
	dat <- ScaleData(object = dat, features = rownames(dat))
	dat <- RunPCA(object = dat, features = VariableFeatures(object = dat))
	dat <- FindNeighbors(object = dat, dims = 1:50, force.recalc = TRUE)
	dat <- FindClusters(object = dat, resolution = 0.6)
	dat <- RunUMAP(object = dat, dims = 1:50)
	dat <- RunTSNE(object = dat, dims = 1:50)
	return(dat)
}

####################################
#AITL1 (initial diagnosis, ID)
####################################
dir.create("./result/preprocess/AITL1")
AITL1 <- read.table(gzfile("./expression_matrix/AITL-1.expression_matrix.txt.gz"), sep = "\t", row.names = 1, header = TRUE)
AITL1.obj <- CreateSeuratObject(AITL1, project="AITL1", min.cells = 3)
AITL1.obj$percent.mt <- PercentageFeatureSet(AITL1.obj, pattern = "^MT-")
ribogenes <- union(rownames(AITL1.obj)[str_detect(rownames(AITL1.obj), "^RPS")],
	rownames(AITL1.obj)[str_detect(rownames(AITL1.obj), "^RPL")])
AITL1.obj$percent.ribo <- PercentageFeatureSet(AITL1.obj, feature = ribogenes)

qcplot <- VlnPlot(AITL1.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)
ggsave("./result/preprocess/AITL1/AITL1_qc.pdf",qcplot,width=30,height=30,units='cm')
pdf("./result/preprocess/AITL1/AITL1_density.pdf", onefile=TRUE)
for (i in c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) {
   plot(density(AITL1.obj@meta.data[[i]]))
}
dev.off()

AITL1.obj1 <- subset(AITL1.obj, subset=nFeature_RNA<6000 & nFeature_RNA > 200 & percent.mt<20)  #31308 4681
saveRDS(AITL1.obj1, "./result/preprocess/AITL1/AITL1.RDS")

####################################
#AITL2(initial diagnosis, ID)
####################################
rm(list=ls())
dir.create("./result/preprocess/AITL2")
AITL2 <- read.table(gzfile("./expression_matrix/AITL-2.expression_matrix.txt.gz"), sep = "\t", row.names = 1, header = TRUE)
AITL2.obj <- CreateSeuratObject(AITL2, project="AITL2", min.cells = 3)
AITL2.obj$percent.mt <- PercentageFeatureSet(AITL2.obj, pattern = "^MT-")
ribogenes <- union(rownames(AITL2.obj)[str_detect(rownames(AITL2.obj), "^RPS")],
	rownames(AITL2.obj)[str_detect(rownames(AITL2.obj), "^RPL")])
AITL2.obj$percent.ribo <- PercentageFeatureSet(AITL2.obj, feature = ribogenes)
qcplot <- VlnPlot(AITL2.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)
ggsave("./result/preprocess/AITL2/AITL2_qc.pdf",qcplot,width=30,height=30,units='cm')
pdf("./result/preprocess/AITL2/AITL2_density.pdf", onefile=TRUE)
for (i in c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) {
   plot(density(AITL2.obj@meta.data[[i]]))
}
dev.off()
AITL2.obj1 <- subset(AITL2.obj, subset = nFeature_RNA<6000 & nFeature_RNA > 200 & percent.mt<20) 
saveRDS(AITL2.obj1, "./result/preprocess/AITL2/AITL2.RDS")


####################################
#LHL(initial diagnosis, ID)
####################################
rm(list=ls())
dir.create("./result/preprocess/LHL")
LHL <- read.table(gzfile("./expression_matrix/LHL_200818_matrix.tsv.gz"), sep = "\t", row.names = 1, header = TRUE)
LHL.obj <- CreateSeuratObject(LHL, project="LHL", min.cells = 3)
LHL.obj$percent.mt <- PercentageFeatureSet(LHL.obj, pattern = "^MT-")
ribogenes <- union(rownames(LHL.obj)[str_detect(rownames(LHL.obj), "^RPS")],
	rownames(LHL.obj)[str_detect(rownames(LHL.obj), "^RPL")])
LHL.obj$percent.ribo <- PercentageFeatureSet(LHL.obj, feature = ribogenes)
qcplot <- VlnPlot(LHL.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)
ggsave("./result/preprocess/LHL/LHL_qc.pdf",qcplot,width=30,height=30,units='cm')
pdf("./result/preprocess/LHL/LHL_density.pdf", onefile=TRUE)
for (i in c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) {
   plot(density(LHL.obj@meta.data[[i]]))
}
dev.off()
LHL.obj1 <- subset(LHL.obj, subset = nFeature_RNA<6000 & nFeature_RNA > 200 & percent.mt<20) 
saveRDS(LHL.obj1, "./result/preprocess/LHL/LHL.RDS")


####################################
#AITL3
####################################
rm(list=ls())
dir.create("./result/preprocess/AITL3")
AITL3 <- read.table(gzfile("./expression_matrix/AITL-3.expression_matrix.txt.gz"), sep = "\t", row.names = 1, header = TRUE)
AITL3.obj <- CreateSeuratObject(AITL3, project="AITL3", min.cells = 3)
AITL3.obj$percent.mt <- PercentageFeatureSet(AITL3.obj, pattern = "^MT-")
ribogenes <- union(rownames(AITL3.obj)[str_detect(rownames(AITL3.obj), "^RPS")],
	rownames(AITL3.obj)[str_detect(rownames(AITL3.obj), "^RPL")])
AITL3.obj$percent.ribo <- PercentageFeatureSet(AITL3.obj, feature = ribogenes)
qcplot <- VlnPlot(AITL3.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)
ggsave("./result/preprocess/AITL3/AITL3_qc.pdf",qcplot,width=30,height=30,units='cm')
pdf("./result/preprocess/AITL3/AITL3_density.pdf", onefile=TRUE)
for (i in c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) {
   plot(density(AITL3.obj@meta.data[[i]]))
}
dev.off()
AITL3.obj1 <- subset(AITL3.obj, subset = nFeature_RNA<6000 & nFeature_RNA > 200 & percent.mt<20) 
saveRDS(AITL3.obj1, "./result/preprocess/AITL3/AITL3.RDS")


####################################
#AITL4
####################################
rm(list=ls())
dir.create("./result/preprocess/AITL4")
AITL4 <- read.table(gzfile("./expression_matrix/AITL_4.expression-matrix.txt.gz"), sep = "\t", row.names = 1, header = TRUE)
AITL4.obj <- CreateSeuratObject(AITL4, project="AITL4", min.cells = 3)
AITL4.obj$percent.mt <- PercentageFeatureSet(AITL4.obj, pattern = "^MT-")
ribogenes <- union(rownames(AITL4.obj)[str_detect(rownames(AITL4.obj), "^RPS")],
	rownames(AITL4.obj)[str_detect(rownames(AITL4.obj), "^RPL")])
AITL4.obj$percent.ribo <- PercentageFeatureSet(AITL4.obj, feature = ribogenes)
qcplot <- VlnPlot(AITL4.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)
ggsave("./result/preprocess/AITL4/AITL4_qc.pdf",qcplot,width=30,height=30,units='cm')
pdf("./result/preprocess/AITL4/AITL4_density.pdf", onefile=TRUE)
for (i in c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) {
   plot(density(AITL4.obj@meta.data[[i]]))
}
dev.off()
AITL4.obj1 <- subset(AITL4.obj, subset = nFeature_RNA<6000 & nFeature_RNA > 200 & percent.mt<20) 
saveRDS(AITL4.obj1, "./result/preprocess/AITL4/AITL4.RDS")


####################################
#AITL7
####################################
rm(list=ls())
dir.create("./result/preprocess/AITL7")
AITL7 <- read.table(gzfile("./expression_matrix/AITL-7.expression_matrix.txt.gz"), sep = "\t", row.names = 1, header = TRUE)
AITL7.obj <- CreateSeuratObject(AITL7, project="AITL7", min.cells = 3)
AITL7.obj$percent.mt <- PercentageFeatureSet(AITL7.obj, pattern = "^MT-")
ribogenes <- union(rownames(AITL7.obj)[str_detect(rownames(AITL7.obj), "^RPS")],
	rownames(AITL7.obj)[str_detect(rownames(AITL7.obj), "^RPL")])
AITL7.obj$percent.ribo <- PercentageFeatureSet(AITL7.obj, feature = ribogenes)
qcplot <- VlnPlot(AITL7.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)
ggsave("./result/preprocess/AITL7/AITL7_qc.pdf",qcplot,width=30,height=30,units='cm')
pdf("./result/preprocess/AITL7/AITL7_density.pdf", onefile=TRUE)
for (i in c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) {
   plot(density(AITL7.obj@meta.data[[i]]))
}
dev.off()

AITL7.obj1 <- subset(AITL7.obj, subset = nFeature_RNA<6000 & nFeature_RNA > 200 & percent.mt<20) 
saveRDS(AITL7.obj1, "./result/preprocess/AITL7/AITL7.RDS")


####################################
# NC1
####################################
rm(list=ls())
dir.create("./result/preprocess/NC1")
NC1 <- read.table(gzfile("./expression_matrix/NC1.expression_matrix.txt.gz"), sep = "\t", row.names = 1, header = TRUE)
NC1.obj <- CreateSeuratObject(NC1, project="NC1", min.cells = 3)
NC1.obj$percent.mt <- PercentageFeatureSet(NC1.obj, pattern = "^MT-")
ribogenes <- union(rownames(NC1.obj)[str_detect(rownames(NC1.obj), "^RPS")],
	rownames(NC1.obj)[str_detect(rownames(NC1.obj), "^RPL")])
NC1.obj$percent.ribo <- PercentageFeatureSet(NC1.obj, feature = ribogenes)
qcplot <- VlnPlot(NC1.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)
ggsave("./result/preprocess/NC1/NC1_qc.pdf",qcplot,width=30,height=30,units='cm')
pdf("./result/preprocess/NC1/NC1_density.pdf", onefile=TRUE)
for (i in c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) {
   plot(density(NC1.obj@meta.data[[i]]))
}
dev.off()
NC1.obj1 <- subset(NC1.obj, subset = nFeature_RNA<6000 & nFeature_RNA > 200 & percent.mt<20) 
saveRDS(NC1.obj1, "./result/preprocess/NC1/NC1.RDS")


####################################
#NC2
####################################
rm(list=ls())
dir.create("./result/preprocess/NC2")
NC2 <- read.table(gzfile("./expression_matrix/NC2.expression_matrix.txt.gz"), sep = "\t", row.names = 1, header = TRUE)
NC2.obj <- CreateSeuratObject(NC2, project="NC2", min.cells = 3)
NC2.obj$percent.mt <- PercentageFeatureSet(NC2.obj, pattern = "^MT-")
ribogenes <- union(rownames(NC2.obj)[str_detect(rownames(NC2.obj), "^RPS")],
	rownames(NC2.obj)[str_detect(rownames(NC2.obj), "^RPL")])
NC2.obj$percent.ribo <- PercentageFeatureSet(NC2.obj, feature = ribogenes)
qcplot <- VlnPlot(NC2.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)
ggsave("./result/preprocess/NC2/NC2_qc.pdf",qcplot,width=30,height=30,units='cm')
pdf("./result/preprocess/NC2/NC2_density.pdf", onefile=TRUE)
for (i in c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) {
   plot(density(NC2.obj@meta.data[[i]]))
}
dev.off()
NC2.obj1 <- subset(NC2.obj, subset = nFeature_RNA<6000 & nFeature_RNA > 200 & percent.mt<20) 
saveRDS(NC2.obj1, "./result/preprocess/NC2/NC2.RDS")

####################################
#NC3
####################################
rm(list=ls())
dir.create("./result/preprocess/NC3")
NC3 <- read.table(gzfile("./expression_matrix/NC3.expression_matrix.txt.gz"), sep = "\t", row.names = 1, header = TRUE)
NC3.obj <- CreateSeuratObject(NC3, project="NC3", min.cells = 3)
NC3.obj$percent.mt <- PercentageFeatureSet(NC3.obj, pattern = "^MT-")
ribogenes <- union(rownames(NC3.obj)[str_detect(rownames(NC3.obj), "^RPS")],
	rownames(NC3.obj)[str_detect(rownames(NC3.obj), "^RPL")])
NC3.obj$percent.ribo <- PercentageFeatureSet(NC3.obj, feature = ribogenes)
qcplot <- VlnPlot(NC3.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)
ggsave("./result/preprocess/NC3/NC3_qc.pdf",qcplot,width=30,height=30,units='cm')
pdf("./result/preprocess/NC3/NC3_density.pdf", onefile=TRUE)
for (i in c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) {
   plot(density(NC3.obj@meta.data[[i]]))
}
dev.off()
NC3.obj1 <- subset(NC3.obj, subset = nFeature_RNA<6000 & nFeature_RNA > 200 & percent.mt<20) 
saveRDS(NC3.obj1, "./result/preprocess/NC3/NC3.RDS")


    



AITL1.obj <- RenameCells(object = AITL1.obj, add.cell.id = "AITL1")
AITL2.obj <- RenameCells(object = AITL2.obj, add.cell.id = "AITL2")
AITL3.obj <- RenameCells(object = AITL3.obj, add.cell.id = "AITL3")
AITL4.obj <- RenameCells(object = AITL4.obj, add.cell.id = "AITL4")
AITL7.obj <- RenameCells(object = AITL7.obj, add.cell.id = "AITL7")
NC1.obj <- RenameCells(object = NC1.obj, add.cell.id = "NC1")
NC2.obj <- RenameCells(object = NC2.obj, add.cell.id = "NC2")
NC3.obj <- RenameCells(object = NC3.obj, add.cell.id = "NC3")

merged <- merge(AITL1.obj, y = c(AITL2.obj, AITL3.obj, AITL4.obj, AITL7.obj, 
	NC1.obj, NC2.obj, NC3.obj), project = "AITL")
library(ggplot2)
library(ggridges)
mycolors = c("deepskyblue", "#46A040" ,"#FFC179","#7D7D7D" ,  "#98D9E9","#F6313E", 
    "turquoise1","maroon1" ,"#FF5A00" ,"#AFAFAF" , "#00AF99","#FFA300" ,"orangered" ,"#0081C9" , 
    "indianred1","#490C65" ,"#BA7FD0","darkgreen", "orange2","lightseagreen","royalblue1","red4","red","darkmagenta","lawngreen")
dat <- as.data.frame(merged@meta.data[, c("orig.ident", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")])
pdf("./result/preprocess/nCount_RNA.pdf")
ggplot(dat, aes(x = log2(nCount_RNA+1), y = orig.ident)) +
  	geom_density_ridges(aes(fill = orig.ident)) + 
	scale_fill_manual(values = mycolors)
dev.off()
pdf("./result/preprocess/nFeature_RNA.pdf")
ggplot(dat, aes(x = nFeature_RNA, y = orig.ident)) +
  	geom_density_ridges(aes(fill = orig.ident)) + 
	scale_fill_manual(values = mycolors)
dev.off()
pdf("./result/preprocess/MT.pdf")
ggplot(dat, aes(x = percent.mt, y = orig.ident)) +
  	geom_density_ridges(aes(fill = orig.ident)) + 
	scale_fill_manual(values = mycolors)
dev.off()
pdf("./result/preprocess/Ribo.pdf")
ggplot(dat, aes(x = percent.ribo, y = orig.ident)) +
  	geom_density_ridges(aes(fill = orig.ident)) + 
	scale_fill_manual(values = mycolors)
dev.off()










AITL1 <- readRDS("./result/preprocess/AITL1/AITL1.RDS")
AITL2 <- readRDS("./result/preprocess/AITL2/AITL2.RDS")
AITL3 <- readRDS("./result/preprocess/AITL3/AITL3.RDS")
AITL4 <- readRDS("./result/preprocess/AITL4/AITL4.RDS")
AITL7 <- readRDS("./result/preprocess/AITL7/AITL7.RDS")
NC1 <- readRDS("./result/preprocess/NC1/NC1.RDS")
NC2 <- readRDS("./result/preprocess/NC2/NC2.RDS")
NC3 <- readRDS("./result/preprocess/NC3/NC3.RDS")

AITL1 <- run_Seurat(AITL1)
AITL1 <- RenameCells(object = AITL1, add.cell.id = "AITL1")
AITL1 <- CellCycleScoring(object = AITL1, g2m.features = cc.genes$g2m.genes, 
	s.features = cc.genes$s.genes)
saveRDS(AITL1, "./result/preprocess/AITL1/AITL1.RDS")

AITL2 <- run_Seurat(AITL2)
AITL2 <- RenameCells(object = AITL2, add.cell.id = "AITL2")
AITL2 <- CellCycleScoring(object = AITL2, g2m.features = cc.genes$g2m.genes, 
	s.features = cc.genes$s.genes)
saveRDS(AITL2, "./result/preprocess/AITL2/AITL2.RDS")

LHL <- run_Seurat(LHL)
LHL <- RenameCells(object = LHL, add.cell.id = "LHL")
LHL <- CellCycleScoring(object = LHL, g2m.features = cc.genes$g2m.genes, 
	s.features = cc.genes$s.genes)
saveRDS(LHL, "./result/preprocess/LHL/LHL.RDS")

AITL3 <- run_Seurat(AITL3)
AITL3 <- RenameCells(object = AITL3, add.cell.id = "AITL3")
AITL3 <- CellCycleScoring(object = AITL3, g2m.features = cc.genes$g2m.genes, 
	s.features = cc.genes$s.genes)
saveRDS(AITL3, "./result/preprocess/AITL3/AITL3.RDS")

AITL4 <- run_Seurat(AITL4)
AITL4 <- RenameCells(object = AITL4, add.cell.id = "AITL4")
AITL4 <- CellCycleScoring(object = AITL4, g2m.features = cc.genes$g2m.genes, 
	s.features = cc.genes$s.genes)
saveRDS(AITL4, "./result/preprocess/AITL4/AITL4.RDS")

AITL7 <- run_Seurat(AITL7)
AITL7 <- RenameCells(object = AITL7, add.cell.id = "AITL7")
AITL7 <- CellCycleScoring(object = AITL7, g2m.features = cc.genes$g2m.genes, 
	s.features = cc.genes$s.genes)
saveRDS(AITL7, "./result/preprocess/AITL7/AITL7.RDS")

NC1 <- run_Seurat(NC1)
NC1 <- RenameCells(object = NC1, add.cell.id = "NC1")
NC1 <- CellCycleScoring(object = NC1, g2m.features = cc.genes$g2m.genes, 
	s.features = cc.genes$s.genes)
saveRDS(NC1, "./result/preprocess/NC1/NC1.RDS")

NC2 <- run_Seurat(NC2)
NC2 <- RenameCells(object = NC2, add.cell.id = "NC2")
NC2 <- CellCycleScoring(object = NC2, g2m.features = cc.genes$g2m.genes, 
	s.features = cc.genes$s.genes)
saveRDS(NC2, "./result/preprocess/NC2/NC2.RDS")

NC3 <- run_Seurat(NC3)
NC3 <- RenameCells(object = NC3, add.cell.id = "NC3")
NC3 <- CellCycleScoring(object = NC3, g2m.features = cc.genes$g2m.genes, 
	s.features = cc.genes$s.genes)
saveRDS(NC3, "./result/preprocess/NC3/NC3.RDS")





