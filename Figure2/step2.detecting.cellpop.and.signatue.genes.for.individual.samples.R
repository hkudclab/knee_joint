#step3.successive-TP.noCC.preprocessing.R

library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggplot2)

load("../seurat/KJ12/KJ125.first.RData")
load("../seurat/KJ13/KJ135.first.RData")
load("../seurat/KJ14/KJ145.first.RData")
load("../seurat/KJ15/KJ155.first.RData")
load("../seurat/KJ18/KJ185.first.RData")
load("../seurat/KJP5/KJP5.first.RData")
load("../seurat/MCP5/MCP5.third.RData")

load("../seurat/KJ12/KJ125.first.markers.RData")
load("../seurat/KJ13/KJ135.first.markers.RData")
load("../seurat/MCP5/MCP5.third.markers.RData")

load("../DATA/KJ.combCore.meta.RData")
CoreMetNonCycling<-droplevels(CoreMet[!CoreMet$seurat_clusters%in%c(3,4), ])
LIDs<-sapply(split(rownames(CoreMetNonCycling), CoreMetNonCycling$orig.ident), function(x)gsub("_.*","",x))

#################################################
tsne12<-read.csv("../DATA/KJ125/projection.csv")
tsne13<-read.csv("../DATA/KJ135/projection.csv")
tsne14<-read.csv("../DATA/KJ145/projection.csv")
tsne15<-read.csv("../DATA/KJ155/projection.csv")
tsne18<-read.csv("../DATA/KJ185/projection.csv")
tsneP5<-read.csv("../DATA/KJP5/projection.csv")
tsneMP5<-read.csv("../DATA/MCP5/projection.csv")
#################################################

newKJ<-function(embed, OBJ, i){
	indObj<-match(colnames(OBJ),embed[,1])
	OBJ@reductions$tsne@cell.embeddings[,1:2]<-as.matrix(embed[indObj,2:3])
	print(head(OBJ@reductions$tsne@cell.embeddings[,1:2]))
	CoreMeti<- CoreMet[grep(paste0("_",i),rownames(CoreMet)),]
	print(dim(CoreMeti))
	ROWNi<-gsub("_[0-9]$","",rownames(CoreMeti))
	print(table(ROWNi%in%colnames(OBJ)))
	print(table(colnames(OBJ)%in%ROWNi))
	OBJ<-subset(OBJ, cells=ROWNi)
	indMeta<-match(colnames(OBJ), ROWNi)
	METACOL<-setdiff(colnames(CoreMeti), colnames(OBJ@meta.data))
	print(METACOL)
	for(mc in METACOL){
		OBJ[[mc]]<-CoreMeti[[mc]][indMeta]
	}
	return(OBJ)
}

KJ125.2<-newKJ(tsne12, KJ125, 1)
KJ135.2<-newKJ(tsne13, KJ135, 2)
KJ145.2<-newKJ(tsne14, KJ145, 3)
KJ155.2<-newKJ(tsne15, KJ155, 4)
KJ185.2<-newKJ(tsne18, KJ185, 5)
KJP5.2<-newKJ(tsneP5, KJP5, 6)
MCP5.2<-MCP5_3# newKJ(tsneMP5, MCP5, 7)
#################################################
NOCC<-function(OBJ, CELLSin){
	print(OBJ)
	OBJ<-subset(OBJ, cells = CELLSin)
	cat("---------\n")
	print(OBJ)

	OBJ<-NormalizeData(OBJ)
	OBJ<-ScaleData(OBJ,fearures=rownames(OBJ))
	OBJ<-FindVariableFeatures(OBJ)
	OBJ<-RunPCA(OBJ)
	OBJ<-RunUMAP(OBJ, reduction= "pca", dims= 1:30)
#	OBJ<-RunTSNE(OBJ, reduction= "pca", dims.use = 1:40, do.fast = T)

	OBJ<-FindNeighbors(OBJ, dims = 1:10)
	OBJ<-FindClusters(OBJ, resolution = 0.5)
	return(OBJ)
}

KJ125.3<-NOCC(KJ125.2, LIDs[["KJ125"]])
KJ135.3<-NOCC(KJ135.2, LIDs[["KJ135"]])
KJ145.3<-NOCC(KJ145.2, LIDs[["KJ145"]])
KJ155.3<-NOCC(KJ155.2, LIDs[["KJ155"]])
KJ185.3<-NOCC(KJ185.2, LIDs[["KJ185"]])
KJP5.3<-NOCC(KJP5.2, LIDs[["KJP5"]])
MCP5.3<-MCP5.2


sample.P5<-RunCCA(KJP5.3, RenameCells(MCP5.3,new.names=gsub("-1","-2",colnames(MCP5.3)) ))
sample.P5[["P5"]]<-paste0(sample.P5@meta.data$orig.ident, "_", sample.P5@meta.data$seurat_clusters)
	sample.P5<-NormalizeData(sample.P5)
	sample.P5<-ScaleData(sample.P5,fearures=rownames(sample.P5))
	sample.P5<-FindVariableFeatures(sample.P5)
	sample.P5<-RunPCA(sample.P5)
	sample.P5<-RunUMAP(sample.P5, reduction= "cca", dims= 1:10)
	sample.P5<-RunTSNE(sample.P5, reduction= "cca", dims= 1:10)

p1<-UMAPPlot(KJ125.3, label=T, label.size=6, pt.size=4) + ggtitle("KJ12.5")
p2<-UMAPPlot(KJ135.3, label=T, label.size=6, pt.size=4) + ggtitle("KJ13.5")
p3<-UMAPPlot(KJ145.3, label=T, label.size=6, pt.size=4) + ggtitle("KJ14.5")
p4<-UMAPPlot(KJ155.3, label=T, label.size=6, pt.size=4) + ggtitle("KJ15.5")
p5<-UMAPPlot(KJ185.3, label=T, label.size=6, pt.size=4) + ggtitle("KJ18.5")
p6<-UMAPPlot(KJP5.3, label=T, label.size=6, pt.size=4) + ggtitle("KJP5.5")
p7<-UMAPPlot(MCP5.3, label=T, label.size=6, pt.size=4) + ggtitle("MCP5.5")
p8<-UMAPPlot(sample.P5, group.by="P5", label=T, label.size=6, pt.size=4) + ggtitle("KJ+MC P5")
TSNEPlot(sample.P5, group.by="P5", split.by="orig.ident", label=T, label.size=6, pt.size=4) + ggtitle("KJ+MC P5")
UMAPPlot(sample.P5, group.by="P5", split.by="orig.ident", label=T, label.size=6, pt.size=4) + ggtitle("KJ+MC P5")
PCAPlot(sample.P5, group.by="P5", split.by="orig.ident", label=T, label.size=6, pt.size=4) + ggtitle("KJ+MC P5")

pdf("individual.TSNE.UMAP.without.cycling.cells.pdf", width=8)
	plot(p1)
	plot(p2)
	plot(p3)
	plot(p4)
	plot(p5)
	plot(p6)
	plot(p7)
	length(unique(Idents(KJ125.3)))
	length(unique(Idents(KJ135.3)))
	length(unique(Idents(KJ145.3)))
	length(unique(Idents(KJ155.3)))
	length(unique(Idents(KJ185.3)))
	length(unique(Idents(KJP5.3)))
	length(unique(Idents(MCP5.3)))
	barplot(rep(1,11),  col=colorRampPalette(brewer.pal(9, "Dark2"))(11))
dev.off()


table(Idents(KJ125.3))
table(Idents(KJ135.3))
table(Idents(KJ145.3))
table(Idents(KJ155.3))
table(Idents(KJ185.3))
table(Idents(KJP5.3))

#save(KJ125.3, KJ135.3, file="Individual.samples.with.GFP.12.13.noCycling.RData")
#save(KJ135.3, KJ145.3, file="Individual.samples.with.GFP.13.14.noCycling.RData")
#save(KJ145.3, KJ155.3, file="Individual.samples.with.GFP.14.15.noCycling.RData")
#save(KJ155.3, KJ185.3, file="Individual.samples.with.GFP.15.18.noCycling.RData")
#save(KJ185.3, KJP5.3, file="Individual.samples.with.GFP.18.p5.noCycling.RData")
#save(KJ185.3, MCP5.3, file="Individual.samples.with.GFP.18.MCp5.noCycling.RData")

FeaturePlot(KJ125.3, features="Mki67", reduction="tsne", pt.size=4)

META125<-KJ125.3@meta.data
META135<-KJ135.3@meta.data
META145<-KJ145.3@meta.data
META155<-KJ155.3@meta.data
META185<-KJ185.3@meta.data
METAP5<-KJP5.3@meta.data

#save(META125, META135, 
	META145, META155,
	META185, METAP5,
	file="step3.noCC.clusterMemb.per.sample.RData")


#################################################
ind125<-sapply(KJ125.2@meta.data[,],class)%in% c("numeric", "integer")
ind135<-sapply(KJ135.2@meta.data[,],class)%in% c("numeric", "integer")
ind145<-sapply(KJ145.2@meta.data[,],class)%in% c("numeric", "integer")
ind155<-sapply(KJ155.2@meta.data[,],class)%in% c("numeric", "integer")
ind185<-sapply(KJ185.2@meta.data[,],class)%in% c("numeric", "integer")
indp5<-sapply(KJP5.2@meta.data[,],class)%in% c("numeric", "integer")


matN<-as.matrix(table(c(colnames(KJ125.2@meta.data[,ind125]),colnames(KJ135.2@meta.data[,ind135]),
	colnames(KJ145.2@meta.data[,ind145]), colnames(KJ155.2@meta.data[,ind155]),
	colnames(KJ185.2@meta.data[,ind185]),colnames(KJP5.2@meta.data[,indp5]))))
matN<-as.matrix(matN[-grep("^RNA|seurat|^n|orig|nn_res|outlier",rownames(matN)),])
rowN2<-rownames(matN)

cbind(colSums(KJ125.2@meta.data[,rowN2]>0,na.rm=T),
	colSums(KJ135.2@meta.data[,rowN2]>0,na.rm=T),
	colSums(KJ145.2@meta.data[,rowN2]>0,na.rm=T),
	colSums(KJ155.2@meta.data[,rowN2]>0,na.rm=T),
	colSums(KJ185.2@meta.data[,rowN2]>0,na.rm=T),
	colSums(KJP5.2@meta.data[,rowN2]>0,na.rm=T))

#################################################
GETSIG<-function(MARKERS){
	sapply(levels(MARKERS$cluster),function(i){
		thisDEG<-MARKERS[MARKERS$cluster==i&MARKERS$p_val_adj<0.05,]
		thisDEG<-thisDEG[(thisDEG[,3]-thisDEG[,4])>0.25 | thisDEG[,2]>0.5 ,] #
		thisDEG<-thisDEG[order(thisDEG[,3]-thisDEG[,4],decreasing=T),]
		as.character(thisDEG$gene)
	})
}
#markers.KJ125.2 <- FindAllMarkers(KJ125.2, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
LSIG125_12.2<-GETSIG(markers.KJ125.2)

#################################################
#markers.KJ125.3 <- FindAllMarkers(KJ125.3, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
#markers.KJ135.3 <- FindAllMarkers(KJ135.3, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
#markers.KJ145.3 <- FindAllMarkers(KJ145.3, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
#markers.KJ155.3 <- FindAllMarkers(KJ155.3, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
#markers.KJ185.3 <- FindAllMarkers(KJ185.3, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
#markers.KJP5.3 <- FindAllMarkers(KJP5.3, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)



LSIG125.3<-GETSIG(markers.KJ125.3)
LSIG135.3<-GETSIG(markers.KJ135.3)
LSIG145.3<-GETSIG(markers.KJ145.3)
LSIG155.3<-GETSIG(markers.KJ155.3)
LSIG185.3<-GETSIG(markers.KJ185.3)
LSIGP5.3<-GETSIG(markers.KJP5.3)

length(LSIG125.3)*length(LSIG135.3)*length(LSIG145.3)*length(LSIG155.3)*length(LSIG185.3)*length(LSIGP5.3)


save(markers.KJ125.3, markers.KJ135.3, markers.KJ145.3, 
	markers.KJ155.3, markers.KJ185.3, markers.KJP5.3, 
	LSIG125.3, LSIG135.3, LSIG145.3, 
	LSIG155.3, LSIG185.3, LSIGP5.3, file="Markers.Individual.samples.with.GFP.noCycling.reso.0.5.RData")

sapply(LSIG125.3, function(x)x[1:40])
sapply(LSIG135.3, function(x)x[1:40])
sapply(LSIG145.3, function(x)x[1:40])
sapply(LSIG155.3, function(x)x[1:40])
sapply(LSIG185.3, function(x)x[1:40])
sapply(LSIGP5.3, function(x)x[1:40])

x<-c(as.vector(sapply(LSIG125.3, function(x)x[1:40])),
	as.vector(sapply(LSIG135.3, function(x)x[1:40])),
	as.vector(sapply(LSIG145.3, function(x)x[1:40])),
	as.vector(sapply(LSIG155.3, function(x)x[1:40])),
	as.vector(sapply(LSIG185.3, function(x)x[1:40])),
	as.vector(sapply(LSIGP5.3, function(x)x[1:40])),
	as.vector(sapply(LSIGP5.3, function(x)x[1:40])))
tabx<-sort(table(x[!is.na(x)]))


#################################################
pdf("New.SeuratClusters.on.Old.TSNE.all.cells.separate.sample2.pdf",width=7.5)
	TSNEPlot(KJ125.2, pt.size = 2,label=T,label.size=12)
	TSNEPlot(KJ135.2, pt.size = 2,label=T,label.size=12)
	TSNEPlot(KJ145.2, pt.size = 2,label=T,label.size=12)
	TSNEPlot(KJ155.2, pt.size = 2,label=T,label.size=12)
	TSNEPlot(KJ185.2, pt.size = 2,label=T,label.size=12)
	TSNEPlot(KJP5.2, pt.size = 2,label=T,label.size=12)
	TSNEPlot(MCP5.3, pt.size = 2,label=T,label.size=12)

	FeaturePlot(KJ125.2, features=c("GDF5pos", "Osr2"),
		pt.size = 2,col=c("grey","#28B34B"), reduction="tsne",
		label=T,label.size=12)
dev.off()

pdf("Dot.plot.for.Fig2D.pdf",width=8.5)
	FEATURES<-c("Gdf5", "Osr2", "Creb5", "Mia", "Col2a1", "Epyc", "Ucma", "Col22a1", 
		"Cxcl12", "Lrrc17", "Col1a1", "Scx", "Tnmd", "Postn", "Cilp", "Prg4")
	DotPlot(KJ125.3, features=FEATURES, cols=c("grey","brown"))
	DotPlot(KJ135.3, features=FEATURES, cols=c("grey","brown"))
	DotPlot(KJ145.3, features=FEATURES, cols=c("grey","brown"))
	DotPlot(KJ155.3, features=FEATURES, cols=c("grey","brown"))
	DotPlot(KJ185.3, features=FEATURES, cols=c("grey","brown"))
	DotPlot(KJP5.3, features=FEATURES, cols=c("grey","brown"))
	DotPlot(MCP5.3, features=FEATURES, cols=c("grey","brown"))

dev.off()

#################################################
FeaturePlot(sample.P5,label=T, label.size=6, features=c("Gdf5","Lgr5","Epyc","Col1a1", "Prg4", "Col22a1"),
		cols=c("lightgrey","brown"),split.by="orig.ident", by.col=F,
		reduction="umap",ncol=4,pt.size=1)

FeaturePlot(KJ125.2,label=T, features=c("Gdf5","Lgr5","EGFP","Col2a1","Col1a1","IRES","LACZ","gpCAG.CreERT2","CREERT2Nelson"),
		cols=c("lightgrey","brown"),reduction="tsne",ncol=4,pt.size=1)


TSNEPlot(KJ125.2, pt.size = 2,label=T,label.size=12)+
	UMAPPlot(KJ125, pt.size = 2,label=T,label.size=12)

TSNEPlot(KJ135, pt.size = 2,label=T,label.size=12)+
	UMAPPlot(KJ135, pt.size = 2,label=T,label.size=12)


TSNEPlot(KJP5.2, pt.size = 2,label=T,label.size=12)+
	UMAPPlot(KJP5, pt.size = 2,label=T,label.size=12)



