library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(igraph)
library(scales)

library(monocle3)

library(slingshot, quietly = FALSE)
library(RColorBrewer)

load("./DATA/monocle3.data.format.RData")
load("./DATA/Individual.samples.with.GFP.12.13.noCycling.RData")

load("./resources/Ligands.mouse.human.gene-symb.RData")
load("./resources/matrisome-mit-mouse/Matrisome.mouse.RData")


#########################################
NAMES<-c("0"="JPC/Chon-prog","1"="Lig","2"="Chon-prog","3"="Cyc","4"="Cyc","5"="JPC-Lig","6"="Lig-prog","7"="GPC-prog","8"="Mes-EMT")
KJ125.3[["BioName"]]<-NAMES[KJ125.3@meta.data$seurat_clusters]
 table(KJ125.3[["BioName"]][,1], KJ125.3@meta.data$seurat_clusters)

KJ125.4<-subset(KJ125.3, cells=colnames(KJ125.3)[KJ125.3@reductions$umap@cell.embeddings[,2]> -4])
UMAPPlot(KJ125.3, pt.size = 2,group.by="BioName", label=T,label.size=6) + NoLegend()
UMAPPlot(KJ125.4, pt.size = 2,group.by="BioName", label=T,label.size=6) + NoLegend()

tagLgr<-table(GetAssayData(KJ125.4)["Lgr5",]>0, KJ125.4@meta.data$BioName)
par(mar=c(8,4,4,2))
barplot(100*tagLgr[2,]/colSums(tagLgr), ylab="Lgr5+ %", las=2)


big_cds<-cds.125[,colnames(KJ125.4)]

table(colnames(big_cds)==colnames(KJ125.4))

METASZ<-KJ125.4@meta.data

big_cds@colData[["orig.ident"]]<-METASZ[["orig.ident"]]
big_cds@colData[["nCount_RNA"]]<-METASZ[["nCount_RNA"]]
big_cds@colData[["nFeature_RNA"]]<-METASZ[["nFeature_RNA"]]
big_cds@colData[["seurat_clusters"]]<-METASZ[["seurat_clusters"]]
big_cds@colData[["CREERT_nelson35"]]<-METASZ[["CREERT_nelson35"]]
big_cds@colData[["BioName"]]<-METASZ[["BioName"]]


#markers.KJ125 <- FindAllMarkers(KJ125.4, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
Idents(KJ125.4)<-"BioName"
FeaturePlot(KJ125.4, features=c("Top2a","Mki67","Birc5"), label=T)
FeaturePlot(KJ125.4, features=c("Epyc","Col2a1","Sox9", "Acan"), label=T, label.size=6, pt.size=1.5, cols=c("grey","brown"))

FeaturePlot(KJ125.4, features=c("Creb3l1","Prg4","Sod3", "Cilp2", "Scx", "Mkx", "Lrrc17"), label=T, label.size=6, pt.size=1, cols=c("grey","brown"))
FeaturePlot(KJ125.4, features=c( "Scx", "Mkx", "Lrrc17", "Gdf5"), label=T, label.size=6, pt.size=1.5, cols=c("grey","brown"))
FeaturePlot(KJ125.4, features=c( "Twist1", "Col1a1", "Col9a1", "Foxp1", "Gdf5", "Creb5"), label=T, label.size=5, pt.size=1, cols=c("grey","brown"), nco=3)


LSIG125<-sapply(levels(markers.KJ125$cluster),function(i){
	thisDEG<-markers.KJ125[markers.KJ125$cluster==i&markers.KJ125$p_val_adj<0.05,]
	thisDEG<-thisDEG[(thisDEG[,3]-thisDEG[,4])>0.25& thisDEG[,2]>0.5 ,] #
	thisDEG<-thisDEG[order(thisDEG[,3]-thisDEG[,4],decreasing=T),]
	as.character(thisDEG$gene)
})
names(LSIG125)<-NAMES[names(LSIG125)]


LSIG125


pdf("step9-kj12.5/feature.plots.pdf", height=8, width=9)
	GENES<-c( "Twist1", "Col1a1", "Col9a1", "Creb5", "Cxcl12", "Foxp1","Birc5", "Mki67", "Gdf5", "Creb5",
		"Osr1", "Prrx1", "Acan", "Cnmd", "Ihh", "Runx2", "Ptch1", "Nog", "Dkk2")
	for(thisG in GENES){
		print(thisG )
		p1<-FeaturePlot(KJ125.4, features=thisG , label=T, label.size=5, pt.size=1, cols=c("grey","brown"))
		plot(p1)
	}
dev.off()

########################
	UMAP<-as.data.frame(KJ125.4@reductions$umap@cell.embeddings)
	thisExpr<-GetAssayData(KJ125.4)["Gdf5",]
		NumStep<-60
		xStep<-seq(min(UMAP$UMAP_1), max(UMAP$UMAP_1), len=NumStep)
		yStep<-seq(min(UMAP$UMAP_2), max(UMAP$UMAP_2), len=NumStep)
		IMGmat<-matrix(0, NumStep, NumStep)
		GridEntropy<-IMGmat
		for(i in 2:NumStep){
			for(j in 2:NumStep){
				indi<-which(UMAP$UMAP_1>xStep[i-1] & UMAP$UMAP_1< xStep[i])
				indj<-which(UMAP$UMAP_2>yStep[j-1] & UMAP$UMAP_2< yStep[j])
				indij<-intersect(indi, indj)
				IMGmat[i,j]<-length(indij)
				if(length(indij)>0){
					GridEntropy[i,j]<-mean(thisExpr[indij])
				}
			}
		}
		image(GridEntropy)

####


PCA<-KJ125.4@reductions$pca@cell.embeddings

lin1 <- getLineages(PCA, KJ125.4@meta.data$BioName, omega = TRUE, start.clus = 'JPC/Chon-prog')
plot(PCA, col = brewer.pal(9,"Set1")[as.factor(KJ125.4@meta.data$BioName)], asp = 1, pch = 16)
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')


########################

big_cds <- preprocess_cds(big_cds, num_dim = 20)
#big_cds<- align_cds(big_cds, alignment_group = "orig.ident")

plot_pc_variance_explained(big_cds)
barplot(apply(big_cds@int_colData@listData$reducedDims@listData$PCA, 2, var))
big_cds<- reduce_dimension(big_cds)

big_cds<- cluster_cells(big_cds)#, reduction_method="UMAP", resolution=2)
#big_cds<- cluster_cells(big_cds, reduction_method="PCA", resolution=2)
#big_cds<- cluster_cells(big_cds, reduction_method="Aligned")

big_cds@int_colData@listData$reducedDims@listData$UMAP <- big_cds@int_colData@listData$reducedDims$PCA[,c(1,2)]

big_cds@int_colData@listData$reducedDims@listData$UMAP <- KJ125.4@reductions$umap@cell.embeddings
big_cds<- learn_graph(big_cds, use_partition = F, close_loop = F) #, learn_graph_control=list(ncenter=8))

big_cds2<-big_cds
big_cds2@int_colData@listData$reducedDims@listData$UMAP <- KJ125.4@reductions$pca@cell.embeddings
big_cds2<- learn_graph(big_cds2, use_partition = F, close_loop = F) #, learn_graph_control=list(ncenter=8))

pdf("step9-kj12.5/pseudotime-and-slingshot-trajectory.pdf", width=9, height=8)
	plot_cells(big_cds, color_cells_by="BioName", reduction_method="UMAP", label_groups_by_cluster=FALSE,  cell_size =2, group_label_size = 6,  graph_label_size=5)

	plot_cells(big_cds2, color_cells_by="BioName", reduction_method="UMAP", label_groups_by_cluster=FALSE,  cell_size =2, group_label_size = 6,  graph_label_size=5)


	M1<-sapply(split(PCA[,1], KJ125.4@meta.data$BioName), mean)
	M2<-sapply(split(PCA[,2], KJ125.4@meta.data$BioName), mean)
	PCA<-KJ125.4@reductions$pca@cell.embeddings[,1:15]
	lin1 <- getLineages(PCA, KJ125.4@meta.data$BioName, omega = TRUE, start.clus = 'JPC/Chon-prog')
	plot(PCA, col =  hue_pal()(8)[as.factor(KJ125.4@meta.data$BioName)], asp = 1, pch = 16)
	lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')
	text(M1,M2, names(M1), cex=2, col="red")

	#crv1 <- getCurves(lin1)
	lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')


dev.off()

pdf("step9-kj12.5/FeaturePlot-PCA.pdf", width=10, height=8)
	FeaturePlot(KJ125.4, features="Lgr5", reduction = "pca", cols=c("grey", "brown"), pt.size = 2,label=T,label.size=6) + NoLegend()
	FeaturePlot(KJ125.4, features="LGRpos_CREERT_pos", reduction = "pca", cols=c("grey", "brown"), pt.size = 2,label=T,label.size=6) + NoLegend()
	FeaturePlot(KJ125.4, features="CREERT_nelson35", reduction = "pca", cols=c("grey", "brown"), pt.size = 2,label=T,label.size=6) + NoLegend()
	FeaturePlot(KJ125.4, features="Gdf5", reduction = "pca", cols=c("grey", "brown"), pt.size = 2,label=T,label.size=6) + NoLegend()
	FeaturePlot(KJ125.4, features="Sox9", reduction = "pca", cols=c("grey", "brown"), pt.size = 2,label=T,label.size=6) + NoLegend()
	FeaturePlot(KJ125.4, features="Scx", reduction = "pca", cols=c("grey", "brown"), pt.size = 2,label=T,label.size=6) + NoLegend()
	FeaturePlot(KJ125.4, features="Creb5", reduction = "pca", cols=c("grey", "brown"), pt.size = 2,label=T,label.size=6) + NoLegend()
	FeaturePlot(KJ125.4, features="Wnt5a", reduction = "pca", cols=c("grey", "brown"), pt.size = 2,label=T,label.size=6) + NoLegend()
	FeaturePlot(KJ125.4, features="Twist1", reduction = "pca", cols=c("grey", "brown"), pt.size = 2,label=T,label.size=6) + NoLegend()
	FeaturePlot(KJ125.4, features="Prrx1", reduction = "pca", cols=c("grey", "brown"), pt.size = 2,label=T,label.size=6) + NoLegend()
	FeaturePlot(KJ125.4, features="Col1a1", reduction = "pca", cols=c("grey", "brown"), pt.size = 2,label=T,label.size=6) + NoLegend()
	FeaturePlot(KJ125.4, features="Col2a1", reduction = "pca", cols=c("grey", "brown"), pt.size = 2,label=T,label.size=6) + NoLegend()
dev.off()


plot_cells(big_cds, genes="Top2a", color_cells_by="BioName", reduction_method="UMAP", label_groups_by_cluster=FALSE,  cell_size =2, group_label_size = 6,  graph_label_size=5)








