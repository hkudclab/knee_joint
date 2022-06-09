library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gplots)
library(ggalluvial)

load("../DATA/Individual.samples.with.GFP.18.P5.RData")
#########################################

#KJ1852<-subset(KJ185,cells=gsub("_.*","",rownames(CoreMet)[CoreMet$orig.ident=="KJ185"]))
#KJP52<-subset(KJP5,cells=gsub("_.*","",rownames(CoreMet)[CoreMet$orig.ident=="KJP5"]))
KJP5.4<-RenameCells(KJP5.3,new.names=gsub("-1","-2",colnames(KJP5.3)) )
set.seed(0)
KJ18P5<-RunCCA(KJ185.3,KJP5.4)


KJ18P5[["OldClust"]]<-paste0(KJ18P5@meta.data$orig.ident,"_",KJ18P5@meta.data$seurat_clusters)
table(KJ18P5[["OldClust"]][,1], KJ18P5@meta.data$orig.ident)


KJ18P5<-NormalizeData(KJ18P5)
KJ18P5<-ScaleData(KJ18P5,fearures=rownames(KJ18P5))
KJ18P5<-FindVariableFeatures(KJ18P5)
KJ18P5<-RunPCA(KJ18P5)
KJ18P5<-RunUMAP(KJ18P5, reduction= "cca", dims= 1:10)
KJ18P5<-RunTSNE(KJ18P5, reduction= "cca", dims.use = 1:40, do.fast = T)
KJ18P5<-FindNeighbors(KJ18P5, reduction= "cca", dims = 1:10)
KJ18P5<-FindClusters(KJ18P5, resolution = 0.25)
KJ18P5.5<-FindClusters(KJ18P5, resolution = 0.5)

table(KJ18P5[["OldClust"]][,1], Idents(KJ18P5))
table(KJ18P5.5[["OldClust"]][,1], Idents(KJ18P5.5))


#save(KJ18P5.5, file="successiveTP/successiveTP.18.P5.RData")

#set.seed(0)
#KJ18P5s<-subset(KJ18P5,cells=sample(colnames(KJ18P5),5e3))
#pdf("tsne.mKJ.18.P5.pdf")
#	TSNEPlot(KJ18P5,label=T,label.size=12, pt.size = 5)+ NoLegend()
#	TSNEPlot(KJ18P5s,label=T,label.size=12, pt.size = 5, split.by="orig.ident")+ NoLegend()
#
#dev.off()


#KJ18P5@reductions$tsne@cell.embeddings


ind185<-which(KJ18P5@meta.data$orig.ident=="KJ185")

MATdist18P5<-as.matrix(dist(KJ18P5.5@reductions$tsne@cell.embeddings[1:5,]))
ALLCLUST<-as.character(sort(unique(KJ18P5.5@meta.data$seurat_clusters)))

MATflow<-matrix(NA,0,3)
colnames(MATflow)<-c("i", "OldClustThisTP", "FLOWto")

for(i in 1:length(ALLCLUST)){
	indi<-which(KJ18P5.5@meta.data$seurat_clusters==ALLCLUST[i])
	MATdist18P5.i<-as.matrix(dist(KJ18P5.5@reductions$tsne@cell.embeddings[indi,]))
	indi2<-which(KJ18P5.5@meta.data$orig.ident[indi]=="KJ185")
	OldClustThisTP<-KJ18P5.5[["OldClust"]][indi, 1][indi2]
	OldClustNextTP<-KJ18P5.5[["OldClust"]][indi, 1][-indi2]

	MATdist18P5.i2<-MATdist18P5.i[indi2, -indi2]
	FLOWto<-apply(MATdist18P5.i2, 1, function(x){
		ORD3<-order(x)[1:3]
		names(which.max(table(OldClustNextTP[ORD3])))
	})
	MATflow<-rbind(MATflow, cbind(rep(i, length(OldClustThisTP)), OldClustThisTP, FLOWto))
}
save(MATflow, file="successiveTP/MATflow.E18.to.EP5.RData")
###############################
library(RColorBrewer)



pdf("successiveTP/alluvial.chart.18P5.pdf")
	ggplot(as.data.frame(table(MATflow[,2], MATflow[,3])),
      	 aes(y = Freq, axis1 = Var1, axis2 = Var2)) +
	  geom_alluvium(aes(fill = Var1), width = 1/12) +
	  geom_stratum(width = 1/12, fill = "grey", color = "white") +
	  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
	  scale_fill_discrete( colorRampPalette(brewer.pal(9, "Set1"))(length(table(MATflow[,2])))) 
dev.off()




