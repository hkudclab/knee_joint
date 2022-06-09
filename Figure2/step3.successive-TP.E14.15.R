library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gplots)
library(ggalluvial)

load("../DATA/Individual.samples.with.GFP.14.15.RData")
#########################################

#KJ1452<-subset(KJ145,cells=gsub("_.*","",rownames(CoreMet)[CoreMet$orig.ident=="KJ145"]))
#KJ1552<-subset(KJ155,cells=gsub("_.*","",rownames(CoreMet)[CoreMet$orig.ident=="KJ155"]))
KJ155.4<-RenameCells(KJ155.3,new.names=gsub("-1","-2",colnames(KJ155.3)) )
KJ1415<-RunCCA(KJ145.3,KJ155.4)


KJ1415[["OldClust"]]<-paste0(KJ1415@meta.data$orig.ident,"_",KJ1415@meta.data$seurat_clusters)
table(KJ1415[["OldClust"]][,1], KJ1415@meta.data$orig.ident)


KJ1415<-NormalizeData(KJ1415)
KJ1415<-ScaleData(KJ1415,fearures=rownames(KJ1415))
KJ1415<-FindVariableFeatures(KJ1415)
KJ1415<-RunPCA(KJ1415)
KJ1415<-RunUMAP(KJ1415, reduction= "cca", dims= 1:10)
KJ1415<-RunTSNE(KJ1415, reduction= "cca", dims.use = 1:40, do.fast = T)
KJ1415<-FindNeighbors(KJ1415, reduction= "cca", dims = 1:10)
KJ1415<-FindClusters(KJ1415, resolution = 0.25)
KJ1415.5<-FindClusters(KJ1415, resolution = 0.5)

table(KJ1415[["OldClust"]][,1], Idents(KJ1415))
table(KJ1415.5[["OldClust"]][,1], Idents(KJ1415.5))


#save(KJ1415.5, file="successiveTP/successiveTP.14.15.RData")

#set.seed(0)
#KJ1415s<-subset(KJ1415,cells=sample(colnames(KJ1415),5e3))
#pdf("tsne.mKJ.14.15.pdf")
#	TSNEPlot(KJ1415,label=T,label.size=12, pt.size = 5)+ NoLegend()
#	TSNEPlot(KJ1415s,label=T,label.size=12, pt.size = 5, split.by="orig.ident")+ NoLegend()
#
#dev.off()


#KJ1415@reductions$tsne@cell.embeddings


ind145<-which(KJ1415@meta.data$orig.ident=="KJ145")

MATdist1415<-as.matrix(dist(KJ1415.5@reductions$tsne@cell.embeddings[1:5,]))
ALLCLUST<-as.character(sort(unique(KJ1415.5@meta.data$seurat_clusters)))

MATflow<-matrix(NA,0,3)
colnames(MATflow)<-c("i", "OldClustThisTP", "FLOWto")

for(i in 1:length(ALLCLUST)){
	indi<-which(KJ1415.5@meta.data$seurat_clusters==ALLCLUST[i])
	MATdist1415.i<-as.matrix(dist(KJ1415.5@reductions$tsne@cell.embeddings[indi,]))
	indi2<-which(KJ1415.5@meta.data$orig.ident[indi]=="KJ145")
	OldClustThisTP<-KJ1415.5[["OldClust"]][indi, 1][indi2]
	OldClustNextTP<-KJ1415.5[["OldClust"]][indi, 1][-indi2]

	MATdist1415.i2<-MATdist1415.i[indi2, -indi2]
	FLOWto<-apply(MATdist1415.i2, 1, function(x){
		ORD3<-order(x)[1:3]
		names(which.max(table(OldClustNextTP[ORD3])))
	})
	MATflow<-rbind(MATflow, cbind(rep(i, length(OldClustThisTP)), OldClustThisTP, FLOWto))
}
save(MATflow, file="successiveTP/MATflow.E14.to.E15.RData")
###############################
pdf("successiveTP/alluvial.chart.1415.pdf")
	ggplot(as.data.frame(table(MATflow[,2], MATflow[,3])),
      	 aes(y = Freq, axis1 = Var1, axis2 = Var2)) +
	  geom_alluvium(aes(fill = Var1), width = 1/12) +
	  geom_stratum(width = 1/12, fill = "grey", color = "white") +
	  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
	  scale_fill_brewer(type = "qual", palette = "Set1") 
dev.off()




