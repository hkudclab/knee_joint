library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gplots)
library(ggalluvial)

load("../DATA/Individual.samples.with.GFP.15.18.RData")
#########################################

#KJ1552<-subset(KJ155,cells=gsub("_.*","",rownames(CoreMet)[CoreMet$orig.ident=="KJ155"]))
#KJ1852<-subset(KJ185,cells=gsub("_.*","",rownames(CoreMet)[CoreMet$orig.ident=="KJ185"]))
KJ185.4<-RenameCells(KJ185.3,new.names=gsub("-1","-2",colnames(KJ185.3)) )
set.seed(0)
KJ1518<-RunCCA(KJ155.3,KJ185.4)


KJ1518[["OldClust"]]<-paste0(KJ1518@meta.data$orig.ident,"_",KJ1518@meta.data$seurat_clusters)
table(KJ1518[["OldClust"]][,1], KJ1518@meta.data$orig.ident)


KJ1518<-NormalizeData(KJ1518)
KJ1518<-ScaleData(KJ1518,fearures=rownames(KJ1518))
KJ1518<-FindVariableFeatures(KJ1518)
KJ1518<-RunPCA(KJ1518)
KJ1518<-RunUMAP(KJ1518, reduction= "cca", dims= 1:10)
KJ1518<-RunTSNE(KJ1518, reduction= "cca", dims.use = 1:40, do.fast = T)
KJ1518<-FindNeighbors(KJ1518, reduction= "cca", dims = 1:10)
KJ1518<-FindClusters(KJ1518, resolution = 0.25)
KJ1518.5<-FindClusters(KJ1518, resolution = 0.5)

table(KJ1518[["OldClust"]][,1], Idents(KJ1518))
table(KJ1518.5[["OldClust"]][,1], Idents(KJ1518.5))


#save(KJ1518.5, file="successiveTP/successiveTP.15.18.RData")

#set.seed(0)
#KJ1518s<-subset(KJ1518,cells=sample(colnames(KJ1518),5e3))
#pdf("tsne.mKJ.15.18.pdf")
#	TSNEPlot(KJ1518,label=T,label.size=12, pt.size = 5)+ NoLegend()
#	TSNEPlot(KJ1518s,label=T,label.size=12, pt.size = 5, split.by="orig.ident")+ NoLegend()
#
#dev.off()


#KJ1518@reductions$tsne@cell.embeddings


ind155<-which(KJ1518@meta.data$orig.ident=="KJ155")

MATdist1518<-as.matrix(dist(KJ1518.5@reductions$tsne@cell.embeddings[1:5,]))
ALLCLUST<-as.character(sort(unique(KJ1518.5@meta.data$seurat_clusters)))

MATflow<-matrix(NA,0,3)
colnames(MATflow)<-c("i", "OldClustThisTP", "FLOWto")

for(i in 1:length(ALLCLUST)){
	indi<-which(KJ1518.5@meta.data$seurat_clusters==ALLCLUST[i])
	MATdist1518.i<-as.matrix(dist(KJ1518.5@reductions$tsne@cell.embeddings[indi,]))
	indi2<-which(KJ1518.5@meta.data$orig.ident[indi]=="KJ155")
	OldClustThisTP<-KJ1518.5[["OldClust"]][indi, 1][indi2]
	OldClustNextTP<-KJ1518.5[["OldClust"]][indi, 1][-indi2]

	MATdist1518.i2<-MATdist1518.i[indi2, -indi2]
	FLOWto<-apply(MATdist1518.i2, 1, function(x){
		ORD3<-order(x)[1:3]
		names(which.max(table(OldClustNextTP[ORD3])))
	})
	MATflow<-rbind(MATflow, cbind(rep(i, length(OldClustThisTP)), OldClustThisTP, FLOWto))
}
save(MATflow, file="successiveTP/MATflow.E15.to.E18.RData")
###############################
library(RColorBrewer)



pdf("successiveTP/alluvial.chart.1518.pdf")
	ggplot(as.data.frame(table(MATflow[,2], MATflow[,3])),
      	 aes(y = Freq, axis1 = Var1, axis2 = Var2)) +
	  geom_alluvium(aes(fill = Var1), width = 1/12) +
	  geom_stratum(width = 1/12, fill = "grey", color = "white") +
	  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
	  scale_fill_discrete( colorRampPalette(brewer.pal(9, "Set1"))(length(table(MATflow[,2])))) 
dev.off()




