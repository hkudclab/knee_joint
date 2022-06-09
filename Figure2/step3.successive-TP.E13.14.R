library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gplots)
library(ggalluvial)

load("../DATA/Individual.samples.with.GFP.13.14.RData")
#########################################

#KJ1352<-subset(KJ135,cells=gsub("_.*","",rownames(CoreMet)[CoreMet$orig.ident=="KJ135"]))
#KJ1452<-subset(KJ145,cells=gsub("_.*","",rownames(CoreMet)[CoreMet$orig.ident=="KJ145"]))
KJ145.4<-RenameCells(KJ145.3,new.names=gsub("-1","-2",colnames(KJ145.3)) )
set.seed(0)
KJ1314<-RunCCA(KJ135.3,KJ145.4)


KJ1314[["OldClust"]]<-paste0(KJ1314@meta.data$orig.ident,"_",KJ1314@meta.data$seurat_clusters)
table(KJ1314[["OldClust"]][,1], KJ1314@meta.data$orig.ident)


KJ1314<-NormalizeData(KJ1314)
KJ1314<-ScaleData(KJ1314,fearures=rownames(KJ1314))
KJ1314<-FindVariableFeatures(KJ1314)
KJ1314<-RunPCA(KJ1314)
KJ1314<-RunUMAP(KJ1314, reduction= "cca", dims= 1:10)
KJ1314<-RunTSNE(KJ1314, reduction= "cca", dims.use = 1:40, do.fast = T)
KJ1314<-FindNeighbors(KJ1314, reduction= "cca", dims = 1:10)
KJ1314<-FindClusters(KJ1314, resolution = 0.25)
KJ1314.5<-FindClusters(KJ1314, resolution = 0.5)

table(KJ1314[["OldClust"]][,1], Idents(KJ1314))
table(KJ1314.5[["OldClust"]][,1], Idents(KJ1314.5))


#save(KJ1314, file="successiveTP/successiveTP.13.14.RData")
#save(KJ1314.5, file="successiveTP/successiveTP.13.14.RData")

#set.seed(0)
#KJ1314s<-subset(KJ1314,cells=sample(colnames(KJ1314),5e3))
#pdf("tsne.mKJ.13.14.pdf")
#	TSNEPlot(KJ1314,label=T,label.size=12, pt.size = 5)+ NoLegend()
#	TSNEPlot(KJ1314s,label=T,label.size=12, pt.size = 5, split.by="orig.ident")+ NoLegend()
#
#dev.off()


#KJ1314@reductions$tsne@cell.embeddings


MATdist1314<-as.matrix(dist(KJ1314.5@reductions$tsne@cell.embeddings[1:5,]))
ALLCLUST<-as.character(sort(unique(KJ1314.5@meta.data$seurat_clusters)))

MATflow<-matrix(NA,0,5)
colnames(MATflow)<-c("i", "cellIDfrom", "cellIDto", "OldClustThisTP", "FLOWto")

LIDto<-list()
LIDto2<-list()

for(i in 1:length(ALLCLUST)){
	indi<-which(KJ1314.5@meta.data$seurat_clusters==ALLCLUST[i])
	MATdist1314.i<-as.matrix(dist(KJ1314.5@reductions$tsne@cell.embeddings[indi,]))
	indi2<-which(KJ1314.5@meta.data$orig.ident[indi]=="KJ135")
	indi2.next<-which(KJ1314.5@meta.data$orig.ident[indi]!="KJ135")
	OldClustThisTP<-KJ1314.5[["OldClust"]][indi, 1][indi2]
	OldClustNextTP<-KJ1314.5[["OldClust"]][indi, 1][indi2.next]

	if(length(indi2)>0 & length(indi2.next)>0){
		MATdist1314.i2<-MATdist1314.i[indi2, indi2.next]
		CELLidsTO<-colnames(KJ1314.5)[indi][indi2.next]
		FLOWto<-apply(MATdist1314.i2, 1, function(x){
			ORD3<-order(x)[1:3]
			names(which.max(table(OldClustNextTP[ORD3])))
		})
		ALLOCATED<-rep(0, ncol(MATdist1314.i2))
		FLOWtoCell<-sapply(1:nrow(MATdist1314.i2), function(j){
			x<-MATdist1314.i2[j, ]
			indUnallocated<-which(ALLOCATED==0)
			if(length(indUnallocated)==0)indUnallocated<-seq_along(ALLOCATED)
			ORD3<-order(x[indUnallocated])[1:3]
			flowToCluster<-names(which.max(table(OldClustNextTP[indUnallocated][ORD3])))

			ord4<-ORD3[which(OldClustNextTP[indUnallocated][ORD3]==flowToCluster)[1]]
			ALLOCATED[indUnallocated[ord4]]<<-1

			CELLidsTO[indUnallocated[ord4]]
		})
		LIDto[[i]]<-FLOWtoCell
		LIDto2[[i]]<-FLOWtoCell
		thisTemp<-cbind(i=rep(i, length(OldClustThisTP)), 
			cellIDfrom=colnames(KJ1314.5)[indi][indi2],
			cellIDto=FLOWtoCell,
			OldClustThisTP, FLOWto)
		MATflow<-rbind(MATflow, thisTemp)
	}
}


save(MATflow, file="successiveTP/MATflow.E13.to.E14.RData")
###############################
pdf("successiveTP/alluvial.chart.1314.pdf")
	ggplot(as.data.frame(table(MATflow[,4], MATflow[,5])),
      	 aes(y = Freq, axis1 = Var1, axis2 = Var2)) +
	  geom_alluvium(aes(fill = Var1), width = 1/12) +
	  geom_stratum(width = 1/12, fill = "grey", color = "white") +
	  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
	  scale_fill_brewer(type = "qual", palette = "Set1") 
dev.off()




