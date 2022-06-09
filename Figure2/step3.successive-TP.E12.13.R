library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gplots)
library(ggalluvial)

load("../DATA/Individual.samples.with.GFP.12.13.RData")
#########################################
load("KJ.combCore.meta.RData")

#KJ1252<-subset(KJ125,cells=gsub("_.*","",rownames(CoreMet)[CoreMet$orig.ident=="KJ125"]))
#KJ1352<-subset(KJ135,cells=gsub("_.*","",rownames(CoreMet)[CoreMet$orig.ident=="KJ135"]))
KJ135.4<-RenameCells(KJ135.3,new.names=gsub("-1","-2",colnames(KJ135.3)) )
set.seed(0)
KJ1213<-RunCCA(KJ125.3,KJ135.4)


KJ1213[["OldClust"]]<-paste0(KJ1213@meta.data$orig.ident,"_",KJ1213@meta.data$seurat_clusters)
table(KJ1213[["OldClust"]][,1], KJ1213@meta.data$orig.ident)


KJ1213<-NormalizeData(KJ1213)
KJ1213<-ScaleData(KJ1213,fearures=rownames(KJ1213))
KJ1213<-FindVariableFeatures(KJ1213)
KJ1213<-RunPCA(KJ1213)
KJ1213<-RunUMAP(KJ1213, reduction= "cca", dims= 1:10)
KJ1213<-RunTSNE(KJ1213, reduction= "cca", dims.use = 1:40, do.fast = T)
KJ1213<-FindNeighbors(KJ1213, reduction= "cca", dims = 1:10)
KJ1213<-FindClusters(KJ1213, resolution = 0.25)
KJ1213.5<-FindClusters(KJ1213, resolution = 0.5)

table(KJ1213[["OldClust"]][,1], Idents(KJ1213))
table(KJ1213.5[["OldClust"]][,1], Idents(KJ1213.5))


save(KJ1213.5, file="successiveTP/successiveTP.12.13.RData")

#set.seed(0)
#KJ1213s<-subset(KJ1213,cells=sample(colnames(KJ1213),5e3))
#pdf("tsne.mKJ.12.13.pdf")
#	TSNEPlot(KJ1213,label=T,label.size=12, pt.size = 5)+ NoLegend()
#	TSNEPlot(KJ1213s,label=T,label.size=12, pt.size = 5, split.by="orig.ident")+ NoLegend()
#
#dev.off()

TSNEPlot(KJ1213.5, split.by="orig.ident", group.by="OldClust",
	label=T, label.size=6, pt.size=4)

#KJ1213@reductions$tsne@cell.embeddings


ind125<-which(KJ1213@meta.data$orig.ident=="KJ125")

MATdist1213<-as.matrix(dist(KJ1213.5@reductions$tsne@cell.embeddings[1:5,]))
ALLCLUST<-as.character(sort(unique(KJ1213.5@meta.data$seurat_clusters)))

MATflow<-matrix(NA,0,5)
colnames(MATflow)<-c("i", "cellIDfrom", "cellIDto", "OldClustThisTP", "FLOWto")

LIDto<-list()
LIDto2<-list()

for(i in 1:length(ALLCLUST)){
	indi<-which(KJ1213.5@meta.data$seurat_clusters==ALLCLUST[i])
	MATdist1213.i<-as.matrix(dist(KJ1213.5@reductions$tsne@cell.embeddings[indi,]))
	indi2<-which(KJ1213.5@meta.data$orig.ident[indi]=="KJ125")
	indi2.next<-which(KJ1213.5@meta.data$orig.ident[indi]!="KJ125")
	OldClustThisTP<-KJ1213.5[["OldClust"]][indi, 1][indi2]
	OldClustNextTP<-KJ1213.5[["OldClust"]][indi, 1][indi2.next]

	if(length(indi2)>0 & length(indi2.next)>0){
		MATdist1213.i2<-MATdist1213.i[indi2, indi2.next]
		CELLidsTO<-colnames(KJ1213.5)[indi][indi2.next]
		FLOWto<-apply(MATdist1213.i2, 1, function(x){
			ORD3<-order(x)[1:3]
			names(which.max(table(OldClustNextTP[ORD3])))
		})
		LIDto[[i]]<-FLOWtoCell
		ALLOCATED<-rep(0, ncol(MATdist1213.i2))
		FLOWtoCell<-sapply(1:nrow(MATdist1213.i2), function(j){
			x<-MATdist1213.i2[j, ]
			indUnallocated<-which(ALLOCATED==0)
			if(length(indUnallocated)==0)indUnallocated<-seq_along(ALLOCATED)
			ORD3<-order(x[indUnallocated])[1:3]
			flowToCluster<-names(which.max(table(OldClustNextTP[indUnallocated][ORD3])))

			ord4<-ORD3[which(OldClustNextTP[indUnallocated][ORD3]==flowToCluster)[1]]
			ALLOCATED[indUnallocated[ord4]]<<-1

			CELLidsTO[indUnallocated[ord4]]
		})
		LIDto2[[i]]<-FLOWtoCell
		thisTemp<-cbind(i=rep(i, length(OldClustThisTP)), 
			cellIDfrom=colnames(KJ1213.5)[indi][indi2],
			cellIDto=FLOWtoCell,
			OldClustThisTP, FLOWto)
		MATflow<-rbind(MATflow, thisTemp)
	}
}


save(MATflow, file="successiveTP/MATflow.E12.to.E13.RData")
###############################
library(RColorBrewer)



pdf("successiveTP/alluvial.chart.1213.pdf")
	ggplot(as.data.frame(table(MATflow[,4], MATflow[,5])),
      	 aes(y = Freq, axis1 = Var1, axis2 = Var2)) +
	  geom_alluvium(aes(fill = Var1), width = 1/12) +
	  geom_stratum(width = 1/12, fill = "grey", color = "white") +
	  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
	  scale_fill_discrete( colorRampPalette(brewer.pal(9, "Set1"))(length(table(MATflow[,2])))) 
dev.off()




