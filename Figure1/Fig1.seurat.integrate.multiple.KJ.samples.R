library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

load("../seurat/KJ12/KJ125.first.RData")
load("../seurat/KJ13/KJ135.first.withCREERT2.RData")
load("../seurat/KJ14/KJ145.first.RData")
load("../seurat/KJ15/KJ155.first.RData")
load("../seurat/KJ18/KJ185.first.RData")
load("../seurat/KJP5/KJP5.first.RData")

#KJ.merge<-merge(KJ125, list(KJ135, KJ145, KJ155, KJ185, KJP5))

KJ.list<-list(KJ125, KJ135, KJ145, KJ155, KJ185, KJP5)
KJ.list <- lapply(X = KJ.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = KJ.list)
KJ.anchors <- FindIntegrationAnchors(object.list = KJ.list, anchor.features = features)
save(KJ.anchors ,file="mKJ.anchors.RData")
KJ.combined <- IntegrateData(anchorset = KJ.anchors)
save(KJ.combined ,file="mKJ.combined.RData")

KJ.combined<- ScaleData(KJ.combined, verbose = FALSE)
KJ.combined<- RunPCA(KJ.combined, npcs = 30, verbose = FALSE)
KJ.combined<- RunUMAP(KJ.combined, reduction = "pca", dims = 1:30)
KJ.combined<<-RunTSNE(KJ.combined, reduction= "pca", dims.use = 1:30, do.fast = T)

save(KJ.combined ,file="mKJ.combined.UMAP.RData")

UMAPPlot(KJ.combined , pt.size = 2,label=T,label.size=12)

indOut<-abs(KJ.combined@reductions$umap@cell.embeddings[,2])>7.5 | KJ.combined@reductions$umap@cell.embeddings[,1]< -7.5
KJ.combined[["outlier"]]<-indOut

plot(KJ.combined@reductions$umap@cell.embeddings,pch=".",
	col=(indOut+1))

KJ.combCore<- subset(KJ.combined, subset = outlier == FALSE)
KJ.combCore<- ScaleData(KJ.combCore, verbose = FALSE)
KJ.combCore<- RunPCA(KJ.combCore, npcs = 30, verbose = FALSE)
KJ.combCore<- RunUMAP(KJ.combCore, reduction = "pca", dims = 1:30)
KJ.combCore<<-RunTSNE(KJ.combCore, reduction= "pca", dims.use = 1:30, do.fast = T)
KJ.combCore<-FindNeighbors(KJ.combCore, reduction="pca", dims = 1:30)
KJ.combCore<-FindClusters(KJ.combCore, resolution = 0.25)

ind135<-which(KJ.combCore@meta.data$orig.ident=="KJ135")
ind1352<-match(gsub("_.*","",colnames(KJ.combCore)[ind135]),colnames(KJ135))
KJ.combCore@meta.data$CREERT_nelson35[ind135]<-KJ135@meta.data$CREERT_2[ind1352]
ind125<-which(KJ.combCore@meta.data$orig.ident=="KJ125")
KJ.combCore@meta.data$CREERT_nelson35[ind125]<-KJ.combCore@meta.data$gpCAG.CreERT2[ind125]
#################################
table(KJ.combCore@meta.data$orig.ident, KJ.combCore@meta.data$CREERT_nelson35>0)


pdf("mKJ.combCore.pdf")
	DimPlot(KJ.combCore, reduction = "umap", group.by = "orig.ident")
	UMAPPlot(KJ.combCore, pt.size = 2,label=T,label.size=12)
	TSNEPlot(KJ.combCore, pt.size = 2,label=T,label.size=12)
	PCAPlot(KJ.combCore, pt.size = 2,label=T,label.size=12)
dev.off()

pdf("mKJ.combCore.split.pdf",width=36)
	UMAPPlot(KJ.combCore, pt.size = 2,split.by="orig.ident",
		label=T,label.size=12)
dev.off()

#################################
COLN2<-gsub("_.*","",colnames(KJ.combCore))

KJ.combCore[["COL2"]]<-0
KJ.combCore[["COL1"]]<-0
KJ.combCore[["SOX9"]]<-0
KJ.combCore[["GDF5"]]<-0
KJ.combCore[["LGR5"]]<-0
KJ.combCore[["MKI67"]]<-0

for(i in 1:6){
	OBJ<-KJ.list[[i]]
	TP<-c("KJ125", "KJ135", "KJ145", "KJ155", "KJ185", "KJP5")[i]
	indTP<-which(KJ.combCore@meta.data$orig.ident==TP)
	indObj<-match(COLN2[indTP],colnames(OBJ))
	
	KJ.combCore[["COL2"]][indTP,1]<-log2(GetAssayData(OBJ, slot="counts")["Col2a1",indObj]+1)
	KJ.combCore[["COL1"]][indTP,1]<-log2(GetAssayData(OBJ, slot="counts")["Col1a1",indObj]+1)
	KJ.combCore[["SOX9"]][indTP,1]<-log2(GetAssayData(OBJ, slot="counts")["Sox9",indObj]+1)
	KJ.combCore[["GDF5"]][indTP,1]<-log2(GetAssayData(OBJ, slot="counts")["Gdf5",indObj]+1)
	KJ.combCore[["LGR5"]][indTP,1]<-log2(GetAssayData(OBJ, slot="counts")["Lgr5",indObj]+1)
	KJ.combCore[["MKI67"]][indTP,1]<-log2(GetAssayData(OBJ, slot="counts")["Mki67",indObj]+1)
	KJ.combCore[["LGR5"]][indTP,1]<-log2(GetAssayData(OBJ, slot="counts")["Lgr5",indObj]+1)


}
#################################
library(scales)

KJ.combCore[["GDF5pos"]]<- KJ.combCore[["GDF5"]][,1]>0
KJ.combCore[["LGRpos"]]<- KJ.combCore[["LGR5"]][,1]>0
KJ.combCore[["CREERT_pos"]]<- KJ.combCore[["CREERT_nelson35"]][,1]>0
KJ.combCore[["LGRpos_CREERT_pos"]]<- KJ.combCore[["CREERT_nelson35"]][,1]>0 |  KJ.combCore[["LGR5"]][,1]>0

tab13<-table(KJ.combCore@meta.data$seurat_clusters,KJ.combCore@meta.data$orig.ident)
ftab13<-apply(tab13,2,function(x)x/sum(x)*100)
ftab14<-apply(tab13[-c(1,5),],2,function(x)x/sum(x)*100)

barplot(ftab13,col=hue_pal()(6))
barplot(ftab14,col=hue_pal()(6)[-c(1,5)])

tabGDF5<-table(paste0(KJ.combCore@meta.data$seurat_clusters,"__",KJ.combCore@meta.data$orig.ident),
	KJ.combCore@meta.data$GDF5pos)
tabGDF52<-as.data.frame(tabGDF5)%>%group_by(Var1)%>%  mutate(Num=sum(Freq))%>%
	mutate(prop = Freq/ sum(Freq) *100) %>%
	mutate(ypos = (cumsum(rev(prop))- 0.5*rev(prop)))%>%
	mutate(ypos2 = prop*0.75+ c(0, cumsum(prop)[-length(prop)]))%>%
	mutate(LAB=rev(paste0(percent(prop/100,0.1),"\n(n=",Freq,")")))%>%
	group_by(Var2)

tabLGRGFP<-table(paste0(KJ.combCore@meta.data$seurat_clusters,"__",KJ.combCore@meta.data$orig.ident),
	KJ.combCore@meta.data$LGRpos_CREERT_pos)
tabLGRGFP2<-as.data.frame(tabLGRGFP)%>%group_by(Var1)%>%  mutate(Num=sum(Freq))%>%
	mutate(prop = Freq/ sum(Freq) *100) %>%
	mutate(ypos = (cumsum(rev(prop))- 0.5*rev(prop)))%>%
	mutate(ypos2 = prop*0.75+ c(0, cumsum(prop)[-length(prop)]))%>%
	mutate(LAB=rev(paste0(percent(prop/100,0.1),"\n(n=",Freq,")")))%>%
	group_by(Var2)

pdf("GDF5.pos.mKJ.pies.pdf")
	ggplot(data.frame(tabGDF52), aes(x="", y=prop, fill=Var2)) +
  		geom_bar(stat="identity", width=1, color="white") +
	  	coord_polar("y", start=0, direction=1) + facet_wrap(vars(Var1)) +
    		theme(axis.text.x=element_blank(),
			panel.grid=element_blank(),axis.ticks = element_blank()) +
		geom_text(aes(y=ypos,label = LAB), size=3, color="white") +
		ggtitle("GDF5 pos in all mKJ\n1:positive\n2:negative")+
		xlab("")+ylab("")

	ggplot(data.frame(tabLGRGFP2), aes(x="", y=prop, fill=Var2)) +
  		geom_bar(stat="identity", width=1, color="white") +
	  	coord_polar("y", start=0, direction=1) + facet_wrap(vars(Var1)) +
    		theme(axis.text.x=element_blank(),
			panel.grid=element_blank(),axis.ticks = element_blank()) +
		geom_text(aes(y=ypos,label = LAB), size=3, color="white") +
		ggtitle("LGR/GFP pos in all mKJ\n1:positive\n2:negative")+
		xlab("")+ylab("")
dev.off()


KJ.combCore[["GFPpos"]]<- KJ.combCore@meta.data$CREERT_nelson35>0 |  KJ.combCore@meta.data$gpCAG.CreERT2>0
KJ.combCore[["GFPpos"]][is.na(KJ.combCore[["GFPpos"]])]<-F
table(KJ.combCore@meta.data$orig.ident,KJ.combCore@meta.data$GFPpos)

save(KJ.combCore, file="KJ.combCore.GFP.RData")
CoreMet<-KJ.combCore@meta.data
save(CoreMet,file="KJ.combCore.meta.RData")

	FeaturePlot(KJ.combCore, features=c("SOX9"),slot = "counts",
		pt.size = 2,split.by="orig.ident", col=c("grey","green"),
		label=T,label.size=12)

pdf("mKJ.combCore.features.split.pdf",width=36)
	FeaturePlot(KJ.combCore, features=c("GDF5"),slot = "counts",
		pt.size = 2,split.by="orig.ident", col=c("grey","green"),
		label=T,label.size=12)
	FeaturePlot(KJ.combCore, features=c("GDF5pos"),slot = "counts",
		pt.size = 2,split.by="orig.ident", col=c("grey","green"),
		label=T,label.size=12)
	FeaturePlot(KJ.combCore, features=c("GDF5"),slot = "counts",
		pt.size = 2,split.by="orig.ident", col=c("grey","#28B34B"),
		label=T,label.size=12)
	FeaturePlot(KJ.combCore, features=c("GDF5pos"),slot = "counts",
		pt.size = 2,split.by="orig.ident", col=c("grey","#28B34B"),
		label=T,label.size=12)

	FeaturePlot(KJ.combCore, features=c("LGRpos"),slot = "counts",
		pt.size = 2,split.by="orig.ident", col=c("grey","#28B34B"),
		label=T,label.size=12)
	FeaturePlot(KJ.combCore, features=c("GFPpos"),slot = "counts",
		pt.size = 2,split.by="orig.ident", col=c("grey","#28B34B"),
		label=T,label.size=12)

	FeaturePlot(KJ.combCore, features=c("CREERT_pos"),slot = "counts",
		pt.size = 2,split.by="orig.ident", col=c("grey","#28B34B"),
		label=T,label.size=12)

	FeaturePlot(KJ.combCore, features=c("LGRpos_CREERT_pos"),slot = "counts",
		pt.size = 2,split.by="orig.ident", col=c("grey","#28B34B"),
		label=T,label.size=12)
dev.off()
	FeaturePlot(subset(KJ.combCore,ident=1), features=c("GDF5"),slot = "counts",col=c("grey","green"),
		pt.size = 2,label=T,label.size=12)

	FeaturePlot(subset(KJ.combCore,ident=1), features=c("GDF5"),slot = "counts",col=c("grey","green"),
		pt.size = 2,label=T,label.size=12)
#save(KJ.combCore,file="mKJ.combCore.RData")

FeaturePlot(KJ.combCore,label=T, features=c("Mki67"),
		cols=c("lightgrey","brown"),pt.size=1/2)
VlnPlot(KJ.combCore,features="Mki67",pt.size=0)
tab12<-table(paste0(KJ.combCore@meta.data$seurat_clusters, ":", KJ.combCore@meta.data$orig.ident), GetAssayData(KJ.combCore)["Mki67",]>5)
par(mar=c(8,4,4,2));barplot(tab12[,2]/rowSums(tab12), las=2)

FeaturePlot(KJ.combCore,label=T, features=c("Epyc","Ucma","Sox9",
		"Col2a1","Mki67","Col1a1","Prrx1"),slot = "counts",
		split.by="orig.ident",reduction="tsne",
		cols=c("lightgrey","brown"),pt.size=1/2)

markers.combCore<- FindAllMarkers(KJ.combCore, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
#save(markers.combCore, file="mKJ.combCore.markers.RData")

LSIG_combCore<-sapply(levels(markers.combCore$cluster),function(i){
	thisDEG<-markers.combCore[markers.combCore$cluster==i&markers.combCore$p_val_adj<0.05,]
	thisDEG<-thisDEG[(thisDEG[,3]-thisDEG[,4])>0.25 | thisDEG[,2]>0.5 ,] #
	thisDEG<-thisDEG[order(thisDEG[,3]-thisDEG[,4],decreasing=T),]
	as.character(thisDEG$gene)
})

VlnPlot(KJ.combCore,features=LSIG_combCore[[6]],pt.size=0)


pdf("heatmap.6.clusters.pdf",height=12)
	DoHeatmap(KJ.combCore,features=unique(unlist(LSIG_combCore)))
dev.off()

KJ.combCore2<- subset(KJ.combined, subset = outlier == FALSE)

VlnPlot(KJ.merge,features=c("Gapdh","Actb"),group.by="orig.ident",pt.size=0)
VlnPlot(KJ.combCore,features=c("Gapdh","Actb"),group.by="orig.ident",pt.size=0)
VlnPlot(KJ.combCore2,features=c("Gapdh","Actb"),group.by="orig.ident",pt.size=0)

VlnPlot(KJ.combCore2,features=unique(unlist(LSIG_combCore))[17:32],group.by="orig.ident",pt.size=0)

KJ.combCore[["Sample_Clust"]]<-paste0(KJ.combCore@meta.data$orig.ident,"__C",Idents(KJ.combCore))
KJ.combCore[["Clust_Sample"]]<-paste0("C",Idents(KJ.combCore),"__",KJ.combCore@meta.data$orig.ident)

KJ.combCore3<- ScaleData(KJ.combCore,features=unique(unlist(LSIG_combCore)), verbose = FALSE)


DotPlot(KJ.combCore3,features=unique(unlist(LSIG_combCore)),group.by="Sample_Clust")
DotPlot(KJ.combCore3,features=unique(unlist(LSIG_combCore)),group.by="Clust_Sample")

d1<-DotPlot(KJ.combined,features=unique(unlist(LSIG_combCore)),group.by="orig.ident")

barplot(rep(1,41),col=KOL40)

names(KOL40)
d1<-DotPlot(KJ.combCore3,features=unique(unlist(LSIG_combCore)),group.by="Clust_Sample")
KOL40<-colorRampPalette(c("#619CFF", "grey","#F8766D"))(length(unique(round(d1$data$avg.exp.scaled*10))))
names(KOL40)<-sort(unique(round(d1$data$avg.exp.scaled*10)))

blank_theme <- theme_minimal()+
  theme(
	axis.text = element_blank(),
  axis.title= element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  panel.spacing = unit(-0.25, "lines"),
  axis.ticks = element_blank(),
  plot.title=element_blank(),
	 strip.text=element_blank()
  )

pdf("array.of.piecharts.pdf",width=24,height=12)
	ggplot(d1$data[,]%>%mutate(GRP=paste0(id,"_",features.plot),
			KOLOR=as.character(round(avg.exp.scaled*10)),
			GRP2=row_number()), 
			aes(x="", y=pct.exp, fill=KOLOR)) +
		geom_bar(stat="identity", width=1) +
		scale_fill_manual(name = "KOLOR",values = KOL40) +
		coord_polar("y", start=0) + facet_wrap(vars(GRP2),nrow=30) +blank_theme
dev.off()


DotPlot(KJ125,features=unique(unlist(LSIG_combCore)),group.by="orig.ident")






