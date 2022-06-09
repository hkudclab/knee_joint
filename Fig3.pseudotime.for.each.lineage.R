library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(igraph)
library(scales)

library(monocle3)

load("../seurat/monocle3.data.format.RData")

load("Individual.samples.with.GFP.12.13.noCycling.RData")
load("Individual.samples.with.GFP.13.14.noCycling.RData")
load("Individual.samples.with.GFP.14.15.noCycling.RData")
load("Individual.samples.with.GFP.15.18.noCycling.RData")
load("Individual.samples.with.GFP.18.p5.noCycling.RData")
load("Individual.samples.with.GFP.18.MCp5.noCycling.RData")

load("../resources/Ligands.mouse.human.gene-symb.RData")
load("../resources/matrisome-mit-mouse/Matrisome.mouse.RData")

load("step8_KJP5-MCP5-merge/SZKJP5.RData")
SZKJP5[["BioName"]]<-"SZ"

KJ125.3[["BioName"]]<-c("0"="JPC/Chon-prog","1"="Lig","2"="Chon-prog","3"="Cyc","4"="Cyc","5"="JPC-Lig","6"="Lig-prog","7"="GPC-prog","8"="Mes-EMT")[KJ125.3@meta.data$seurat_clusters]
 table(KJ125.3[["BioName"]][,1], KJ125.3@meta.data$seurat_clusters)
KJ135.3[["BioName"]]<-c("0"="GPC-chon","1"="JPC","2"="JPC/AC-prog1","3"="AC-prog","4"="Lig-prog","5"="Cyc","6"="JPC/AC-prog2","7"="Mes-Lig")[KJ135.3@meta.data$seurat_clusters]
 table(KJ135.3[["BioName"]][,1], KJ135.3@meta.data$seurat_clusters)
KJ145.3[["BioName"]]<-c("0"="GPC-chon","1"="AC-prog","2"="AC/SZ-prog","3"="Lig/CL-prog","4"="JPC","5"="Lig-prog","6"="Lig","7"="Cyc-Mes")[KJ145.3@meta.data$seurat_clusters]
 table(KJ145.3[["BioName"]][,1], KJ145.3@meta.data$seurat_clusters)
KJ155.3[["BioName"]]<-c("0"="GPC-chon","1"="AC-prog","2"="CL-prog","3"="GPC","4"="Lig","5"="MClig-prog","6"="Chon","7"="Cyc","8"="SZ","9"="Lig")[KJ155.3@meta.data$seurat_clusters]
 table(KJ155.3[["BioName"]][,1], KJ155.3@meta.data$seurat_clusters)
KJ185.3[["BioName"]]<-c("0"="GPC-chon","1"="CL","2"="GPC","3"="Chon","4"="MCchon","5"="AC","6"="Lig","7"="MClig","8"="SZ","9"="Cyc", "10"="Chon")[KJ185.3@meta.data$seurat_clusters]
 table(KJ185.3[["BioName"]][,1], KJ185.3@meta.data$seurat_clusters)
KJP5.3[["BioName"]]<-c("0"="GPC-chon","1"="enthesis","2"="AC","3"="CL","4"="MCchon","5"="AC/MC-chon","6"="GPC","7"="GPC","8"="MClig/SZ")[KJP5.3@meta.data$seurat_clusters]
 table(KJP5.3[["BioName"]][,1], KJP5.3@meta.data$seurat_clusters)
MCP5.3[["BioName"]]<-c("0"="MC-chon","1"="mMC","2"="MC/AC","3"="MC-lig","4"="yMC","5"="SZ","6"="MC-chon")[MCP5.3@meta.data$seurat_clusters]
 table(MCP5.3[["BioName"]][,1], MCP5.3@meta.data$seurat_clusters)

UMAPPlot(KJ125.3, label=T, label.size=6, pt.size=4) + ggtitle("KJ12.5")


indGPC<-list(7, 0, 0, 0, 0, c(0,6), NA)
indSZ<-list(0, c(1, 2,6), 2, 8, 8, 8, 5)
indSZ<-list(0, c(2,6), 2, 8, 8, 8, 5)

indAC<-list(0, c(2,3,6), 1, 1, 5, 2, 2)
indMCchon<-list(0, c(2,3,6), 1, 1, 4, 4, c(0,6))
indCL<-list(c(0,6,8), c(1,4), c(3,4,5), 2, 1, 3, c(1,4))
indMClig<-list(c(6,8), 4, c(5), c(5), 7, 8, 3)

METALL<-data.frame("orig.ident"="", "nCount_RNA"=0, "nFeature_RNA"=0, "seurat_clusters"=0, "CREERT_nelson35"=0, "BioName"="", "combID"="")
METALL<-METALL[NULL, ]
for(OBJ in list(KJ125.3, KJ135.3, KJ145.3, KJ155.3, KJ185.3, KJP5.3, MCP5.3)){
	cat("########################################\n")
	#tab1<-table(OBJ[["BioName"]][,1], GetAssayData(OBJ)["Lgr5",]>0 | OBJ@meta.data$CREERT_nelson35>0)
	#tab2<-table(OBJ[["BioName"]][,1], GetAssayData(OBJ)["Gdf5",]>0 )
	tab1<-table(paste0(OBJ[["BioName"]][,1],"_", OBJ@meta.data$seurat_clusters), GetAssayData(OBJ)["Lgr5",]>0 | OBJ@meta.data$CREERT_nelson35>0)
	tab2<-table(paste0(OBJ[["BioName"]][,1],"_", OBJ@meta.data$seurat_clusters), GetAssayData(OBJ)["Gdf5",]>0 )
	mat1<-cbind(tab1, percent(tab1[,2]/rowSums(tab1),0.1), tab2, percent(tab2[,2]/rowSums(tab2),0.1))
	TMP<-data.frame(OBJ@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA", "seurat_clusters", "CREERT_nelson35", "BioName")], combID=colnames(OBJ))
	METALL<-rbind(METALL, TMP)
	print(mat1)
}
METALL[["combID2"]]<-paste0(METALL$combID,"_",as.integer(METALL$orig.ident))
table(METALL$orig.ident=="KJP5" & METALL$combID%in%colnames(SZKJP5))
METALL[METALL$orig.ident=="KJP5" & METALL$combID%in%colnames(SZKJP5),"BioName"]<-"SZ"

branches<-list("GPC"=indGPC, "SZ"=indSZ, "AC"=indAC, "MCchon"=indMCchon, "CL"=indCL, "MClig"=indMClig)
thisbranch<-"MClig"
listpop<-branches[[thisbranch]]

CDS.all<-list(cds.125[,colnames(KJ125.3)], cds.135[,colnames(KJ135.3)], cds.145[,colnames(KJ145.3)], cds.155[,colnames(KJ155.3)], 
	cds.185[,colnames(KJ185.3)], cds.p5[,colnames(KJP5.3)], cds.mcp5[,colnames(MCP5.3)])
LISTall<-list(KJ125.3, KJ135.3, KJ145.3, KJ155.3, KJ185.3, KJP5.3, MCP5.3)


SZ.list<-list()
CDS.list<-list()

for(i in 1:7){
	INDi<-listpop[[i]]
	if(!is.na(INDi[1])){
		SZ.list[[i]]<-subset(LISTall[[i]], ident=INDi)
		if(thisbranch=="SZ"&&i==6)
			SZ.list[[i]]<-SZKJP5
		CDS.list[[i]]<-CDS.all[[i]][,colnames(SZ.list[[i]])]
	}
}

#####################################
#####################################
#METASZ<-SZ.combined@meta.data
METASZ<-droplevels(METALL[as.integer(METALL$orig.ident)%in% which(sapply(CDS.list, length)>0), ])

CDS.list<-CDS.list[sapply(CDS.list, length)>0]
big_cds <- combine_cds(CDS.list)

METASZ[["combID2"]]<-paste0(METASZ$combID,"_",as.integer(METASZ$orig.ident))

table(colnames(big_cds) %in% METASZ[["combID2"]])
METASZ<-METASZ[match(colnames(big_cds), METASZ[["combID2"]]),]

#big_cds@colData[["orig.ident"]]<-gsub(".*-","",colnames(big_cds))
#big_cds<- align_cds(big_cds, alignment_group = "orig.ident")
#big_cds2<-big_cds
#big_cds<-big_cds[,rownames(METASZ)]
#table(colnames(SZ.combined)==colnames(big_cds))

big_cds@colData[["orig.ident"]]<-METASZ[["orig.ident"]]
big_cds@colData[["nCount_RNA"]]<-METASZ[["nCount_RNA"]]
big_cds@colData[["nFeature_RNA"]]<-METASZ[["nFeature_RNA"]]
big_cds@colData[["seurat_clusters"]]<-METASZ[["seurat_clusters"]]
big_cds@colData[["CREERT_nelson35"]]<-METASZ[["CREERT_nelson35"]]
big_cds@colData[["BioName"]]<-METASZ[["BioName"]]
big_cds@colData[["BioName2"]]<-paste0(METASZ[["orig.ident"]],"_",METASZ[["BioName"]])
indnonSZ<-which(big_cds@colData[["BioName2"]]!="KJP5_SZ")


big_cds<-big_cds[,indnonSZ]
big_cds <- preprocess_cds(big_cds, num_dim = 20)
big_cds<- align_cds(big_cds, alignment_group = "orig.ident")

plot_pc_variance_explained(big_cds)
big_cds<- reduce_dimension(big_cds)

big_cds<- cluster_cells(big_cds)#, reduction_method="UMAP", resolution=2)
#big_cds<- cluster_cells(big_cds, reduction_method="PCA", resolution=2)
#big_cds<- cluster_cells(big_cds, reduction_method="Aligned")

RSQR<-matrix(NA,20,2)
for(i in 1:20){
	resi<-summary(lm(big_cds@int_colData@listData$reducedDims$PCA[,i]~big_cds@colData[["orig.ident"]]))
	RSQR[i,]<-c(resi$r.squared,resi$adj.r.squared)
}

big_cds@int_colData@listData$reducedDims@listData$UMAP <- big_cds@int_colData@listData$reducedDims$PCA[,c(1,4)]

big_cds<- learn_graph(big_cds, use_partition = F, close_loop = F) #, learn_graph_control=list(ncenter=13))

plot_cells(big_cds, color_cells_by="BioName2", reduction_method="UMAP", label_groups_by_cluster=FALSE,  cell_size =2, group_label_size = 6,  graph_label_size=5)


###########
###########

plot_cells(big_cds, reduction_method="Aligned", label_groups_by_cluster=FALSE,  cell_size =2, group_label_size = 6,  graph_label_size=5)

plot_cells(big_cds, color_cells_by="BioName2", reduction_method="PCA", label_groups_by_cluster=FALSE,  cell_size =2, group_label_size = 6,  graph_label_size=5)

pdf("step7-shockwaves/branch.MCLig.pseudotime.pdf", width=8)
	plot_cells(big_cds, color_cells_by="BioName2", 
		reduction_method="UMAP", label_groups_by_cluster=FALSE,  cell_size =2, group_label_size = 6,  graph_label_size=5)
dev.off()

###########
big_cds<- learn_graph(big_cds, use_partition = F, close_loop = F)
big_cds<- order_cells(big_cds)

bymedian <- reorder(big_cds@colData$orig.ident, pseudotime(big_cds), median)
bymedian2 <- reorder(big_cds@colData$BioName2, pseudotime(big_cds), median)

par(mar=c(12,4,4,2))
boxplot(pseudotime(big_cds)~bymedian2, las=2, xlab="")

plot(pseudotime(big_cds), GetAssayData(SZ.combined, assay = "RNA",slot="counts")["Lgr5",])
plot(pseudotime(big_cds), GetAssayData(SZ.combined, assay = "RNA",slot="counts")["Gdf5",])



#####################################

#SZ.list<-list(jpcchon12, jpc13, chonsz14, sz15, sz18, szp5)
SZ.list<-SZ.list[sapply(SZ.list, length)>0]
SZ.list <- lapply(X = SZ.list, FUN = function(x) {
	x <- NormalizeData(x)
	x <- FindVariableFeatures(x, selection.method = "vst")
})

features <- SelectIntegrationFeatures(object.list = SZ.list)
SZ.anchors <- FindIntegrationAnchors(object.list = SZ.list, anchor.features = features)
SZ.combined <- IntegrateData(anchorset = SZ.anchors,  k.weight = 1/10)

SZ.combined[["BioName2"]]<-paste0(SZ.combined@meta.data$orig.ident,"_",SZ.combined@meta.data$BioName)

SZ.combined<- ScaleData(SZ.combined, verbose = FALSE)
SZ.combined<- RunPCA(SZ.combined, npcs = 30, verbose = FALSE)
SZ.combined<- RunUMAP(SZ.combined, reduction = "pca", dims = 1:30)
SZ.combined<-RunTSNE(SZ.combined, reduction= "pca", dims.use = 1:30, do.fast = T)
SZ.combined<-FindNeighbors(SZ.combined, reduction="pca", dims = 1:30)
SZ.combined<-FindClusters(SZ.combined, resolution = 0.25)

library(entropy)
 table(colnames(big_cds)==colnames(SZ.combined))

ENTRO<-c()
for(i in 1:dim(SZ.combined)[2]){
	if(i%%200==0)cat(i,"\t",date(),"\n")
	y<-as.matrix(GetAssayData(SZ.combined)[,i])[,1]
	z<-log2(y[y>0])
	ENTRO[i]<-entropy(discretize(z,20))
}

plot(ENTRO[indnonSZ], pseudotime(big_cds))

bymedian3 <- reorder(big_cds@colData$BioName2, ENTRO[indnonSZ], median)
boxplot(ENTRO[indnonSZ]~bymedian3 , las=2, xlab="")
#####################################
library(ggridges)

bymedian2 <- reorder(big_cds@colData$BioName2, pseudotime(big_cds), mean)

DAT.pseudo<-data.frame(bymedian,bymedian2 , BioName2=big_cds@colData$BioName2, pseudotime=pseudotime(big_cds))
pdf("step7-shockwaves/ridge.pseudo.MCLIG.pop.pdf", width=10)
	ggplot(DAT.pseudo, aes(x = pseudotime, y = bymedian)) +
		geom_density_ridges(scale = 4) + 
		#stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7) +
		scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
		scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
		coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
		theme_ridges()

	ggplot(DAT.pseudo, aes(x = pseudotime, y = bymedian2 )) +
		geom_density_ridges(scale = 4) + 
		#stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7) +
		scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
		scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
		coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
		theme_ridges()
dev.off()
#####################################

save(big_cds, file=paste0("step7-shockwaves/MCLig.big_cds.RData"))

#####################################

UMAPPlot(SZ.combined , pt.size = 2,label=T,label.size=12) + NoLegend()

TSNEPlot(SZ.combined , pt.size = 2,label=T,label.size=12) + NoLegend()

UMAPPlot(SZ.combined , pt.size = 2,label=T,label.size=12, group.by="BioName2")

table(SZ.combined@meta.data$seurat_clusters, SZ.combined@meta.data$orig.ident)

table(SZ.combined@meta.data$seurat_clusters, SZ.combined@meta.data$orig.ident)
table(GetAssayData(SZ.combined, assay = "RNA",slot="counts")["Lgr5",]>0, SZ.combined@meta.data$orig.ident)
SZ.combined[["LGR5pos"]]<-GetAssayData(SZ.combined, assay = "RNA",slot="counts")["Lgr5",]>0
UMAPPlot(SZ.combined , pt.size = 4,label=T,label.size=12, group.by="LGR5pos")
UMAPPlot(SZ.combined , pt.size = 4,label=T,label.size=12)

FeaturePlot(SZ.combined, features="Lgr5",slot="counts", pt.size=4)
tabLgr5<-table(GetAssayData(SZ.combined, assay = "RNA",slot="counts")["Lgr5",]>0, SZ.combined@meta.data$seurat_clusters)

tabLgr5<-t(table(GetAssayData(SZ.combined, assay = "RNA",slot="counts")["Lgr5",]>0, paste0(SZ.combined@meta.data$orig.ident,"_",SZ.combined@meta.data$BioName)))
cbind(as.matrix(tabLgr5), tabLgr5[,2]/rowSums(tabLgr5))

barplot(apply(table(SZ.combined@meta.data$seurat_clusters, SZ.combined@meta.data$orig.ident),2,function(x)x/sum(x)), col=seq(7))
 

#####################################

GETSIG<-function(MARKERS){
	sapply(levels(MARKERS$cluster),function(i){
		thisDEG<-MARKERS[MARKERS$cluster==i&MARKERS$p_val_adj<0.05,]
		thisDEG<-thisDEG[(thisDEG[,3]-thisDEG[,4])>0.25 | thisDEG[,2]>0.5 ,] #
		thisDEG<-thisDEG[order(thisDEG[,3]-thisDEG[,4],decreasing=T),]
		as.character(thisDEG$gene)
	})
}
#markers.SZ <- FindAllMarkers(SZ.combined, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
LSIG_SZ<-GETSIG(markers.SZ)
##############################
##############################
#########################################################
NumWave<-5
NumPoint=1000

DNORM<-dnorm(seq(-6,6,len=(NumWave-1)*NumPoint*2+1))
plot(DNORM)


DROMmat<-matrix(NA, (NumWave-1)*NumPoint, NumWave)
for(i in 1:NumWave){
	Ni<- NumWave*NumPoint - i*NumPoint 
	DNORMi<-DNORM[-seq(Ni)][1:nrow(DROMmat)]
	DROMmat[,i]<-DNORMi
}
matplot(DROMmat, type="l", lty=1)

DATnorm<-data.frame(x=rep(seq(nrow(DROMmat)),NumWave),
	y=as.vector(DROMmat), 
	z=rep(seq(NumWave),rep(nrow(DROMmat),NumWave)))
ggplot(DATnorm, aes(x, y, group=z))+
	geom_point(aes(x,y)) +
	facet_wrap(vars(z), ncol=1)

#plot(DNORM)

PSEUDO<-round(pseudotime(big_cds)/ max(pseudotime(big_cds))*nrow(DROMmat))
PSEUDO[PSEUDO==0]<-1
range(PSEUDO)
hist(PSEUDO)

####################################################
####################################################
load("../resources/TF.CD.receptors.human.mouse.RData")
load("../resources/human.mouse.receptors.RData")
big_cds_500<-big_cds[which(rowSums(exprs(big_cds))>500),]
big_cds_TF<-big_cds[which(rowData(big_cds)$gene_short_name%in%TFm& rowSums(exprs(big_cds))>500),]
big_cds_Receptors<-big_cds[which(rowData(big_cds)$gene_short_name%in%mouse.receptors& rowSums(exprs(big_cds))>500),]
big_cds_Ligands<-big_cds[which(rowData(big_cds)$gene_short_name%in%mouse.Ligands & rowSums(exprs(big_cds))>500),]
big_cds_ECM<-big_cds[which(rowData(big_cds)$gene_short_name%in%unique(unlist(LgMatrisome)) & rowSums(exprs(big_cds))>500),]

LCOEF_TF<-list()
LCOEF_Ligands<-list()
LCOEF_Receptors<-list()
LCOEF_ECM<-list()

LCOEF_ALL<-list()

for(i in 1:NumWave){
	cat(i, ":\t", date(), "\n")
	big_cds_TF@colData$DNORM<-DROMmat[,i][PSEUDO+1]

	TF_fits <- fit_models(big_cds_TF, model_formula_str = "~DNORM")
	#TF_coefs <- coefficient_table(TF_fits)
	TF_coefs <- select(coefficient_table(TF_fits), -one_of("model", "model_summary")) %>% filter(term=="DNORM")

	LCOEF_TF[[i]]<-TF_coefs
}

for(i in 1:NumWave){
	cat(i, ":\t", date(), "\n")
	big_cds_Ligands@colData$DNORM<-DROMmat[,i][PSEUDO+1]

	Ligands_fits <- fit_models(big_cds_Ligands, model_formula_str = "~DNORM")
	#Ligands_coefs <- coefficient_table(Ligands_fits)
	Ligands_coefs <- select(coefficient_table(Ligands_fits), -one_of("model", "model_summary")) %>% filter(term=="DNORM")
	LCOEF_Ligands[[i]]<-Ligands_coefs
}
for(i in 1:NumWave){
	cat(i, ":\t", date(), "\n")
	big_cds_Receptors@colData$DNORM<-DROMmat[,i][PSEUDO+1]

	Receptors_fits <- fit_models(big_cds_Receptors, model_formula_str = "~DNORM")
	#Receptors_coefs <- coefficient_table(Receptors_fits)
	Receptors_coefs <- select(coefficient_table(Receptors_fits), -one_of("model", "model_summary")) %>% filter(term=="DNORM")
	LCOEF_Receptors[[i]]<-Receptors_coefs
}

for(i in 1:NumWave){
	cat(i, ":\t", date(), "\n")
	big_cds_ECM@colData$DNORM<-DROMmat[,i][PSEUDO+1]

	ECM_fits <- fit_models(big_cds_ECM, model_formula_str = "~DNORM")
	#ECM_coefs <- coefficient_table(ECM_fits)
	ECM_coefs <- select(coefficient_table(ECM_fits), -one_of("model", "model_summary")) %>% filter(term=="DNORM")
	LCOEF_ECM[[i]]<-ECM_coefs
}

for(i in 1:NumWave){
	cat(i, ":\t", date(), "\n")
	big_cds_500@colData$DNORM<-DROMmat[,i][PSEUDO+1]

	ALL_fits <- fit_models(big_cds_500, model_formula_str = "~DNORM")
	#ALL_coefs <- coefficient_table(ALL_fits)
	ALL_coefs <- select(coefficient_table(ALL_fits), -one_of("model", "model_summary")) %>% filter(term=="DNORM")
	LCOEF_ALL[[i]]<-ALL_coefs
}

save(LCOEF_TF,
	LCOEF_Ligands,
	LCOEF_Receptors,
	LCOEF_ECM,
	LCOEF_ALL,  file=paste0("step7-shockwaves/MCLig.LCOEF.",NumWave,".waves.RData"))

##################

LHEAD_TF<-list()
LHEAD_Ligands<-list()
LHEAD_Receptors<-list()
LHEAD_ECM<-list()
LHEAD_stage<-list()

LHEAD_ALL<-list()
for(i in 1:NumWave){
	TF_coefs<-LCOEF_TF[[i]]
	TF_coefs_top<-TF_coefs%>%filter(status=="OK" & q_value<0.01 & term=="DNORM"& (estimate)>1) %>% arrange(desc(estimate))
	LHEAD_TF[[i]]<-TF_coefs_top$gene_short_name

	Ligands_coefs<-LCOEF_Ligands[[i]]
	Ligands_coefs_top<-Ligands_coefs%>%filter(status=="OK"& q_value<0.01 & term=="DNORM" & (estimate)>1) %>% arrange(desc(estimate))
	LHEAD_Ligands[[i]]<-Ligands_coefs_top$gene_short_name

	Receptors_coefs<-LCOEF_Receptors[[i]]
	Receptors_coefs_top<-Receptors_coefs%>%filter(status=="OK"& q_value<0.01 & term=="DNORM" & (estimate)>1) %>% arrange(desc(estimate))
	LHEAD_Receptors[[i]]<-Receptors_coefs_top$gene_short_name

	ECM_coefs<-LCOEF_ECM[[i]]
	ECM_coefs_top<-ECM_coefs%>%filter(status=="OK"& q_value<0.01 & term=="DNORM" & (estimate)>1) %>% arrange(desc(estimate))
	LHEAD_ECM[[i]]<-ECM_coefs_top$gene_short_name
	
	LHEAD_stage[[i]]<-c(LHEAD_TF[[i]], LHEAD_Ligands[[i]], LHEAD_ECM[[i]])
}

for(i in 1:NumWave){
	ALL_coefs<-LCOEF_ALL[[i]]
	ALL_coefs_top<-ALL_coefs%>%filter(status=="OK"& q_value<0.01 & term=="DNORM" & (estimate)>1) %>% arrange(desc(estimate))
	LHEAD_ALL[[i]]<-setdiff(ALL_coefs_top$gene_short_name, "conf")
}

sapply(LHEAD_TF, paste, collapse=", ")
sapply(LHEAD_Ligands, paste, collapse=", ")
sapply(LHEAD_ECM, paste, collapse=", ")

str(LHEAD_TF)
str(LHEAD_Ligands)
str(LHEAD_ECM)

plot(GetAssayData(SZ.combined,assay = "RNA", slot="counts")["Col2a1",,drop=F], 
	GetAssayData(SZ.combined)["Col2a1",,drop=F], log="x")

###############################################
OUTDAT<-data.frame("Wave"=seq(NumWave), 
	"Ligands"=sapply(LHEAD_Ligands, function(x){paste0(paste0(x, collapse=","),"\n(n=",length(x),")")}),
	"Receptors"=sapply(LHEAD_Receptors, function(x){paste0(paste0(x, collapse=","),"\n(n=",length(x),")")}),
	"TFs"=sapply(LHEAD_TF, function(x){paste0(paste0(x, collapse=","),"\n(n=",length(x),")")}),
	"ECM"=sapply(LHEAD_ECM, function(x){paste0(paste0(x, collapse=","),"\n(n=",length(x),")")})
)
readr::write_csv(OUTDAT, paste0("step7-shockwaves/MCLig.",NumWave,".waves.wave-associated-genes.csv"))
###############################################
ttrust<-readr::read_tsv("P:/OneDrive - The University Of Hong Kong/resources/TF-targtet-TTRUST/trrust_rawdata.mouse.tsv")
TTRUST<-ttrust[match(unique(paste0(ttrust$TF, "_", ttrust$target)), paste0(ttrust$TF, "_", ttrust$target)),]



LNET<-list()
for(i in 1:NumWave){
	targets<-unique(unlist(LHEAD_stage[seq(i:min(NumWave, i+1))]))
	ttrusti<-TTRUST[TTRUST$target%in%targets,]
	ttrusti2<-ttrusti[ttrusti$TF%in%LHEAD_TF[[i]], ]
	intersect(ttrusti2$TF, LHEAD_TF[[i]])
	ttrusti2[["stage"]]<-i
	LNET[[i]]<-ttrusti2
}

g2 <- graph(edges=as.vector(t(TTRUST[,1:2])),  directed=T ) 

g1 <- graph( edges=as.vector(t(ttrusti[,1:2])), n=nrow(ttrusti), directed=T ) 
res1<-shortest_paths(g1, from = "Chd7", to = "Ezh2")
res1<-distances(g1, v = "Chd7", to = "Ezh2",   mode = "out")
res1<-distances(g1, v = "Ezh2", to = "Chd7",   mode = "out")
res1<-distances(g2, v = intersect(LHEAD_TF[[i]], TTRUST$target), to = intersect(targets, TTRUST$target),   mode = "out")
res1[res1==Inf]<-NA
sort( rowSums(res1<3,  na.rm=T))


 attr(res1$vpath[[1]], "names")


plot(g1, edge.arrow.size=.5, vertex.color="gold", vertex.size=4) 

sapply(names(LNET[[1]]), function(x)unlist(sapply(LNET, function(y)y[[x]])))
###############################################
pseudo.MCLig<-pseudotime(big_cds)
MCLig.combined<-subset(SZ.combined,, cells=colnames(SZ.combined,)[indnonSZ])

save(MCLig.combined, pseudo.MCLig, LHEAD_Ligands, LHEAD_Receptors, LHEAD_TF, LHEAD_ECM, 
	big_cds, file="step7-shockwaves/combined.and.pseudo.and.wave-sigs.MCLig.RData")
save(pseudo.MCLig, LIST1, LIST2,
	big_cds, file="step7-shockwaves/cds.and.pseudo.and.wave-sigs.MCLig.RData")

source("../helpers/step7.helper.TRENDPLOT.R")

pdf(paste0("step7-shockwaves/MCLig.",NumWave,".waves.trendline.loess.pdf"), width=12)
	LIST1<-list("Ligands"=LHEAD_Ligands, "Receptors"=LHEAD_Receptors, "TF"=LHEAD_TF, "ECM"=LHEAD_ECM)
	LIST2<-list("Ligands"=LCOEF_Ligands, "Receptors"=LCOEF_Receptors, "TF"=LCOEF_TF, "ECM"=LCOEF_ECM)
	TRENDPLOT2(big_cds,LIST1,LIST2,NumWave)
dev.off()



