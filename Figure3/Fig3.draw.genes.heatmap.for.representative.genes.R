library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(igraph)
library(scales)
library(ggridges)

library(monocle3)

set.seed(0)
HUE20<-sample(hue_pal()(20))

display.brewer.all()

Wilson<-readxl::read_excel("step7-shockwaves/20220513 Wilson selected 5 genes based on trendlines.xlsx", sheet = "Sheet1")


FnData<-list("AC"="../DATA/cds.and.pseudo.and.wave-sigs.AC.RData", 
		"SZ"="../DATA/cds.and.pseudo.and.wave-sigs.SZ.RData", 
		"Liga"="../DATA/cds.and.pseudo.and.wave-sigs.Lig.RData", 
		"MCLig"="../DATA/cds.and.pseudo.and.wave-sigs.MCLig.RData")
struct<-c("AC"="AC", "SZ"="SZ", "Liga"="Lig", "MCLig"="MCLig")




pdfout<-paste0("step7-shockwaves/STEP-7e2.heatmap.wilson.selected.full-lists.per.gene.3-waves.4structs.pdf")
pdf(pdfout, height=12)
	for(STRUCT in c("AC", "SZ", "Liga", "MCLig")){
		load(FnData[[STRUCT]])
		Wilson2<-Wilson[,grep(paste0("^",STRUCT, "|wave"),colnames(Wilson))]

		Wilson3<-Wilson2 %>% 
  			pivot_longer(!waves, names_to = "GeneCat", values_to = "GeneSymb") %>% filter(!is.na(GeneSymb))

		OBJ<-big_cds

		thisNumWave<-c("AC"=4, "SZ"=3, "Liga"=5, "MCLig"=3)[STRUCT]

		fnCOEFin<-paste0("step7-shockwaves/",struct[STRUCT] ,".LCOEF.",thisNumWave,".waves.RData")
		cat("[",STRUCT,"], ",NumWave,fnCOEFin,":",file.exists(fnCOEFin),"\n")
		rm(list=c("LHEAD_ALL", "LHEAD_ECM", "LHEAD_Ligands", "LHEAD_Receptors",  "LHEAD_TF", 
				 "LCOEF_ALL", "LCOEF_ECM", "LCOEF_Ligands", "LCOEF_Receptors", "LCOEF_TF"))
		load(fnCOEFin)
		LIST2<-list("Ligands"=LCOEF_Ligands, "Receptors"=LCOEF_Receptors, "TF"=LCOEF_TF, "ECM"=LCOEF_ECM, "ALL_else"=LCOEF_ALL, "ALL"=LCOEF_ALL)
		LIST1<-lapply(LIST2, function(x){
			lapply(x, function(y){
				(y%>%filter(estimate>2 & q_value < 0.01) %>% arrange(desc(estimate)))$gene_short_name
			})
		})
		LWAVES<-LIST1$ALL

		bymedian2 <- reorder(OBJ@colData$BioName2, pseudotime(big_cds), mean)
		DAT.pseudo<-data.frame(bymedian2 , BioName2=OBJ@colData$BioName2, PSEDUOT=pseudotime(big_cds))
		gp<-ggplot(DAT.pseudo, aes(x = PSEDUOT, y = bymedian2 )) +
			geom_density_ridges(scale = 4) + 
			scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
			scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
			coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
			theme_ridges()
		#plot(gp)

	
##
		genes<-intersect(Wilson3$GeneSymb,unlist(LWAVES))
		ids<-rowData(OBJ)[match(genes, rowData(OBJ)[,2]),1]
		COUNTS<-exprs(OBJ)[ids,,drop=F]
		COUNTS2<-apply(COUNTS,1,function(x)x/size_factors(OBJ))
		COUNTS3<-apply(COUNTS2,2,function(x){
			sapply(split(x, OBJ@colData$BioName2), mean, na.rm=T)
		})
		COUNTS4<-t(log2(COUNTS3+1))
		rownames(COUNTS4)<-genes
		COUNTS5<-t(apply(COUNTS4,1,function(x)x/max(x)))
		MEMB<-sapply(genes, function(x){
			round(mean(which(sapply(LWAVES, function(y){
				x%in%y
			}))))
		})
		ORDm<-order(MEMB)
		PSTcor<-cor(COUNTS2,pseudotime(OBJ))
		Lind<-sapply(sort(unique(MEMB)), function(i){
			indi<-which(MEMB==i)
			indi[order(PSTcor[indi])]
		})
		ORDt<-unlist(Lind)
		blues_fun <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))

		gplots::heatmap.2(COUNTS5[ORDt, ], mar=c(16,8), trace="none", dendrogram="none", Colv= F, Rowv= F, col=blues_fun(n=11), 
			main=paste0(STRUCT))
		diff(PSTcor[ORDt])
		image(t(as.matrix(rev(MEMB[ORDt]))), col=hue_pal()(max(MEMB)))
	}

dev.off()






