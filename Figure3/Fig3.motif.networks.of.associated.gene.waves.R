library(RcisTarget)
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(igraph)
library(scales)
library(ggridges)
library(RColorBrewer)

library(packcircles)
library(ggforce)

library(gridExtra)
library(grid)

library(arcdiagram)

library(RCy3)

load("P:/OneDrive - The University Of Hong Kong/resources/TF.CD.receptors.human.mouse.expanded.RData")


data(motifAnnotations_mgi)
motifRankings <- importRankings("P:/databases/scenic/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
motifRankings <- importRankings("P:/databases/scenic/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")


FnData<-list("AC"="step7-shockwaves/combined.and.pseudo.and.wave-sigs.AC.RData", 
		"SZ"="step7-shockwaves/combined.and.pseudo.and.wave-sigs.SZ.RData", 
		"Liga"="step7-shockwaves/combined.and.pseudo.and.wave-sigs.Lig.RData", 
		"MCLig"="step7-shockwaves/combined.and.pseudo.and.wave-sigs.MCLig.RData")

Wilson<-readxl::read_excel("step7-shockwaves/20220513 Wilson selected 5 genes based on trendlines.xlsx", sheet = "Sheet1")
struct<-c("AC"="AC", "SZ"="SZ", "Liga"="Lig", "MCLig"="MCLig")

load("step7-shockwaves/cds.and.pseudo.and.wave-sigs.AC.RData")

L.DATNET<-list()
for(STRUCT in c("AC", "SZ", "Liga", "MCLig")){
	#load(FnData[[STRUCT]])
	#STRUCT<-"SZ"

	thisNumWave<-c("AC"=4, "SZ"=3, "Liga"=5, "MCLig"=3)[STRUCT]

	fnCOEFin<-paste0("step7-shockwaves/",struct[STRUCT] ,".LCOEF.",thisNumWave,".waves.RData")
	cat("[",STRUCT,"], ",thisNumWave,fnCOEFin,":",file.exists(fnCOEFin),"\n")
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

	#motifEnrichmentTable_wGenes <- cisTarget(LHEAD_ALL, motifRankings,
	  #                         motifAnnot=motifAnnotations_mgi)

	names(LWAVES)<-paste0("wave", seq(length(LWAVES)))
	motifEnrichmentTable_wGenes <- cisTarget(LWAVES, motifRankings,  nesThreshold = 2,
							   motifAnnot=motifAnnotations_mgi)

	Ltargets<-strsplit(motifEnrichmentTable_wGenes$enrichedGenes, ";")


	ALLTFs<-gsub(" .*|;","",motifEnrichmentTable_wGenes$TF_highConf)
	AllTargets<-strsplit(motifEnrichmentTable_wGenes$enrichedGenes, ";")

	#intersect(setdiff(ALLTFs, ""), unique(unlist(LHEAD_TF)))
	LTF0<-split(ALLTFs, 
		motifEnrichmentTable_wGenes$geneSet)
	LTF1<-sapply(LTF0, function(x)setdiff(x,""))
	LTF2<-sapply(seq(length(LTF1)), function(i){sort(intersect(LTF1[[i]], unlist(LIST1$TF[seq(i)])))})

	Ltargets<-split(motifEnrichmentTable_wGenes$enrichedGenes, motifEnrichmentTable_wGenes$geneSet)
	LAUC<-split(motifEnrichmentTable_wGenes$AUC, motifEnrichmentTable_wGenes$geneSet)

	Ltargetlen<-sapply(seq(length(LTF2)), function(i){
		sapply(LTF2[[i]], function(x){
			indi<-which(LTF0[[i]]%in%x)
			length(sort(unique(unlist(strsplit(Ltargets[[i]][indi], ";")))))
		})
	})
	Lauc<-sapply(seq(length(LTF2)), function(i){
		sapply(LTF2[[i]], function(x){
			indi<-which(LTF0[[i]]%in%x)
			mean(LAUC[[i]][indi])
		})
	})

	Lprior<-sapply(seq(length(LTF2)), function(i){
		pos<-c()
		for(j in 1:length(LTF2[[i]])){
			etf <-  LTF2[[i]][j]
			pos[j]<-which(sapply(LIST1$TF[seq(i)], function(tfs){
				etf%in%tfs
			}))
		}
		pos
	})
	waves<-rep(seq(length(Ltargetlen)), sapply(Ltargetlen, length))
	DATNET<-data.frame(TFenriched=names(unlist(Ltargetlen)),AUC=unlist(Lauc), numTar=unlist(Ltargetlen), from=unlist(Lprior), to=waves)
	L.DATNET[[STRUCT]]<-DATNET
}

#save(L.DATNET, file="step7-shockwaves/networks/L.DATNET.RData")



nosology<-readxl::read_excel("skeletal_dysplaysia_nosology_2019.xlsx", sheet = "Sheet2")
nosology<-nosology[!is.na(nosology$Genes), ]
LNos<-strsplit(nosology$Genes, ",")
for(i in  which(is.na(nosology$nameOfDisorder))){
	nosology$nameOfDisorder[i]<-nosology$nameOfDisorder[i-1]
}
disorders<-rep(nosology$nameOfDisorder, sapply(LNos, length))
NOSOLOGY<-data.frame(disorders, gene=unlist(LNos))
NOSOLOGY[!grepl("[A-Z]",NOSOLOGY[, 2]),]
NOSOLOGY2<-NOSOLOGY[grepl("[A-Z]",NOSOLOGY[, 2]), ]
table(NOSOLOGY2[,2]%in%human.mouse.2020$Gene.name)
setdiff(NOSOLOGY2[,2], human.mouse.2020$Gene.name)

indHM<-match(NOSOLOGY2[,2], human.mouse.2020$Gene.name)
NOSOLOGY2[is.na(indHM),2]
MGenes<-as.character(human.mouse.2020$Mouse.gene.name[indHM])
MGenes[which(MGenes=="")]<-str_to_title(NOSOLOGY2[,2][which(MGenes=="")])

NOSOLOGY2[["mouseGenes"]]<-MGenes

pdf(paste0("step7-shockwaves/networks/arc.diagram.FOUR-STRUCTURES.nosology.only.pdf"),  height=8,  width=24)
for(STRUCT in c("AC", "SZ", "Liga", "MCLig")){
	#STRUCT<-"Liga"
	DATNET<-L.DATNET[[STRUCT]]
######################
	gobjects<-list()
	#pdf(paste0("step7-shockwaves/networks/arc.diagram.",STRUCT,".pdf"),  height=8,  width=24)
	for(j in 1:max(DATNET$from)){
		DATNET3<-DATNET[DATNET$from==j,]

		startingTFs<-sort(sapply(split(DATNET3$numTar, DATNET3$TFenriched), sum),decreasing=T)
		
		packing <- circleProgressiveLayout(sqrt(startingTFs), sizetype='area')
		data <- cbind(TF=names(startingTFs), SIZE=log2(startingTFs), packing)
		
		dat.gg <- circleLayoutVertices(packing, npoints=50)
		if(j==1)
			DIA<-diff(range(data$x))*2

		data$x<-data$x + DIA*(DATNET3$from[1]-1)
		dat.gg$x<-dat.gg$x + DIA*(DATNET3$from[1]-1)

		Mx<-mean(data$x)
		My<-mean(data$y)

		for(i in 1:nrow(DATNET3)){
			xi<-data[DATNET3[i,1],"x"]
			yi<-data[DATNET3[i,1],"y"]

			x2<-Mx + DIA+ (DATNET3[i,"to"]-DATNET3[i,"from"])*DIA/2
			y2<-My
			
			vecX<- x2-xi
			vecY<- y2-yi
			angle<-  pi/3 #runif(1, min=pi/4, max=pi/3)
			newX <- vecX* cos(angle) - vecY* sin(angle)
			newY <- vecX* sin(angle) + vecY* cos(angle)
			x0<- xi + newX 
			y0<- yi+ newY 

			x02<-(xi+x2)/2
			y02<-(yi+y2)/2

			ratio<- 3/10 #runif(1)

			x00<-x02*ratio + x0*(1-ratio )
			y00<-y02*ratio + y0*(1-ratio )

			#x0<-Mx+DATNET3[i,"to"]*DIA/2
			#y0<- (yi+My)/2 - ((Mx-x0)^2 -(xi-x0)^2)/(yi-My)/2

			radi<- sqrt((x00-xi)^2 + (y00-yi)^2)
			deltaY<-yi-y00
			deltaX<-xi-x00
			tani<- deltaY/deltaX
			startRad<-  atan(tani)

			deltaY<-My-y00
			deltaX<-x2-x00
			tani<- deltaY/deltaX
			endRad<- atan(tani) 

			tmpi<-c(x0=x0,y0=y0, x2=x2, y2=y2, x02=x02, y02=y02, x00=x00, y00=y00, start=startRad, end=endRad, radius=radi, numTar=DATNET3[i,"numTar"], AUC=DATNET3[i,"AUC"])
			if(i==1)ARCS<-as.data.frame(t(as.matrix(tmpi)))
			if(i!=1)ARCS<-rbind(ARCS, tmpi)
		}
		ARCS<-data.frame(ARCS, ID=match(DATNET3$TFenriched, data$TF))
		if(j==1)
			RANGE<-range(c(dat.gg$x, ARCS$x2, max(ARCS$x2)+2.5*DIA))*1.1
		X2<-unique(ARCS$x2)
		Y2<-ARCS$y2[match(X2, ARCS$x2)]
		#DATA2<-data.frame(X2, Y2, texts=paste0("targets\nof\nwave-", -1+j+seq(length(X2))))
		DATA2<-data.frame(X2, Y2, texts=paste0("wave-", -1+j+seq(length(X2))))
		DATA2[nrow(DATA2)+1, 1:2]<-c(X2=Mx, Y2=max(data$y)*1.2)
		DATA2[nrow(DATA2), 3]<-paste0("TFs associated with wave-",j," and \nlinked to enriched motifs of subsequent waves")

		RANGEy<-range(c(dat.gg$y, ARCS$y2))*1.1
		data2<-data[which(data$TF %in% MGenes),,drop=F]
		g1<-ggplot() + 
			geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)),  alpha = 0.6) +
			#scale_size_continuous(range = c(3,5)) + 
			theme_void() + 
			theme(legend.position="none") +
			coord_equal()+ xlim(RANGE) + #ylim(RANGEy) +
			#geom_point(data=ARCS, aes(x0, y0)) +
			#geom_point(data=ARCS, aes(x2, y2), col="red") +
			#geom_point(data=ARCS, aes(x02, y02), col="green") +
			#geom_point(data=ARCS, aes(x00, y00), col="blue") +
			geom_arc(data=ARCS[,], aes(x0 = x00 , y0 = y00 , r = radius, start =3/2* pi  - start , end = pi/2 - end ,
				size=AUC*40, color=as.factor(ID)), alpha=0.5)+
			#geom_text(data = data, aes(x, y, size=6, label = TF), lineheight = .6) +
			scale_size_continuous(range = c(1,5)) + 
			geom_text(data = DATA2, aes(x=X2, y=Y2, label = texts), lineheight = .6) +
			ggtitle(paste0(STRUCT," [wave-",j,"] (total:",max(DATNET$from)," waves)" ))

		if(nrow(data2)>0){
			g1<-g1+geom_text(data = data2, aes(x, y, size=6, label = TF), lineheight = .6) 
		}
		plot(g1)

		if(nrow(data2)>0){
			pushViewport(viewport(x=0.75, y=.25,height=.5))
	
			NOSOLOGY3<-NOSOLOGY2[which(MGenes %in% data$TF), , drop=F]
			grid.table(NOSOLOGY3)
			popViewport()
		}
	}
}
dev.off()

		if(length(gobjects)==3)
			ggarrange(gobjects[[1]],gobjects[[2]],gobjects[[3]],ncol =1)+
				theme(plot.margin = margin(0.1,0.1,2,0.1, "cm")) 

		if(length(gobjects)==4)
			ggarrange(gobjects[[1]],gobjects[[2]],gobjects[[3]],gobjects[[4]],ncol =2)+
				theme(plot.margin = margin(0.1,0.1,2,0.1, "cm")) 
		if(length(gobjects)==5)
			ggarrange(gobjects[[1]],gobjects[[2]],gobjects[[3]],gobjects[[4]],gobjects[[5]],nrow =1)+
				theme(plot.margin = margin(0.1,0.1,2,0.1, "cm")) 

}





