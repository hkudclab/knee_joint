library(ggplot2)
library(ggalluvial)
library(RColorBrewer)


files.successive<-list.files(path="successiveTP", pattern="successiveTP.noCC.*.RData", full.names=T)
for(i in 1:length(files.successive))
	load(files.successive[i])
load("../DATA/KJ.combCore.meta.RData")
LIDs<-sapply(split(rownames(CoreMetNonCycling), CoreMetNonCycling$orig.ident), function(x)gsub("_.*","",x))


GETMATflow<-function(OBJ, thisTP="KJ125"){
	ind125<-which(OBJ@meta.data$orig.ident==thisTP)

	MATdist1213<-as.matrix(dist(OBJ@reductions$tsne@cell.embeddings[1:5,]))
	ALLCLUST<-as.character(sort(unique(OBJ@meta.data$seurat_clusters)))

	MATflow<-matrix(NA,0,5)
	colnames(MATflow)<-c("i", "cellIDfrom", "cellIDto", "OldClustThisTP", "FLOWto")

	LIDto<-list()
	LIDto2<-list()

	for(i in 1:length(ALLCLUST)){
		indi<-which(OBJ@meta.data$seurat_clusters==ALLCLUST[i])
		MATdist1213.i<-as.matrix(dist(OBJ@reductions$tsne@cell.embeddings[indi,]))
		indi2<-which(OBJ@meta.data$orig.ident[indi]==thisTP)
		indi2.next<-which(OBJ@meta.data$orig.ident[indi]!=thisTP)
		OldClustThisTP<-OBJ[["OldClust"]][indi, 1][indi2]
		OldClustNextTP<-OBJ[["OldClust"]][indi, 1][indi2.next]

		if(length(indi2)>0 & length(indi2.next)>0){
			MATdist1213.i2<-MATdist1213.i[indi2, indi2.next, drop=F]
			CELLidsTO<-colnames(OBJ)[indi][indi2.next]
			FLOWto<-apply(MATdist1213.i2, 1, function(x){
				ORD3<-order(x)[1:3]
				names(which.max(table(OldClustNextTP[ORD3])))
			})
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
			LIDto[[i]]<-FLOWtoCell
			LIDto2[[i]]<-FLOWtoCell
			thisTemp<-cbind(i=rep(i, length(OldClustThisTP)), 
				cellIDfrom=colnames(OBJ)[indi][indi2],
				cellIDto=FLOWtoCell,
				OldClustThisTP, FLOWto)
			MATflow<-rbind(MATflow, thisTemp)
		}
	}
	return(MATflow)
}

MATflow1213<-GETMATflow(KJ1213.5, thisTP="KJ125")
MATflow1314<-GETMATflow(KJ1314.5, thisTP="KJ135")
MATflow1415<-GETMATflow(KJ1415.5, thisTP="KJ145")
MATflow1518<-GETMATflow(KJ1518.5, thisTP="KJ155")
MATflow18P5<-GETMATflow(KJ18P5.5, thisTP="KJ185")
MATflow18MCP5<-GETMATflow(KJ18MCP5.5, thisTP="KJ185")


save(MATflow1213, 
	MATflow1314,
	MATflow1415,
	MATflow1518,
	MATflow18P5,
	file="successiveTP/MATflow.noCC.E12.to.P5.RData")

#############################################
pdf("successiveTP/alluvial.chart.successive.timepoints.noCC.pdf", width=6)
	LMAT<-list("E12.5 -> E13.5"=MATflow1213, 
			"E13.5 -> E14.5"=MATflow1314, 
			"E14.5 -> E15.5"=MATflow1415, 
			"E15.5 -> E18.5"=MATflow1518, 
			"E18.5 -> EP5.MC"=MATflow18MCP5,
			"E18.5 -> EP5.5"=MATflow18P5)

	LMAT2<-sapply(LMAT, function(x){
		X4<-sapply(strsplit(x[,4],"_"), function(y){
			paste0(y[1],"_",formatC(as.integer(y[2]), width = 2, format = "d", flag = "0"))
		})
		#print(head(X4))
		X5<-sapply(strsplit(x[,5],"_"), function(y){
			paste0(y[1],"_",formatC(as.integer(y[2]), width = 2, format = "d", flag = "0"))
		})
		newDat<-data.frame(x[,1:3], X4, X5)
		print(head(newDat))
		return(newDat)
	})


	for(i in 1:6){
		DATi<-as.data.frame(table(LMAT2[[i]][,4], LMAT2[[i]][,5]))
		p1<-ggplot(DATi,
      			aes(y = Freq, axis1 = Var1, axis2 = Var2)) +
		  	geom_alluvium(aes(fill = Var1), width = 1/12) +
		 	geom_stratum(width = 1/12, fill = "grey", color = "white") +
		  	geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
		  	scale_fill_discrete( colorRampPalette(brewer.pal(9, "Set1"))(length(table(MATflow[,2])))) +
			ggtitle(names(LMAT)[i])
		plot(p1)
	}
dev.off()





