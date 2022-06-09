TRENDPLOT<-function(CDS, LIST1, LIST2, NumWave){
	for(i in 1:NumWave){
		for(j in 1:4){

			genei<-LIST1[[j]][[i]]
			print(genei)
			if(length(genei)==0){
				plot(0, main=paste0("[Wave-",i,"]; [", names(LIST2)[j],"]"), type="n")
				mtext(paste0("average signals of [",length(genei),"] genes"))
				next
			}
			IDi<-rowData(CDS)[which(rowData(CDS)[,2] %in% genei),1]
			yy<-rowMeans(apply(exprs(CDS)[IDi,,drop=F],1,function(x)x/size_factors(CDS)), na.rm=T)

			yMat<-apply(exprs(CDS)[IDi,,drop=F],1,function(x)x/size_factors(CDS))
			yy<-rowMeans(apply(yMat,2,function(x)x/max(x)), na.rm=T)

			pseu<-pseudotime(CDS)
			fitL<-loess(yy~pseu,  span = 0.5)

			predi<-predict(fitL, newdata=sort(pseu), se=T)

			plot(pseu, yy, pch=18, main=paste0("Total: ",NumWave," waves. This Wave-",i,". [", names(LIST1)[j],"]"))
			lines(sort(pseu), predi$fit, col="green", lwd=4)

			COEFij<-LIST2[[j]][[i]]
			estimate<-COEFij$estimate

			plot(estimate, -log10(COEFij$q_value), xlab="Estimate", ylab="-log10(q value)", main=paste0("[Wave-",i,"]; [", names(LIST2)[j],"]"),
				xlim=c(0,max(estimate)))
			mtext(paste0("average signals of [",length(IDi),"] genes"))
			legend("topleft", pch=16, col=c("#F8766D","#619CFF"), 
				legend=c("top positively associated genes", "top negatively associated genes (but will not be shown as we are concerned with the pos only)"))
			indpos<-which(COEFij$estimate>5 | (COEFij$estimate>1 & COEFij$q_value<1e-20))
			indneg<-which(COEFij$estimate< -5 |(COEFij$estimate< -1 & COEFij$q_value<1e-20))
			if(length(indpos)>0){
				points(COEFij$estimate[indpos], -log10(COEFij$q_value[indpos]), pch=16, col="#F8766D")
				text(COEFij$estimate[indpos], -log10(COEFij$q_value[indpos]), COEFij$gene_short_name[indpos], col=2, xpd=T)
			}
			if(length(indneg)>0){
				points(COEFij$estimate[indneg], -log10(COEFij$q_value[indneg]), pch=16, col="#619CFF")
				text(COEFij$estimate[indneg], -log10(COEFij$q_value[indneg]), COEFij$gene_short_name[indneg], col="blue", xpd=T)
			}

		}

	}
	for(genei in c("Lgr5", "Gdf5", "Sox9", "Scx")){
		IDi<-rowData(CDS)[which(rowData(CDS)[,2] %in% genei),1]
		yy<-exprs(CDS)[IDi,]/size_factors(CDS)

		pseu<-pseudotime(CDS) 
		yy<-yy
		fitL<-loess(yy~pseu,  span = 0.5)

		predi<-predict(fitL, newdata=sort(pseu), se=T)
		plot(pseu, yy, pch=18, main=genei)
		lines(sort(pseu), predi$fit, col="green", lwd=4)

		yy2<-yy[yy>0]
		pseu2<-pseu[yy>0]
		fitL2<-loess(yy2~pseu2,  span = 0.5)
		predi<-predict(fitL2, newdata=sort(pseu2), se=T)
		lines(sort(pseu2), predi$fit, col="magenta", lwd=4)
		legend("top", col=c("green","magenta"),lty=1, lwd=2, legend=c("fitted with zeros", "fitted with positive only"))
	}

}

TRENDPLOT2<-function(CDS, LIST1, LIST2, NumWave){
	for(j in 1:4){
		for(i in 1:NumWave){


			genei<-LIST1[[j]][[i]]
			genei<-LIST1[[j]][[i]]
			print(genei)


			IDi<-rowData(CDS)[which(rowData(CDS)[,2] %in% genei),1]
			if(length(IDi)==0){
				plot(0, main=paste0("[Wave-",i,"]; [", names(LIST2)[j],"]"), type="n")
				mtext(paste0("average signals of [",length(IDi),"] genes"))
				next
			}
			yy<-rowMeans(apply(exprs(CDS)[IDi,,drop=F],1,function(x)x/size_factors(CDS)), na.rm=T)

			yMat<-apply(exprs(CDS)[IDi,,drop=F],1,function(x)x/size_factors(CDS))
			yy<-rowMeans(apply(yMat,2,function(x)x/max(x)), na.rm=T)

			pseu<-pseudotime(CDS)
			fitL<-loess(yy~pseu,  span = 0.5)

			predi<-predict(fitL, newdata=sort(pseu), se=T)

			plot(pseu, yy, pch=18, main=paste0("Total: ",NumWave," waves. This Wave-",i,". [", names(LIST1)[j],"]"))
			lines(sort(pseu), predi$fit, col="green", lwd=4)

			if(0){
				COEFij<-LIST2[[j]][[i]]
				estimate<-COEFij$estimate

				plot(estimate, -log10(COEFij$q_value), xlab="Estimate", ylab="-log10(q value)", main=paste0("[Wave-",i,"]; [", names(LIST2)[j],"]"),
					xlim=c(0,max(estimate)))
				mtext(paste0("average signals of [",length(IDi),"] genes"))
				legend("topleft", pch=16, col=c("#F8766D","#619CFF"), 
					legend=c("top positively associated genes", "top negatively associated genes (but will not be shown as we are concerned with the pos only)"))
				indpos<-which(COEFij$estimate>5 | (COEFij$estimate>1 & COEFij$q_value<1e-20))
				indneg<-which(COEFij$estimate< -5 |(COEFij$estimate< -1 & COEFij$q_value<1e-20))
				if(length(indpos)>0){
					points(COEFij$estimate[indpos], -log10(COEFij$q_value[indpos]), pch=16, col="#F8766D")
					text(COEFij$estimate[indpos], -log10(COEFij$q_value[indpos]), COEFij$gene_short_name[indpos], col=2, xpd=T)
				}
				if(length(indneg)>0){
					points(COEFij$estimate[indneg], -log10(COEFij$q_value[indneg]), pch=16, col="#619CFF")
					text(COEFij$estimate[indneg], -log10(COEFij$q_value[indneg]), COEFij$gene_short_name[indneg], col="blue", xpd=T)
				}
			}
		}

	}
	for(genei in c("Lgr5", "Gdf5", "Sox9", "Scx")){
		IDi<-rowData(CDS)[which(rowData(CDS)[,2] %in% genei),1]
		yy<-exprs(CDS)[IDi,]/size_factors(CDS)

		pseu<-pseudotime(CDS) 
		yy<-yy
		fitL<-loess(yy~pseu,  span = 0.5)

		predi<-predict(fitL, newdata=sort(pseu), se=T)
		plot(pseu, yy, pch=18, main=genei)
		lines(sort(pseu), predi$fit, col="green", lwd=4)

		yy2<-yy[yy>0]
		pseu2<-pseu[yy>0]
		fitL2<-loess(yy2~pseu2,  span = 0.5)
		predi<-predict(fitL2, newdata=sort(pseu2), se=T)
		lines(sort(pseu2), predi$fit, col="magenta", lwd=4)
		legend("top", col=c("green","magenta"),lty=1, lwd=2, legend=c("fitted with zeros", "fitted with positive only"))
	}

}

TRENDPLOT_loess<-function(CDS, LIST2, NumWave){
	LIST1<-lapply(LIST2, function(x){
		lapply(x, function(y){
			(y%>%filter(estimate>2 & q_value < 0.01) %>% arrange(desc(estimate)))$gene_short_name
		})
	})
	
	for(i in 1:NumWave){
		LIST1[[5]][[i]]<-setdiff(LIST1[[5]][[i]], c(LIST1[[1]][[i]], LIST1[[2]][[i]], LIST1[[3]][[i]], LIST1[[4]][[i]]))
	}
	for(j in 6:6){
		Lpred<-list()
		Lyy<-list()
		cat(j,"\t",date(),"\n")
		print(sapply(LIST1[[j]], length))
		for(i in 1:NumWave){
			genei<-LIST1[[j]][[i]]

			IDi<-rowData(CDS)[which(rowData(CDS)[,2] %in% genei),1]
			if(length(IDi)==0){
				#plot(0, main=paste0("[Wave-",i,"]; [", names(LIST2)[j],"]"), type="n")
				#mtext(paste0("average signals of [",length(IDi),"] genes"))
				next
			}
			yy<-rowMeans(apply(exprs(CDS)[IDi,,drop=F],1,function(x)x/size_factors(CDS)), na.rm=T)

			yMat<-apply(exprs(CDS)[IDi,,drop=F],1,function(x)x/size_factors(CDS))
			yy<-rowMeans(apply(yMat,2,function(x)x/max(x)), na.rm=T)

			pseu<-pseudotime(CDS)
			fitL<-loess(yy~pseu,  span = 0.5)

			predi<-predict(fitL, newdata=sort(pseu), se=T)
			Lpred[[i]]<-predi
			Lyy[[i]]<-yy

		}
		Lpred2<-Lpred[sapply(Lpred, length)>0]
		NumWave2<-sum(sapply(Lpred, length)>0)
		DATrib<-data.frame(waves=as.character(rep(seq(NumWave)[sapply(Lpred, length)>0], rep(length(pseu),NumWave2))),
			PST=rep(sort(pseu), NumWave2), 
			fit=as.vector(sapply(Lpred2, function(x)x$fit)),
			yy=unlist(Lyy),
			sefit=as.vector(sapply(Lpred2, function(x)x$se.fit)))


		g1<-ggplot(DATrib, aes(x=PST, y = yy))
		for(i in 1:NumWave)
			g1<-g1 + geom_ribbon(data=DATrib%>%filter(waves==i),
					aes(ymin = fit - 10*sefit, ymax = fit + 10*sefit,fill=waves), alpha=0.25) +
  				geom_line(data=DATrib%>%filter(waves==i),
					aes(y = fit, color=waves), size=1)# + geom_point(aes(color=waves))
			#facet_wrap(vars(waves),ncol =1)
		str1<-paste0("[",names(LIST1)[j],"] in [",STRUCT,"]")
		str2<-paste0("[",paste0("wave-",seq(NumWave),": ",sapply(LIST1[[j]], length), collapse=" genes], ["),"]")
		g1<-g1+ggtitle(paste0(str1, "\n", str2))
		plot(g1)

		str1<-sapply(LIST1[[j]], function(x){
			paste0(x, collapse=", ")
		})
		matGenes<-cbind("structure"=STRUCT, "gene category"=names(LIST2)[j], 
			"waves"=paste0("wave-", seq(NumWave)), 
			"num genes"=sapply(LIST1[[j]], length), "genes (in descending order of assoc)"=str1)

		tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
		matGenes[,5]= str_wrap(matGenes[,5], 90)
		if(sum(sapply(LIST1[[j]], length))>250){
			tt <- gridExtra::ttheme_default(
    				core = list(fg_params=list(cex = 1/2)),
    				colhead = list(fg_params=list(cex = 1.0)),
    				rowhead = list(fg_params=list(cex = 1.0)))
			matGenes[,5]= str_wrap(matGenes[,5], 150)
		}
		if(sum(sapply(LIST1[[j]], length))>800){
			tt <- gridExtra::ttheme_default(
    				core = list(fg_params=list(cex = 1/3)),
    				colhead = list(fg_params=list(cex = 1/2)),
    				rowhead = list(fg_params=list(cex = 1/2)))
			matGenes[,5]= str_wrap(matGenes[,5], 240)
		}
		if(sum(sapply(LIST1[[j]], length))>1200){
			tt <- gridExtra::ttheme_default(
    				core = list(fg_params=list(cex = 0.3)),
    				colhead = list(fg_params=list(cex = 1/3)),
    				rowhead = list(fg_params=list(cex = 1/3)))
			matGenes[,5]= str_wrap(matGenes[,5], 350)
		}
		grid.newpage()
		grid.table(matGenes, theme=tt )

	}
	for(genei in c("Lgr5", "Gdf5", "Sox9", "Scx")){
		IDi<-rowData(CDS)[which(rowData(CDS)[,2] %in% genei),1]
		yy<-exprs(CDS)[IDi,]/size_factors(CDS)

		pseu<-pseudotime(CDS) 
		yy<-yy
		#fitL<-loess(yy~pseu,  span = 0.5)

		#predi<-predict(fitL, newdata=sort(pseu), se=T)

		#DATall<-data.frame(PST=sort(pseu),yy, fit=predi$fit, sefit=predi$se.fit, zeros="Y")
		#g1<-ggplot(DATall, aes(x=PST, y=yy))+geom_ribbon(data=DATall,
		#		aes(ymin = fit - 3*sefit, ymax = fit + 3*sefit, fill=zeros),  alpha=0.25) +
  		#	geom_line(data=DATall,aes(y = fit, color=zeros), size=1) #+
			#geom_point(color="grey")

		yy2<-yy[yy>0]
		pseu2<-pseu[yy>0]
		fitL2<-loess(yy2~pseu2,  span = 0.5)
		predi2<-predict(fitL2, newdata=sort(pseu2), se=T)

		DATpos<-data.frame(PST=sort(pseu2), yy=yy2[order(pseu2)], fit=predi2$fit, sefit=predi2$se.fit, zeros="N")
		g1<-ggplot(DATpos, aes(x=PST, y=yy))
		g1<-g1+geom_ribbon(data=DATpos,
				aes(ymin = fit - 3*sefit, ymax = fit + 3*sefit, fill=zeros), alpha=0.25) +
  			geom_line(data=DATpos,aes(y = fit, color=zeros), size=1) + #geom_point(color="grey") +
			ggtitle(genei)
		plot(g1)
	}
}


TRENDPLOT_gam<-function(CDS, LIST2, NumWave, BS="tp"){
	SEpt<-qnorm(0.995)
	LIST1<-lapply(LIST2, function(x){
		lapply(x, function(y){
			(y%>%filter(estimate>2 & q_value < 0.01) %>% arrange(desc(estimate)))$gene_short_name
		})
	})
	
	for(i in 1:NumWave){
		LIST1[[5]][[i]]<-setdiff(LIST1[[5]][[i]], c(LIST1[[1]][[i]], LIST1[[2]][[i]], LIST1[[3]][[i]], LIST1[[4]][[i]]))
	}
	for(j in 1:5){
		#break
		Lpred<-list()
		for(i in 1:NumWave){
			genei<-LIST1[[j]][[i]]

			IDi<-rowData(CDS)[which(rowData(CDS)[,2] %in% genei),1]
			if(length(IDi)==0){
				#plot(0, main=paste0("[Wave-",i,"]; [", names(LIST2)[j],"]"), type="n")
				#mtext(paste0("average signals of [",length(IDi),"] genes"))
				next
			}
			yy<-rowMeans(apply(exprs(CDS)[IDi,,drop=F],1,function(x)x/size_factors(CDS)), na.rm=T)

			yMat<-apply(exprs(CDS)[IDi,,drop=F],1,function(x)x/size_factors(CDS))
			yy<-rowMeans(apply(yMat,2,function(x)x/max(x)), na.rm=T)

			pseu<-pseudotime(CDS)
			#fitL<-loess(yy~pseu,  span = 0.5)

			fitGAM<-gam(yy~s(pseu, bs=BS),  span = 0.5, family=quasipoisson())
			prediGAM<-predict(fitGAM, newdata=data.frame(pseu=sort(pseu)),  se=T)

			#predi<-predict(fitL, newdata=sort(pseu), se=T)
			Lpred[[i]]<-prediGAM

		}
		Lpred2<-Lpred[sapply(Lpred, length)>0]
		NumWave2<-sum(sapply(Lpred, length)>0)
		DATrib<-data.frame(waves=as.character(rep(seq(NumWave)[sapply(Lpred, length)>0], rep(length(pseu),NumWave2))),
			PST=rep(sort(pseu), NumWave2), 
			fit=as.vector(sapply(Lpred2, function(x)x$fit)),
			sefit=as.vector(sapply(Lpred2, function(x)x$se.fit)))


		g1<-ggplot(DATrib, aes(x=PST))
		for(i in 1:NumWave)
			g1<-g1 + geom_ribbon(data=DATrib%>%filter(waves==i),
					aes(ymin = exp(fit - SEpt*sefit), ymax = exp(fit + SEpt*sefit),fill=waves), alpha=0.25) +
  				geom_line(data=DATrib%>%filter(waves==i),
					aes(y = exp(fit), color=waves), size=1)
			#facet_wrap(vars(waves),ncol =1)
		str1<-paste0("[",names(LIST1)[j],"] in [",STRUCT,"]")
		str2<-paste0("[",paste0("wave-",seq(NumWave),": ",sapply(LIST1[[j]], length), collapse=" genes], ["),"]")
		g1<-g1+ggtitle(paste0(str1, "\n", str2))
		plot(g1)

		str1<-sapply(LIST1[[j]], function(x){
			paste0(x, collapse=", ")
		})
		matGenes<-cbind("structure"=STRUCT, "gene category"=names(LIST2)[j], 
			"waves"=paste0("wave-", seq(NumWave)), 
			"num genes"=sapply(LIST1[[j]], length), "genes (in descending order of assoc)"=str1)

		tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
		matGenes[,5]= str_wrap(matGenes[,5], 90)
		if(sum(sapply(LIST1[[j]], length))>250){
			tt <- gridExtra::ttheme_default(
    				core = list(fg_params=list(cex = 1/2)),
    				colhead = list(fg_params=list(cex = 1.0)),
    				rowhead = list(fg_params=list(cex = 1.0)))
			matGenes[,5]= str_wrap(matGenes[,5], 150)
		}
		if(sum(sapply(LIST1[[j]], length))>800){
			tt <- gridExtra::ttheme_default(
    				core = list(fg_params=list(cex = 1/3)),
    				colhead = list(fg_params=list(cex = 1.0)),
    				rowhead = list(fg_params=list(cex = 1.0)))
			matGenes[,5]= str_wrap(matGenes[,5], 200)
		}
		if(sum(sapply(LIST1[[j]], length))>1200){
			tt <- gridExtra::ttheme_default(
    				core = list(fg_params=list(cex = 0.3)),
    				colhead = list(fg_params=list(cex = 1/3)),
    				rowhead = list(fg_params=list(cex = 1/3)))
			matGenes[,5]= str_wrap(matGenes[,5], 350)
		}
		grid.newpage()
		grid.table(matGenes, theme=tt )

	}
	GENES<-c("Lgr5", "Gdf5", "Sox9", "Scx")
	for(i in 1:4){
		genei<- GENES[i]
		cat(i,"\t",genei,"\n")
		IDi<-rowData(CDS)[which(rowData(CDS)[,2] %in% genei),1]
		yy<-exprs(CDS)[IDi,]/size_factors(CDS)

		pseu<-pseudotime(CDS) 
		yy<- yy/max(yy)

		fitGAM<-gam(yy~s(pseu, bs=BS),  span = 0.5, family=quasipoisson())
		prediGAM<-predict(fitGAM, newdata=data.frame(pseu=sort(pseu)),  se=T)

		DATall<-data.frame(PST=sort(pseu),yy, fit=prediGAM$fit, sefit=prediGAM$se.fit, genes=genei)
		if(i==1)
			g1<-ggplot(DATall, aes(x=PST, y=yy))

		g1<- g1 + geom_ribbon(data=DATall,
				aes(ymin = exp(fit - SEpt*sefit), ymax = exp(fit + SEpt*sefit), fill=genes),  alpha=0.25) +
  			geom_line(data=DATall,aes(y = exp(fit), color=genes), size=1)
	}
	g1<- g1 +ylim(c(0,0.5)) + xlim(range(pseudotime(CDS)))
	g1<-g1+ggtitle(paste0("[", paste0(GENES, collapse=", "),"]\n fitted with zero values"))
	plot(g1)

	for(i in 1:4){
		genei<- GENES[i]
		IDi<-rowData(CDS)[which(rowData(CDS)[,2] %in% genei),1]
		yy<-exprs(CDS)[IDi,]/size_factors(CDS)
		yy<- yy/max(yy)

		pseu<-pseudotime(CDS) 
		yy2<-yy[yy>0]
		pseu2<-pseu[yy>0]
		#fitL2<-loess(yy2~pseu2,  span = 0.5)
		#predi2<-predict(fitL2, newdata=data.frame(pseu2=sort(pseu2)), se=T)

		fitGAM2<-gam(yy2~s(pseu2, bs=BS),  span = 0.5, family=quasipoisson())
		prediGAM2<-predict(fitGAM2, newdata=data.frame(pseu2=sort(pseu2)),  se=T)
		DATpos<-data.frame(PST=sort(pseu2), yy=yy2[order(pseu2)], fit=prediGAM2$fit, sefit=prediGAM2$se.fit, genes=genei)

		if(i==1)
			g1<-ggplot(DATpos, aes(x=PST, y=yy))

		g1<-g1+geom_ribbon(data=DATpos,
				aes(ymin = exp(fit - SEpt*sefit), ymax = exp(fit + SEpt*sefit), fill=genes), alpha=0.25) +
  			geom_line(data=DATpos,aes(y = exp(fit), color=genes), size=1) 
	}
	g1<- g1 +ylim(c(0,0.5)) + xlim(range(pseudotime(CDS)))
	g1<-g1+ggtitle(paste0("[", paste0(GENES, collapse=", "),"]\n fitted with non-zero values (i.e. positive values), only"))
	plot(g1)

}
