load("successiveTP/MATflow.noCC.E12.to.P5.kNN.RData")

library(dplyr)
library(tidyverse)

library(ggplot2)
library(ggalluvial)
data(majors)
################################
GETWEIGHT<-function(MATflow.kNN, back=F){
	from<-5; to<-6
	if(back){from<-6; to<-5;}
	MATOUT<-as.matrix(table(MATflow.kNN[,from], MATflow.kNN[,to]))
	MATOUT[,]<-0
	for(i in 1:nrow(MATflow.kNN)){
		ROWN<-MATflow.kNN[i,from]
		COLN<-MATflow.kNN[i,to]
		MATOUT[ROWN, COLN]<-MATOUT[ROWN, COLN]+MATflow.kNN$WEIGHTs[i]
	}
	as.data.frame(as.table(MATOUT))%>%filter(Freq>0)
}
data.frame(table(MATflow1213.kNN[,5], MATflow1213.kNN[,6]))%>%filter(Freq>0)

data1213<-GETWEIGHT(MATflow1213.kNN)
data1314<-GETWEIGHT(MATflow1314.kNN)
data1415<-GETWEIGHT(MATflow1415.kNN)
data1518<-GETWEIGHT(MATflow1518.kNN)
data18P5<-GETWEIGHT(MATflow18P5.kNN)

data1213.back<-GETWEIGHT(MATflow1213.kNN.back, back=T)
data1314.back<-GETWEIGHT(MATflow1314.kNN.back, back=T)
data1415.back<-GETWEIGHT(MATflow1415.kNN.back, back=T)
data1518.back<-GETWEIGHT(MATflow1518.kNN.back, back=T)
data18P5.back<-GETWEIGHT(MATflow18P5.kNN.back, back=T)

LDATA<-list(data1213, data1314, data1415, data1518, data18P5)
getTransitDatFwd<-function(LDATA, CUTOFF=0.1){
	DATA<-LDATA[[1]]
	currentPtr<-2
	for(i in 2:5){
		print(i)
		Vari<-paste0("Var",i)
		Vari.p1<-paste0("Var",i+1)
		thisData<-LDATA[[i]]
		#DATA<-data.frame(DATA, CNT=as.vector(table(thisData$Var1)[as.character(DATA[[Vari]])]))%>%
		#	uncount(CNT)

		indFreq<-which(colnames(DATA)=="Freq")
		thisFreq<-paste0(colnames(DATA)[indFreq], i)
		colnames(DATA)[indFreq]<-thisFreq
		DATA[[Vari.p1]]<-""
		DATA[["Freq"]]<-0

		LEVELSi<-as.character(unique(DATA[[Vari]]))
		for(j in 1:length(LEVELSi)){
			indj<-which(thisData$Var1==LEVELSi[j])
			thisData[indj,]

			which(DATA[[Vari]]==LEVELSi[j])
			DATA2<-DATA
			for(k in which(DATA2[[Vari]]==LEVELSi[j])){
				if(k%%100==0)cat(k," ")
				if(k%%1000==0)cat("\n")
				datk<-data.frame(DATA2[k,],CNT=length(indj))%>%uncount(CNT)
				datk[[Vari.p1]]<-as.character(thisData$Var2[indj])
				datk$Freq<-round(datk[[thisFreq]]* thisData$Freq[indj]/sum(thisData$Freq[indj]), 2)
				DATA<-DATA[-k,]
				DATA<-rbind(datk, DATA)
			}
		}
		print(table(DATA$Freq>CUTOFF))
		DATA<-DATA[DATA$Freq>CUTOFF, ]
	}
	return(DATA)
}
LDATA.back<-list(data18P5.back, data1518.back, data1415.back, data1314.back, data1213.back)
getTransitDatBack<-function(LDATA.back, CUTOFF=0.1){
	DATA<-LDATA.back[[1]]
	currentPtr<-2
	for(i in 2:5){
		print(i)
		Vari<-paste0("Var",i)
		Vari.p1<-paste0("Var",i+1)
		thisData<-LDATA.back[[i]]
		#DATA<-data.frame(DATA, CNT=as.vector(table(thisData$Var1)[as.character(DATA[[Vari]])]))%>%
		#	uncount(CNT)

		indFreq<-which(colnames(DATA)=="Freq")
		thisFreq<-paste0(colnames(DATA)[indFreq], i)
		colnames(DATA)[indFreq]<-thisFreq
		DATA[[Vari.p1]]<-""
		DATA[["Freq"]]<-0

		LEVELSi<-as.character(unique(DATA[[Vari]]))
		for(j in 1:length(LEVELSi)){
			indj<-which(thisData$Var1==LEVELSi[j])
			thisData[indj,]

			which(DATA[[Vari]]==LEVELSi[j])
			DATA2<-DATA
			for(k in which(DATA2[[Vari]]==LEVELSi[j])){
				if(k%%100==0)cat(k," ")
				if(k%%1000==0)cat("\n")
				datk<-data.frame(DATA2[k,],CNT=length(indj))%>%uncount(CNT)
				datk[[Vari.p1]]<-as.character(thisData$Var2[indj])
				datk$Freq<-round(datk[[thisFreq]]* thisData$Freq[indj]/sum(thisData$Freq[indj]), 2)
				DATA<-DATA[-k,]
				DATA<-rbind(datk, DATA)
			}
		}
		print(table(DATA$Freq>CUTOFF))
		DATA<-DATA[DATA$Freq>CUTOFF, ]
	}
	return(DATA)
}
####################################


DATA0.1<-getTransitDatFwd(LDATA)
DATA<-getTransitDatFwd(LDATA, CUTOFF=5)

GETpropa<-function(DATAin){
	MATOUT<-as.matrix(table(DATAin$Var6, DATAin$Var1))
	MATOUT[,]<-0
	print(MATOUT)
	for(i in 1:nrow(DATAin)){
		ROWN<-DATAin$Var6[i]
		COLN<-DATAin$Var1[i]
		MATOUT[ROWN, COLN]<-MATOUT[ROWN, COLN]+DATAin$Freq[i]
	}
	as.table(MATOUT)
}
forward<-GETpropa(DATA0.1)
pdf("forward.prediction.pdf")
	for(i in 1:ncol(forward))
		pie(forward[,i])
dev.off()



DATA_long0 <- to_lodes_form(droplevels(DATA[DATA$Freq>0.2,]%>%mutate(GRP=Var1=="KJ125_0")),
                              key = "TP",
                              axes =c(1,2,4,6,8,10))%>%mutate(timepoint=paste0("E", gsub("_[0-9]+|KJ","",stratum)))%>%
				mutate(CLUST=gsub("^.*_","",stratum))

DATA_long8 <- to_lodes_form(droplevels(DATA[DATA$Freq>0.2,]%>%mutate(GRP=Var1=="KJ125_8")),
                              key = "TP",
                              axes =c(1,2,4,6,8,10))%>%mutate(timepoint=paste0("E", gsub("_[0-9]+|KJ","",stratum)))%>%
				mutate(CLUST=gsub("^.*_","",stratum))


DATA_long56 <- to_lodes_form(droplevels(DATA[DATA$Freq>0.2,]%>%mutate(GRP=Var1=="KJ125_5"|Var1=="KJ125_6")),
                              key = "TP",
                              axes =c(1,2,4,6,8,10))%>%mutate(timepoint=paste0("E", gsub("_[0-9]+|KJ","",stratum)))%>%
				mutate(CLUST=gsub("^.*_","",stratum))

pdf("all-in-one.kNN.forward.pdf", width=16)
	ggplot(data = DATA_long0,
		   aes(x = timepoint, stratum = stratum, alluvium = alluvium,
			   y = Freq, label = CLUST)) +
		geom_alluvium(aes(fill=GRP)) +
		geom_stratum(aes()) +
		geom_text(stat = "stratum") +
		theme_minimal()

	ggplot(data = DATA_long8,
		   aes(x = timepoint, stratum = stratum, alluvium = alluvium,
			   y = Freq, label = CLUST)) +
		geom_alluvium(aes(fill=GRP)) +
		geom_stratum(alpha = .5) +
		geom_stratum() + geom_text(stat = "stratum") +
		theme_minimal()

	ggplot(data = DATA_long56,
		   aes(x = timepoint, stratum = stratum, alluvium = alluvium,
			   y = Freq, label = CLUST)) +
		geom_alluvium(aes(fill=GRP)) +
		geom_stratum(alpha = .5) +
		geom_stratum() + geom_text(stat = "stratum") +
		theme_minimal()
dev.off()

############################################################

DATAback0.1<-getTransitDatBack(LDATA.back)
DATAback<-getTransitDatBack(LDATA.back, CUTOFF=5)


DATAback_long0 <- to_lodes_form(droplevels(DATAback[DATAback$Freq>0.2,]%>%mutate(GRP=Var1=="KJP5_2")),
                              key = "TP",
                              axes =c(1,2,4,6,8,10))%>%mutate(timepoint=paste0("E", gsub("_[0-9]+|KJ","",stratum)))%>%
				mutate(CLUST=gsub("^.*_","",stratum))

DATAback_long8 <- to_lodes_form(droplevels(DATAback[DATAback$Freq>0.2,]%>%mutate(GRP=Var1=="KJP5_8")),
                              key = "TP",
                              axes =c(1,2,4,6,8,10))%>%mutate(timepoint=paste0("E", gsub("_[0-9]+|KJ","",stratum)))%>%
				mutate(CLUST=gsub("^.*_","",stratum))

DATAback_longSZ <- to_lodes_form(droplevels(DATAback[DATAback$Freq>0.2,]%>%mutate(GRP=Var2=="KJ185_8")),
                              key = "TP",
                              axes =c(1,2,4,6,8,10))%>%mutate(timepoint=paste0("E", gsub("_[0-9]+|KJ","",stratum)))%>%
				mutate(CLUST=gsub("^.*_","",stratum))


DATAback_long4 <- to_lodes_form(droplevels(DATAback[DATAback$Freq>0.2,]%>%mutate(GRP=Var1=="KJP5_4")),
                              key = "TP",
                              axes =c(1,2,4,6,8,10))%>%mutate(timepoint=paste0("E", gsub("_[0-9]+|KJ","",stratum)))%>%
				mutate(CLUST=gsub("^.*_","",stratum))

pdf("all-in-one.kNN.backward.pdf", width=16)
	ggplot(data = DATAback_long0,
		   aes(x = timepoint, stratum = stratum, alluvium = alluvium,
			   y = Freq, label = CLUST)) +
		geom_alluvium(aes(fill=GRP)) +
		geom_stratum(aes()) +
		geom_text(stat = "stratum") +
		theme_minimal()

	ggplot(data = DATAback_long8,
		   aes(x = timepoint, stratum = stratum, alluvium = alluvium,
			   y = Freq, label = CLUST)) +
		geom_alluvium(aes(fill=GRP)) +
		geom_stratum(alpha = .5) +
		geom_stratum() + geom_text(stat = "stratum") +
		theme_minimal()

	ggplot(data = DATAback_longSZ ,
		   aes(x = timepoint, stratum = stratum, alluvium = alluvium,
			   y = Freq, label = CLUST)) +
		geom_alluvium(aes(fill=GRP)) +
		geom_stratum(alpha = .5) +
		geom_stratum() + geom_text(stat = "stratum") +
		theme_minimal()


	ggplot(data = DATAback_long4,
		   aes(x = timepoint, stratum = stratum, alluvium = alluvium,
			   y = Freq, label = CLUST)) +
		geom_alluvium(aes(fill=GRP)) +
		geom_stratum(alpha = .5) +
		geom_stratum() + geom_text(stat = "stratum") +
		theme_minimal()
dev.off()



