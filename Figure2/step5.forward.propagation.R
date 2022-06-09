#step5.reconstruct.tree.noCC.R

load("successiveTP/MATflow.noCC.E12.to.P5.RData")

library(dplyr)
library(tidyverse)

library(ggplot2)
library(ggalluvial)
data(majors)
################################

data1213<-data.frame(table(MATflow1213[,4], MATflow1213[,5]))%>%filter(Freq>0)
data1314<-data.frame(table(MATflow1314[,4], MATflow1314[,5]))%>%filter(Freq>0)
data1415<-data.frame(table(MATflow1415[,4], MATflow1415[,5]))%>%filter(Freq>0)
data1518<-data.frame(table(MATflow1518[,4], MATflow1518[,5]))%>%filter(Freq>0)
data18P5<-data.frame(table(MATflow18P5[,4], MATflow18P5[,5]))%>%filter(Freq>0)


LDATA<-list(data1213, data1314, data1415, data1518, data18P5)

DATA<-data1213
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
	print(table(DATA$Freq>0.1))
	DATA<-DATA[DATA$Freq>0.1, ]
}

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


pdf("all-in-one.pdf", width=16)
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


