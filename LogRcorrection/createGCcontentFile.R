args <- commandArgs(TRUE)

library(plyr)

data<-read.table(args[1])
name<-args[2]

agg<-ddply(data, .(V4, V1, V5), function(z) return(t(z$V8)))
cols<-ifelse(unique(data$V6)/1000000>=1, paste(unique(data$V6)/1000000, "M", sep=""), unique(data$V6))

rownames(agg)<-agg[,1]
agg[,1]<-NULL

colnames(agg)<-c("Chr","Position",cols)
agg<-agg[order(agg$Chr, as.numeric(as.character(agg$Position))),]
write.table(agg,paste(name,".txt",sep=""),quote=F, row.names=T, sep="\t")
