args <- commandArgs(TRUE)

data<-read.table(args[1])
name<-args[2]

agg<-aggregate(data$V8,by=list(data$V4,gsub("chr","", data$V1),data$V5), FUN=t)
agg<-as.data.frame(as.matrix(agg))
rownames(agg)<-agg[,1]
agg[,1]<-NULL

cols<-ifelse(unique(data$V6)/1000000>=1, paste(unique(data$V6)/1000000, "M", sep=""), unique(data$V6))
colnames(agg)<-c("Chr","Position",cols)
agg<-agg[order(agg$Chr, agg$Position),]
write.table(agg,paste(name,".txt",sep=""),quote=F, row.names=T, sep="\t")