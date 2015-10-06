data<-read.table()

agg<-aggregate(data$V8,by=list(data$V4,gsub("chr","", data$V1),data$V5), FUN=t)
agg<-as.data.frame(as.matrix(agg))
colnames(agg)<-c("Name","Chr","Position",unique(data$V6))
agg<-agg[order(agg$Chr, agg$Position),]
write.table(agg,"GC_Illumina_Human_CoreExome_24_v1.0.txt",row.names=F, quote=F, sep="\t")