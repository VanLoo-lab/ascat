args <- commandArgs(TRUE)

library(plyr)

data<-read.table(args[1])
name<-args[2]

agg<-ddply(data, .(V4, V1, V5), function(z) return(t(round((z$V10+z$V11)/(z$V9+z$V10+z$V11+z$V12),6))))
cols<-unique(data$V6)
cols2=cols
cols2[which(cols<1e3)]=paste0(cols[which(cols<1e3)],'bp')
cols2[which(cols>=1e3 & cols<1e6)]=paste0(cols[which(cols>=1e3 & cols<1e6)]/1e3,'kb')
cols2[which(cols>=1e6)]=paste0(cols[which(cols>=1e6)]/1e6,'M')

rownames(agg)<-agg[,1]
agg[,1]<-NULL

colnames(agg)<-c("Chr","Position",cols2)
agg<-agg[order(agg$Chr, as.numeric(as.character(agg$Position))),]
write.table(agg,paste(name,".txt",sep=""),quote=F, row.names=T, sep="\t")
