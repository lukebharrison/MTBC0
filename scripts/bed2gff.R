## MTBC0 helper script
## Small R script to convert a bed file created with haliftover to GFF3 format 
## This is most helpful when the bed file was created from a GFF3 file using gff2bed
## The GFF3 header should be added manually
in.bed<-read.csv("../MTBC0_NC000962.3.gff3.to.bed",sep="\t",header=FALSE)
coord.1<-in.bed[,2]+1
num.records<-dim(in.bed)[1]
out.bed<-rep("MTBC0",times=num.records)
dim(out.bed)<-c(num.records,1)
out.bed<-cbind(out.bed,in.bed[,7],in.bed[,8],coord.1,in.bed[,3],in.bed[,4],in.bed[,6],in.bed[,9],in.bed[,10])
write.table(out.bed,file="MTBC0_NC000962.3.gff3",sep="\t",quote=F,row.names=F,col.names=F)
