# R Script to compare the MTBC0 to H37Rv and Comas et al ancestor and create bedtracks, with calculated 
# nucleotide identity in 100 bp windows

library(seqinr)
library(ape)

### Read input files
comas<-read.fasta("../../other.reference.sequences/Comas_MTB_anc.fasta")
mtbc0v1.1<-read.fasta("../../MTBC0_v1.1.fasta")
h37rv<-read.fasta("../../other.reference.sequences/H37Rv.fasta")
h37rv.mtbc0.maf<-read.table("H37Rv.on.MTBC0v1.1.maf",header=F)

identity.matrix<-rbind(rep(NA,times=length(mtbc0v1.1[[1]])),rep(NA,times=length(mtbc0v1.1[[1]])))
reconstructed.comas<-comas
reconstructed.h37rv<-h37rv
reconstructed.comas$MTB_anc<-mtbc0v1.1$MTBC0
attr(reconstructed.comas$MTB_anc,"name")<-attr(comas$MTB_anc,"name")
attr(reconstructed.comas$MTB_anc,"Annot")<-attr(comas$MTB_anc,"Annot")
reconstructed.h37rv$NC_000962.3<-mtbc0v1.1$MTBC0
attr(reconstructed.h37rv$NC_000962.3,"name")<-attr(h37rv$NC_000962.3,"name")
attr(reconstructed.h37rv$NC_000962.3,"Annot")<-attr(h37rv$NC_000962.3,"Annot")
rownames(identity.matrix) <- c("Comas","H37Rv")

aligned.length<-0
for(i in which(h37rv.mtbc0.maf[,1] %in% "Anc01")) {
 mtbc0.seq<-strsplit(h37rv.mtbc0.maf[i,4],split="")[[1]]
 h37rv.seq<-strsplit(h37rv.mtbc0.maf[i+1,4],split="")[[1]]
 h37rv.pos<-h37rv.mtbc0.maf[i+1,2]
 for(j in 1:length(mtbc0.seq)) {
  mtbc0.pos<-h37rv.mtbc0.maf[i,2]+j
  if(h37rv.seq[j] != "-") {
   h37rv.pos<-h37rv.pos+1
   if(toupper(mtbc0.seq[j]) == toupper(comas[[1]][h37rv.pos])) identity.matrix[1,mtbc0.pos]<-1 else  identity.matrix[1,mtbc0.pos]<-0
   if(toupper(mtbc0.seq[j]) == toupper(h37rv.seq[j])) identity.matrix[2,mtbc0.pos]<-1 else  identity.matrix[2,mtbc0.pos]<-0
  } else {
   ## GAP
   identity.matrix[1:2,mtbc0.pos]<-"-"
  }
 }
aligned.length<-aligned.length+length(mtbc0.seq)
}
comas.identity<-data.frame(chr="MTBC0",start=0:(dim(identity.matrix)[2]-1),end=1:dim(identity.matrix)[2],identity=identity.matrix[1,])
window.size<-100
comas.window.identity<-data.frame(chr="MTBC0",start=seq(0,(dim(identity.matrix)[2]-1),by=window.size),end=(seq(0,(dim(identity.matrix)[2]),by=window.size))+window.size-1,identity=NA)
pb<-txtProgressBar(min=0,max=nrow(comas.window.identity),style=3)
for(i in 1:nrow(comas.window.identity)) {
 t.values<-comas.identity[seq(comas.window.identity[i,2],comas.window.identity[i,3],by=1),4]
 t.sum<-mean(as.integer(t.values),na.rm=T)
 comas.window.identity[i,"identity"]<-t.sum
 setTxtProgressBar(pb,i)
}
close(pb)
comas.window.identity<-comas.window.identity[-which(is.na(comas.window.identity[,"identity"])),]
comas.identity<-comas.identity[-which(comas.identity[,4] %in% "-"),]

## merge overlapping regions
comas.merged<-data.frame(chr = character(0), start = numeric(0), end = numeric(0), identity = numeric(0))
current_region <- NULL
pb<-txtProgressBar(min=0,max=nrow(comas.identity),style=3)
for(i in 1:nrow(comas.identity)) {
 if (is.null(current_region)) {
  current_region <- comas.identity[i, ]
 } else {
  if(comas.identity[i, "chr"] == current_region$chr && comas.identity[i, "identity"] == current_region$identity && comas.identity[i, "start"] <= current_region$end) {
   current_region$end <- max(current_region$end, comas.identity[i, "end"])
  } else {
   comas.merged <- rbind(comas.merged, current_region)
   current_region <- comas.identity[i, ]
  }
 }
 setTxtProgressBar(pb,i)
}
close(pb)
if(!is.null(current_region)) {
  comas.merged <- rbind(comas.merged, current_region)
}
options(scipen = 999)
## Write out the bedtrack file
write("track type=bedGraph",file="Comas.vs.mtbc0.identity.windowed.bedgraph")
write.table(comas.window.identity,file="Comas.vs.mtbc0.identity.windowed.bedgraph",col.names=F,row.names=F,quote=F,append=TRUE,sep="\t")



h37rv.identity<-data.frame(chr="MTBC0",start=0:(dim(identity.matrix)[2]-1),end=1:dim(identity.matrix)[2],identity=identity.matrix[2,])
window.size<-100
h37rv.window.identity<-data.frame(chr="MTBC0",start=seq(0,(dim(identity.matrix)[2]-1),by=window.size),end=(seq(0,(dim(identity.matrix)[2]),by=window.size))+window.size-1,identity=NA)
pb<-txtProgressBar(min=0,max=nrow(h37rv.window.identity),style=3)
for(i in 1:nrow(h37rv.window.identity)) {
 #t.pos<-which(h37rv.identity[,2] %in% seq(h37rv.window.identity[i,2],h37rv.window.identity[i,3],by=1))
 t.values<-h37rv.identity[seq(h37rv.window.identity[i,2],h37rv.window.identity[i,3],by=1),4]
 t.sum<-mean(as.integer(t.values),na.rm=T)
 h37rv.window.identity[i,"identity"]<-t.sum
 setTxtProgressBar(pb,i)
}
close(pb)
h37rv.window.identity<-h37rv.window.identity[-which(is.na(h37rv.window.identity[,"identity"])),]
h37rv.identity<-h37rv.identity[-which(h37rv.identity[,4] %in% "-"),]

## merge overlapping regions
h37rv.merged<-data.frame(chr = character(0), start = numeric(0), end = numeric(0), identity = numeric(0))
current_region <- NULL
pb<-txtProgressBar(min=0,max=nrow(h37rv.identity),style=3)
for(i in 1:nrow(h37rv.identity)) {
 if (is.null(current_region)) {
  current_region <- h37rv.identity[i, ]
 } else {
  if(h37rv.identity[i, "chr"] == current_region$chr && h37rv.identity[i, "identity"] == current_region$identity && h37rv.identity[i, "start"] <= current_region$end) {
   current_region$end <- max(current_region$end, h37rv.identity[i, "end"])
  } else {
   h37rv.merged <- rbind(h37rv.merged, current_region)
   current_region <- h37rv.identity[i, ]
  }
 }
 setTxtProgressBar(pb,i)
}
close(pb)
if(!is.null(current_region)) {
  h37rv.merged <- rbind(h37rv.merged, current_region)
}

## Write out the bedtrack file
write("track type=bedGraph",file="H37Rv.vs.mtbc0.identity.windowed.bedgraph")
write.table(h37rv.window.identity,file="H37Rv.vs.mtbc0.identity.windowed.bedgraph",col.names=F,row.names=F,quote=F,append=TRUE,sep="\t")







