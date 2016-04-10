library(diptest)
library(monocle)
library(destiny)
library(RColorBrewer)
library(gplots)

colorL=c("red","purple","blue","yellow","green","orange","brown","gray","black","coral","beige","cyan","pink","khaki","magenta")
set.seed(1)
fpkm <- read.table("FPKM_OY.txt",header=T, row.names=1)

youngA <- fpkm[,grep("YoungA",colnames(fpkm))]
oldA <- fpkm[,grep("OldA",colnames(fpkm))]

youngB <- fpkm[,grep("YoungB",colnames(fpkm))]
oldB <- fpkm[,grep("OldB",colnames(fpkm))]

mT = c()
for (i in 1:nrow(fpkm)) {
	gene <- rownames(fpkm)[i]
	if ((sum(youngA[i,])==0) && (sum(oldA[i,])==0)) next
	if (length(grep("OTT",gene))==0){
		tA = t.test(youngA[i,],oldA[i,])
		if(tA$p.value<0.01){
			if(tA$estimate[2]!=0){
				FC = log(tA$estimate[1] / tA$estimate[2],2)
				mT = rbind(mT,c(gene,tA$p.value,FC,tA$estimate[1],tA$estimate[2]))
			}
			else mT = rbind(mT,c(gene,tA$p.value,100,tA$estimate[1],tA$estimate[2]))
		}
	}
}

fname = "Diff.YO_A.dat"
write.table(mT,fname,sep="\t",quote=F,row.names=F,col.names=F)

mT = c()

for (i in 1:nrow(fpkm)) {
	gene <- rownames(fpkm)[i]
	if ((sum(youngB[i,])==0) && (sum(oldB[i,])==0)) next
	if (length(grep("OTT",gene))==0){
		tB = t.test(youngB[i,],oldB[i,])
		if(tB$p.value<0.01){
			if(tA$estimate[2]!=0){
				FC = log(tB$estimate[1] / tB$estimate[2],2)
				mT = rbind(mT,c(gene,tB$p.value,FC,tB$estimate[1],tB$estimate[2]))
			}
			else mT = rbind(mT,c(gene,tB$p.value,100,tB$estimate[1],tB$estimate[2]))
		}
	}
}

fname = "Diff.YO_B.dat"
write.table(mT,fname,sep="\t",quote=F,row.names=F,col.names=F)

stop("bawk")
