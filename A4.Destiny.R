#A = read.table("../Vdat.sel.txt",header=T,row.names=1)
#A = t(A)
#write.table(A,"FPKM.txt",row.names=T, col.names=T, sep="\t",quote=FALSE)
library(diptest)
library(monocle)
library(destiny)
library(RColorBrewer)
library(gplots)
colorL=c("red","purple","blue","yellow","green","orange","brown","gray","black","coral","beige","cyan","pink","khaki","magenta")
set.seed(1)
fpkm <- read.table("FPKM.xls",header=T, row.names=1)
o <- order(rowSums(fpkm[,2:ncol(fpkm)]),decreasing=TRUE)
fpkm <- fpkm[o,]
print(dim(fpkm))
d <- duplicated(fpkm$"Gene")
fpkm<-fpkm[!d,]
rownames(fpkm) <- fpkm$"Gene"
print(dim(fpkm))
fpkm <- fpkm[,2:ncol(fpkm)]
#d2 <- duplicated(fpkm[1:ncol(fpkm)])
#fpkm<-fpkm[!d2,]
print(dim(fpkm))
write.table(fpkm,"Exp.txt",col.names=T, row.names=T, sep="\t", quote=F)

Type = c()
NewCol <- c()
for (i in 1:ncol(fpkm)){
	#print(colnames(fpkm)[i])
	Ins = fpkm["INS",i]
	Gcg = fpkm["GCG",i]
	Sst = fpkm["SST",i]
	Ghrl = fpkm["GHRL",i]
	PPY  = fpkm["PPY",i]
	#PRSS1 = fpkm["PRSS1",i]
	hormone = c(Ins,Gcg,Sst,Ghrl,PPY)
	s_hormone = sort(hormone,decreasing=T)
	if (s_hormone[1]<100) type="X"
	else if (s_hormone[1]>20*s_hormone[2]){
		if (s_hormone[1]==Ins) type="B"
		else if (s_hormone[1]==Gcg) type="A"
		else if (s_hormone[1]==Sst) type="S"
		else if (s_hormone[1]==PPY) type="P"
	}
	else { 
		type = "X"
		#print(s_hormone) 
	}
	if ((Ins>1000) & (Gcg>1000)) type="D"
	if ((Ins<1000) & (type=="B")) type="X"
	if ((Gcg<1000) & (type=="A")) type="X"
	Type = c(Type,type)
	mcol <- strsplit(colnames(fpkm)[i], "_")[[1]]
	mmcol <- paste(mcol[1],type,"_",mcol[2], sep="")
	NewCol <- c(NewCol, mmcol)
}
fpkm2 <- fpkm
colnames(fpkm2)<-NewCol
write.table(fpkm2,"ExpExt.txt",col.names=T, row.names=T, sep="\t", quote=F)
colO <- grep("X",colnames(fpkm2))
fpkm3 <- fpkm2[,c(colO)]
write.table(fpkm3,"FPKM_X.txt",col.names=T, row.names=T, sep="\t", quote=F)
#colA <- grep("A_",colnames(fpkm2))
#fpkmA <- fpkm2[,c(colA)]
#colB <- grep("B_",colnames(fpkm2))
#fpkmB <- fpkm2[,c(colB)]
#write.table(fpkmA,"FPKM_A.txt",col.names=T, row.names=T, sep="\t", quote=F)
#write.table(fpkmB,"FPKM_B.txt",col.names=T, row.names=T, sep="\t", quote=F)


colO <- grep("Old",colnames(fpkm2))
colY <- grep("Young",colnames(fpkm2))
fpkm2 <- fpkm2[,c(colO,colY)]
write.table(fpkm2,"FPKM_OY.txt",col.names=T, row.names=T, sep="\t", quote=F)

#stop("D")
sample <- read.delim("Vdat.xls",row.names=1,header=T)
sample <- cbind(sample,Type)

pd <- new("AnnotatedDataFrame",data=sample)
HSKJ <- newCellDataSet(as.matrix(fpkm), phenoData=pd)
Exp <- log((exprs(HSKJ)+1),2)

Old = HSKJ[,pData(HSKJ)[,1]=="Old"] # 168 samples
Young = HSKJ[,pData(HSKJ)[,1]=="Young"] # 144 samples
T2D_Old = HSKJ[,pData(HSKJ)[,1]=="T2D"] #190
T1D = HSKJ[,pData(HSKJ)[,1]=="T1D"] #96  
T2D_OldC = HSKJ[,pData(HSKJ)[,1]=="T2CD"] #54 

Old_A = Old[,pData(Old)["Type"]=="A"]
Young_A = Young[,pData(Young)["Type"]=="A"]

Old_B = Old[,pData(Old)["Type"]=="B"]
Young_B = Young[,pData(Young)["Type"]=="B"]



for (i in 5:5){
	N.cluster = 6
	if (i==1) { M = HSKJ[,(pData(HSKJ)[,1]=="Old") | (pData(HSKJ)[,1]=="Young") ] 
		title = "Old_Young"
		Exp = log(exprs(M)+1,2)
		A = log( exprs( HSKJ[,(pData(HSKJ)[,1]=="Old")])+1,2)
		B = log( exprs( HSKJ[,(pData(HSKJ)[,1]=="Young")])+1,2 )
		labA = "Old"; labB="Young"
		N.cluster =4 
	}
	if (i==2) { M = HSKJ[,(pData(HSKJ)[,1]=="T2D") | (pData(HSKJ)[,1]=="Old") ] 
		title = "T2D"
		Exp = log(exprs(M)+1,2)
		A = log( exprs( HSKJ[,(pData(HSKJ)[,1]=="T2D")])+1,2)
		B = log( exprs( HSKJ[,(pData(HSKJ)[,1]=="Old")])+1,2)
		labA = "T2D Old"; labB = "Old"
		N.cluster = 4 
	}
	if (i==3) { M = HSKJ[,(pData(HSKJ)[,1]=="T1D") | (pData(HSKJ)[,1]=="Old") ] 
		title = "T1D"
		Exp = log(exprs(M)+1,2)
		A = log( exprs( HSKJ[,(pData(HSKJ)[,1]=="T1D")])+1,2)
		B = log( exprs( HSKJ[,(pData(HSKJ)[,1]=="Old")])+1,2)
		labA = "T1D Old"; labB ="Old"
		N.cluster = 3 
	}
	if (i==4) { M = HSKJ[,(pData(HSKJ)[,1]=="T1D") | (pData(HSKJ)[,1]=="T2D") ] 
		title = "T1vs2"
		Exp = log(exprs(M)+1,2)
		A = log( exprs( HSKJ[,(pData(HSKJ)[,1]=="T1D")])+1,2)
		B = log( exprs( HSKJ[,(pData(HSKJ)[,1]=="T2D")])+1,2)
		labA = "T1D"; labB ="T2D"
		N.cluster = 2 
	}
	if (i==5) { 
		title = "YoungOldAlpha"
		Exp = cbind(exprs(Young_A),exprs(Old_A))
		A = log( exprs(Young_A)+1,2)
		B = log( exprs(Old_A)+1,2)
		labA = "Young Alpha"; labB ="Old Alpha"
		N.cluster = 2 
		ExpB = cbind(exprs(Young_B),exprs(Old_B))
		YB = log( exprs(Young_B)+1,2)
		OB = log( exprs(Old_B)+1,2)
	}

	print(i) 


#	dm <- DiffusionMap(M)
#	pngfile = paste("Destiny/PCA",title,"png",sep=".")
#	png(pngfile,height=1000, width=1000)
#	plot(dm,pch=20,col.by="cond")
#	dev.off()
#	housekeepers <- c("GAPDH","ACTB","RPL5","RPL11","RPL14","RPL27","RPL32","RPL34","RPL35","RPS13")
#	overlap = intersect(housekeepers,rownames(exprs(M)))
#	normalizations <- colMeans(exprs(M)[overlap,])
#	normalized.ct = M
#	exprs(normalized.ct) <- exprs(M)-normalizations
#
#	print("  a0")
#	dm <- DiffusionMap(normalized.ct)
#	pngfile = paste("Destiny/PCA.norm",title,"png",sep=".")
#	png(pngfile,height=1000, width=1000)
#	plot(dm,pch=20,col.by="cond")
#	dev.off()
#
#	pngfile = paste("Destiny/PCA.GCG",title,"png",sep=".")
#	png(pngfile,height=1000, width=1000)
#	plot(dm,1:2, pch=20, cex=Exp["GCG",] ,col.by="cond")
#	dev.off()
#	
#	my <- slot(dm,"eigenvectors")
#	rownames(my) <- rownames(pData(normalized.ct))
#	hc <- hclust(dist(cbind(my[,1],my[,2])), method='ward.D2')
#	pngfile = paste("Destiny/Dendro",title,"png",sep=".")
#	png(pngfile,width=1500, height=500)
#	plot(hc, main="dendro", cex.main=2)
#	rect.hclust(hc, k=N.cluster, border='red')
#	dev.off()
#	print("  a")
	# ANOVA
	mT=c()
	for (i in 1:nrow(Exp)) {
		gene <- rownames(Exp)[i]
		if ((sum(A[i,])==0) && (sum(B[i,])==0)) next
		if (length( grep("OTT",gene) )==0){
			tA = t.test(A[i,], B[i,])
			if(tA$p.value<0.01)	{
				if (tA$estimate[2]!=0) {
					FC = log( tA$estimate[1] / tA$estimate[2],2)
					mT = rbind(mT,c(gene,tA$p.value,FC,tA$estimate[1],tA$estimate[2]))
				}	
				else mT = rbind(mT,c(gene,tA$p.value,100,tA$estimate[1],tA$estimate[2]))
			}	
		}

	}



	fname = paste("Destiny/Diff",title,"dat",sep=".")
	write.table(mT,fname,sep="\t",quote=F, row.names=F,col.names=F)
	readT <- read.table(fname,header=F, row.names=1)
	readT <- readT[(readT[,2]>1) | (readT[,2]< -1) , ]
	#sig_mT <- exprs(M[rownames(readT),])
	sig_mT <- cbind(A,B)[rownames(readT),]
	
	sig_d <- as.dist(1-cor(t(sig_mT)))
	sig_h <- hclust(sig_d, method="ward.D2")
	sig_dend = as.dendrogram(sig_h)
	bk <- seq(-3, 3, by=0.1)
	cluster=vector(mode="character",length=nrow(sig_mT))
	rowColor=vector(mode="character",length=nrow(sig_mT))
	#rCol = vector(mode="character",length=ncol(sig_mT))
	rCol = c(rep(col2hex("khaki"),ncol(A)), rep(col2hex("magenta"),ncol(B)))
	remove_gene <- paste("Destiny/Results/",title,"*",sep="")
	unlink(remove_gene)
	for (j in 1:N.cluster) {
   		sub.tf=cutree(sig_h,k=N.cluster)==j
		rowColor[sub.tf]=colorL[j]
		clustername = paste("Destiny/Results/",title,".",colorL[j],j,".txt",sep="")
	    print(clustername)
	    write.table(names(sub.tf[sub.tf==TRUE]),clustername,quote=F,col.names=F,row.names=F)
	    cluster[sub.tf]=as.character(j)
	}       
	print("	b")
	cluster <- as.matrix(cluster)
	sig_mT <- t(scale(t(sig_mT)))
	sig_mT[sig_mT>3]=3
	sig_mT[sig_mT< -3]= -3
	rownames(cluster)<-rownames(sig_mT)
	pngname = paste("Destiny/heatmap",title,"png",sep=".")
	png(pngname,width=1000,height=800)
	heatmap.2(as.matrix(sig_mT), Rowv=sig_dend, Colv=F, scale="none",dendrogram="row",trace="none",labCol=NULL,
	col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=1, ColSideColors=rCol,  RowSideColors=rowColor, #hline=max(breaks),
	 cexCol=0.1,cexRow=0.01,margins=c(0,20))#colsep=1:ncol(my), rowsep=1:nrow(my), sepwidth=c(0,5,0.5), sepcolor="black")# trace=T,  hline=seq(0,1,1/(ncol(my))) )
	mylegend=c(labA, labB)
	 legend("top",mylegend,cex=2.5,col=c("khaki","magenta"), pch=c(15,15), bty="n",horiz=TRUE)
	 dev.off()

	mode = c()
	modeYO1 = c()
	modeYO2 = c()
	modeYO3 = c()
	modeYO4 = c()

	for (i in 1:nrow(Exp)) {
		gene <- rownames(Exp)[i]

		bi_A = dip.test(A[i,])$p.value
		bi_B = dip.test(B[i,])$p.value

		S_bi <- c(bi_A, bi_B)

		if((S_bi[1]<0.001) & (S_bi[2]>0.3)){
			modeYO1 = rbind(modeYO1,c(gene, bi_A, bi_B))
		}
	}

	for (i in 1:nrow(Exp)) {
		gene <- rownames(Exp)[i]

		bi_A = dip.test(A[i,])$p.value
		bi_B = dip.test(B[i,])$p.value

		S_bi <- c(bi_A, bi_B)

		if((S_bi[1]>0.3) & (S_bi[2]<0.001)){
			modeYO2 = rbind(modeYO2,c(gene, bi_A, bi_B))
		}
	}

	for (i in 1:nrow(ExpB)) {
		gene <- rownames(ExpB)[i]

		bi_A = dip.test(YB[i,])$p.value
		bi_B = dip.test(OB[i,])$p.value

		S_bi <- c(bi_A, bi_B)

		if((S_bi[1]<0.001) & (S_bi[2]>0.3)){
			modeYO3 = rbind(modeYO3,c(gene, bi_A, bi_B))
		}
	}

	for (i in 1:nrow(ExpB)) {
		gene <- rownames(ExpB)[i]

		bi_A = dip.test(YB[i,])$p.value
		bi_B = dip.test(YB[i,])$p.value

		S_bi <- c(bi_A, bi_B)

		if((S_bi[1]>0.3) & (S_bi[2]<0.001)){
			modeYO4 = rbind(modeYO4,c(gene, bi_A, bi_B))
		}
	}

	modeTotal = rbind(modeYO1,modeYO2,modeYO3,modeYO4)

	write.table(modeYO1,"Destiny/Bi_YO_A_OX.txt",sep="\t",quote=F,row.names=F,col.names=F)
	write.table(modeYO2,"Destiny/Bi_YO_A_XO.txt",sep="\t",quote=F,row.names=F,col.names=F)
	write.table(modeYO3,"Destiny/Bi_YO_B_OX.txt",sep="\t",quote=F,row.names=F,col.names=F)
	write.table(modeYO4,"Destiny/Bi_YO_B_XO.txt",sep="\t",quote=F,row.names=F,col.names=F)
	write.table(modeTotal,"Destiny/Bi_Total.txt",sep="\t",quote=F,row.names=F,col.names=F)
	
	modA = modeYO1
	modB = modeYO2
	EmodA = A[rownames(modA),]
	EmodB = B[rownames(modB),]
	oA <- order(rowSums(EmodA),decreasing=TRUE)
	oB <- order(rowSums(EmodB),decreasing=TRUE)
	EmodA <- EmodA[oA,]
	EmodB <- EmodB[oB,]
	mybreak = seq(0,20,0.2)
	histA = c()
	histB = c()
	for (i in 1:nrow(modA)) {
		mmA <- hist(A[rownames(modA)[i],],breaks=mybreak)
		histA <- rbind(histA, mmA$density[1:81])
	}

	outname = "Mode.YO_A_OX.txt"
	write.table(EmodA, outname, sep="\t", row.names=T, col.names=F, quote=F)
	pngname = "Mode.YO_A_OX.png"
	par(mfrow=c(1,1),mar=c(5,10,5,2))
	image(t(histA), col=colorpanel(length(bk)-1,"white","blue","blue"), axes=F, xlab="log2 FPKM", cex.lab=2)
	axis(1, at=seq(0,1,1/16) , labels=c(0:16), las=1, tick=FALSE, cex.axis=1.2)
	axis(2, at=seq(0,1,1/(nrow(histA)-1)) , labels=rownames(modA), las=1, tick=TRUE,cex.axis=1.2)
	dev.off()
	stop("bawk")


	#mfile = read.table("Destiny/Bi_Total.txt",row.names=1)

	#modA = mFile[mFile[1]>mFile[2],]
	#modB = mFile[mFile[1]<mFile[2],]
	#EmodA = A[rownames(modA),]
	#EmodB = B[rownames(modB),]
	#oA 

	for (i in 1:nrow(Exp)) {
		gene <- rownames(Exp)[i]
	
		bi_A = dip.test(A[i,])$p.value
		bi_B = dip.test(B[i,])$p.value
	
		S_bi <- sort(c(bi_A, bi_B))
		if((S_bi[1]<0.001) & (S_bi[2]>0.3)) {
			mode = rbind(mode,c(gene, bi_A, bi_B))
		}
	}
	fname = paste("Destiny/Mode",title,"dat",sep=".")
	write.table(mode,fname,sep="\t",quote=F, row.names=F,col.names=F)
	mFile = read.table("Destiny/Bi_Total.txt",row.names=1)

	modA = mFile[ mFile[1]>mFile[2],]
	modB = mFile[ mFile[1]< mFile[2],]
	EmodA = A[rownames(modA),]
	EmodB = B[rownames(modB),]
	oA <- order(rowSums(EmodA),decreasing=TRUE)
	oB <- order(rowSums(EmodB),decreasing=TRUE)
	EmodA <- EmodA[oA,]
	EmodB <- EmodB[oB,]
	mybreak = seq(0,20,0.2)
	histA = c()
	histB = c()
	for (i in 1:nrow(modA)) {
		mmA <- hist(A[rownames(modA)[i],],breaks=mybreak)
		mmB <- hist(B[rownames(modA)[i],],breaks=mybreak)
		histA <- rbind(histA, mmA$density[1:81])
		histB <- rbind(histB, mmB$density[1:81])
	}
	#outname = paste("Destiny/Mode",title,1,"txt",sep=".")
	outname = "EmodA.txt"
	write.table(EmodA,outname,sep="\t",row.names=T,col.names=F,quote=F)
	pngname = paste("Destiny/Mode",title,1,"png",sep=".")
	png(pngname,width=1000,height=1300)
	par(mfrow=c(1,2),mar=c(5,10,5,2))
	image(t(histA), col=colorpanel(length(bk)-1,"white","blue","blue"), axes=F, xlab="log2 FPKM", cex.lab=2)
	#axis(1, at=seq(0,1,1/(ncol(histA)-1)), labels=colnames(tSNP), las=2, tick=FALSE,cex.axis=1.2
	title(labA, cex.main=2)
	axis(1, at=seq(0,1,1/16) , labels=c(0:16), las=1, tick=FALSE, cex.axis=1.2)
	axis(2, at=seq(0,1,1/(nrow(histA)-1)) , labels=rownames(modA), las=1, tick=TRUE,cex.axis=1.2)
	#image(t(histB), col=colorpanel(length(bk)-1,"white","blue","blue"))
	image(t(histB), col=colorpanel(length(bk)-1,"white","blue","blue"), axes=F, xlab="log2 FPKM", cex.lab=2)
	title(labB, cex.main=2)
	axis(1, at=seq(0,1,1/16) , labels=c(0:16), las=1, tick=FALSE, cex.axis=1.2)
	dev.off()
	histA = c()
	histB = c()
	for (i in 1:nrow(modB)) {
		mmA <- hist(A[rownames(modB)[i],],breaks=mybreak)
		mmB <- hist(B[rownames(modB)[i],],breaks=mybreak)
		histA <- rbind(histA, mmA$density[1:81])
		histB <- rbind(histB, mmB$density[1:81])
	}
	outname = paste("Destiny/Mode",title,2,"txt",sep=".")
	write.table(EmodB,outname,sep="\t",row.names=T,col.names=F,quote=F)
	pngname = paste("Destiny/Mode",title,2,"png",sep=".")
	png(pngname,width=1500,height=1300)
	par(mfrow=c(1,2),mar=c(5,10,5,2))
	image(t(histA), col=colorpanel(length(bk)-1,"white","blue","blue"), axes=F, xlab="log2 FPKM", cex.lab=2)
	title(labA, cex.main=2)
	#axis(1, at=seq(0,1,1/(ncol(histA)-1)), labels=colnames(tSNP), las=2, tick=FALSE,cex.axis=1.2
	axis(1, at=seq(0,1,1/16) , labels=c(0:16), las=1, tick=FALSE, cex.axis=1.2)
	axis(2, at=seq(0,1,1/(nrow(histA)-1)) , labels=rownames(modB), las=1, tick=TRUE,cex.axis=1.2)
	#image(t(histB), col=colorpanel(length(bk)-1,"white","blue","blue"))
	image(t(histB), col=colorpanel(length(bk)-1,"white","blue","blue"), axes=F, xlab="log2 FPKM", cex.lab=2)
	title(labB, cex.main=2)
	axis(1, at=seq(0,1,1/16) , labels=c(0:16), las=1, tick=FALSE, cex.axis=1.2)
	dev.off()
	print("  c")

#	nf <- layout(matrix(c(1,2),1,2,byrow=T))
#	heatmap.2(as.matrix(histA), Rowv=NULL, Colv=F, scale="row",dendrogram=NULL,trace="none",labCol=NULL,
#	col=colorpanel(length(bk)-1,"white","white","blue"), key=T,keysize=1,# ColSideColors=rCol,  #RowSideColors=rowColor, #hline=max(breaks),
#	 cexCol=0.1,cexRow=0.01,margins=c(0,20))#colsep=1:ncol(my), rowsep=1:nrow(my), sepwidth=c(0,5,0.5), sepcolor="black")# trace=T,  hline=seq(0,1,1/(ncol(my))) )
#	#mylegend=c(labA, labB)
#	 #legend("top",mylegend,cex=2.5,col=c("khaki","magenta"), pch=c(15,15), bty="n",horiz=TRUE)
#	dev.off()

}


