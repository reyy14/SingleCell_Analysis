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
d <- duplicated(fpkm$"Gene")
fpkm<-fpkm[!d,]
rownames(fpkm) <- fpkm$"Gene"
fpkm <- fpkm[,2:ncol(fpkm)]

write.table(fpkm,"Exp.txt",col.names=T, row.names=T, sep="\t", quote=F)

sample <- read.delim("Vdat.xls",row.names=1,header=T)
#gene_ann <- read.delim("../attribute.txt",row.names=1)
#pd <- new("AnnotatedDataFrame",data=sample)
#fd <- new("AnnotatedDataFrame",data=gene_ann)
#HSKJ <- newCellDataSet(as.matrix(fpkm), phenoData=pd, featureData=fd)
#HSKJ <- newCellDataSet(as.matrix(fpkm), phenoData=pd)
#eKRT7 <- log(exprs(HSKJ)["KRT7",]+1,2)+0.1
#Exp <- log(exprs(HSKJ)+1)

Type = c()
NewCol <- c()
for (i in 1:ncol(fpkm)){
	Ins = fpkm["INS",i]
	Gcg = fpkm["GCG",i]
	Sst = fpkm["SST",i]
	Ghrl = fpkm["GHRL",i]
	PPY = fpkm["PPY",i]
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
	}
	if ((Ins>1000) & (Gcg>1000)) type="D"
	if ((Ins<1000) & (type=="B")) type="X"
	if ((Gcg<1000) & (type=="A")) type="X"
	Type = c(Type,type)
	mcol <- strsplit(colnames(fpkm)[i], "_")[[1]]
	mmcol <- paste(type,"_",mcol[1],"_",mcol[2], sep="")
	NewCol <- c(NewCol, mmcol)
}

fpkm2 <- fpkm
colnames(fpkm2)<-NewCol
write.table(fpkm2,"ExpExt.txt",col.names=T, row.names=T, sep="\t", quote=F)
colO <- grep("X",colnames(fpkm2))
fpkm3 <- fpkm2[,c(colO)]
write.table(fpkm3,"FPKM_X.txt",col.names=T, row.names=T, sep="\t", quote=F)

colO <- grep("Old",colnames(fpkm2))
colY <- grep("Young",colnames(fpkm2))
fpkm2 <- fpkm2[,c(colO,colY)]
write.table(fpkm2,"FPKM_OY.txt",col.names=T, row.names=T, sep="\t", quote=F)

sample <- cbind(sample,Type)
pd <- new("AnnotatedDataFrame",data=sample)
HSKJ <- newCellDataSet(as.matrix(fpkm), phenoData=pd)
Exp <- log((exprs(HSKJ)+1))


Old = HSKJ[,pData(HSKJ)[,1]=="Old"] # 168 samples
Young = HSKJ[,pData(HSKJ)[,1]=="Young"] # 144 samples
T2D_Old = HSKJ[,pData(HSKJ)[,1]=="T2DOld"] #96  
T1D = HSKJ[,pData(HSKJ)[,1]=="T1D"] #96  
T2D_OldC = HSKJ[,pData(HSKJ)[,1]=="T2DOldC"] #54  

Old_A = Old[,pData(Old)["Type"]=="A"]
Young_A = Young[,pData(Young)["Type"]=="A"]
Old_B = Old[,pData(Old)["Type"]=="B"]
Young_B = Young[,pData(Young)["Type"]=="B"]

#al <- grep("A_", NewCol)
#be <- grep("B_", NewCol)

#fpkm2 <- fpkm
#colnames(fpkm2)<-NewCol

#HSKJ <- newCellDataSet(as.matrix(fpkm2), phenoData=pd)

P <- read.table("P.gene.txt",row.names=1)
for (i in 1:5) {

	if (i==2) {
		#My <- Young
		My <- Young_A
		#sig_mT <- log(exprs(Young[rownames(P),])+1)
		sig_mT <- log(exprs(Young_A[rownames(P),])+1)
		nOldA <- length(grep("Young_E", colnames(My)))
		nOldB <- length(grep("Young_F", colnames(My)))
		title = "Young_A"
	}
	if (i==3) {
		#My <- Old
		My <- Old_A
		sig_mT <- log(exprs(Old_A[rownames(P),])+1)
		nOldA <- length(grep("Old_A", colnames(My)))
		nOld <- length(grep("Old", colnames(My)))
		nOldB <- nOld-nOldA
		title = "Old_A"
	}
	if (i==4) {
		#My <- T2D_Old
		My <- Young_B
		sig_mT <- log(exprs(Young_B[rownames(P),])+1)
		nOldA <- length(grep("Young_E", colnames(My)))
		nOldB <- length(grep("Young_F", colnames(My)))
		title = "Young_B"
	}
	if (i==5) {
		#My <- T1D
		My <- Old_B
		sig_mT <- log(exprs(Old_B[rownames(P),])+1)
		nOldA <- length(grep("Old_A", colnames(My)))
		nOld <- length(grep("Old", colnames(My)))
		nOldB <- nOld-nOldA
		title = "Old_B"
	}
	if (i==6) {
		My <- T2D_OldC
		sig_mT <- log(exprs(T2D_OldC[rownames(P),])+1)
		nOldA <- length(grep("T2DOldC_", colnames(My)))
		title = "T2DOldC"
	}
	if (i!=1) { 
		rCol = c(rep(col2hex("khaki"),nOldA) ,rep(col2hex("magenta"),nOldB))
		#else rCol = rep(col2hex("khaki"),nOldA)
	}

	

	if (i==1) {
		My <- cbind(Young,Old,T2D_Old,T1D,T2D_OldC)
		sig_mT <- log(exprs(HSKJ[rownames(P),])+1)
		title = "All"
		rCol = c(rep(col2hex("khaki"),ncol(Young)), rep(col2hex("magenta"),ncol(Old)),rep(col2hex("gray"),ncol(T2D_Old)), rep(col2hex("coral"),ncol(T2D_OldC)), rep(col2hex("cyan"),ncol(T1D)))
	}

	if (i==1) {
		sig_d <- as.dist(1-cor(t(sig_mT)))
		sig_h <- hclust(sig_d, method="ward.D2")
		sig_dend = as.dendrogram(sig_h)
	}
	#sig_d <- as.dist(1-cor(t(sig_mT)))
	#sig_h <- hclust(sig_d, method="ward.D2")
	#sig_dend = as.dendrogram(sig_h)
	for (j in 1:ncol(sig_mT)){
		sig_mT[1,j] = sig_mT[1,j]+0.001  
	}
	bk <- seq(-3, 3, by=0.1)
	cluster=vector(mode="character",length=nrow(sig_mT))
	rowColor=vector(mode="character",length=nrow(sig_mT))
	remove_gene <- paste("Destiny/Panc/",title,"*",sep="")
	N.cluster=3	
	
	sig_d2 <- as.dist(1-cor((sig_mT)))
	sig_h2 <- hclust(sig_d2, method="ward.D2")
	sig_dend2 = as.dendrogram(sig_h2)
	
	#unlink(remove_gene)
	for (j in 1:N.cluster) {
   		sub.tf=cutree(sig_h,k=N.cluster)==j
		rowColor[sub.tf]=colorL[j]
		clustername = paste("Destiny/",title,".",colorL[j],j,".txt",sep="")
	    print(clustername)
	    write.table(names(sub.tf[sub.tf==TRUE]),clustername,quote=F,col.names=F,row.names=F)
	    cluster[sub.tf]=as.character(j)
	}       
	cluster <- as.matrix(cluster)
	sig_mT <- t(scale(t(sig_mT)))
	sig_mT[sig_mT>3]=3
	sig_mT[sig_mT< -3]= -3
	rownames(cluster)<-rownames(sig_mT)

	C.cluster=10
	if (i==1) {
		clusterC=vector(mode="character",length=ncol(sig_mT))
		for (j in 1:C.cluster) {
   			sub.tf=cutree(sig_h2,k=C.cluster)==j
			rCol[sub.tf]=colorL[j]
			clustername = paste("Destiny/Column",title,".",colorL[j],j,".txt",sep="")
		    print(clustername)
		    write.table(names(sub.tf[sub.tf==TRUE]),clustername,quote=F,col.names=F,row.names=F)
		    clusterC[sub.tf]=as.character(j)
		 
			print(c("Young", 100*length(grep("Young",names(sub.tf[sub.tf==TRUE])))/144))
			print(c("T1D", 100*length(grep("T1D",names(sub.tf[sub.tf==TRUE])))/96))
			print(c("T2DOld_", 100*length(grep("T2DOld_",names(sub.tf[sub.tf==TRUE])))/190))
			print(c("T2DOldC_", 100*length(grep("T2DOldC",names(sub.tf[sub.tf==TRUE])))/144))
			print(c("Old_", 100*(length(grep("Old_",names(sub.tf[sub.tf==TRUE])))-length(grep("T2DOld_",names(sub.tf[sub.tf==TRUE])))) /168))
		
			mylen = length(names(sub.tf[sub.tf==TRUE])) # length(names(sub.tf))
			print(c("| Young", 100*length(grep("Young",names(sub.tf[sub.tf==TRUE])))/mylen))
			print(c("| T1D", 100*length(grep("T1D",names(sub.tf[sub.tf==TRUE])))/mylen))
			print(c("| T2DOld_", 100*length(grep("T2DOld_",names(sub.tf[sub.tf==TRUE])))/mylen))
			print(c("| T2DOldC_", 100*length(grep("T2DOldC",names(sub.tf[sub.tf==TRUE])))/mylen))
			print(c("| Old_", 100*(length(grep("Old_",names(sub.tf[sub.tf==TRUE])))-length(grep("T2DOld_",names(sub.tf[sub.tf==TRUE])))) /mylen))
		}

		clusterC <- as.matrix(clusterC)
	}
	cluster <- as.matrix(cluster)
	sig_mT <- t(scale(t(sig_mT)))
	sig_mT[sig_mT>3]=3
	sig_mT[sig_mT< -3]= -3
	rownames(cluster)<-rownames(sig_mT)

	
	#unlink("Destiny/Panc/heatm*")
	pngname = paste("Destiny/Panc",title,"png",sep=".")
	png(pngname,width=1000,height=800)
	if (i==1) {
		heatmap.2(as.matrix(sig_mT), Rowv=sig_dend, Colv=sig_dend2, scale="row",dendrogram="both",trace="none",labCol=NULL,
	col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=1, ColSideColors=rCol,  RowSideColors=rowColor, #hline=max(breaks),
	 cexCol=0.1,cexRow=2, margins=c(0,20))#
	}
	else {
		heatmap.2(as.matrix(sig_mT), Rowv=sig_dend, Colv=sig_dend2, scale="row",dendrogram="both",trace="none",labCol=NULL,
	col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=1, ColSideColors=rCol,  RowSideColors=rowColor, #hline=max(breaks),
	 cexCol=0.1,cexRow=2, margins=c(0,20))#colsep=1:ncol(my), rowsep=1:nrow(my), sepwidth=c(0,5,0.5), sepcolor="black")# trace=T,  hline=seq(0,1,1/(ncol(my))) )
	#mylegend=c(labA, labB)
	 #legend("top",mylegend,cex=2.5,col=c("khaki","magenta"), pch=c(15,15), bty="n",horiz=TRUE)
	 }
	 dev.off()
}
