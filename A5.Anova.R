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
d2 <- duplicated(fpkm[1:ncol(fpkm)])
fpkm<-fpkm[!d2,]

sample <- read.delim("Vdat.xls",row.names=1,header=T)
pd <- new("AnnotatedDataFrame",data=sample)
HSKJ <- newCellDataSet(as.matrix(fpkm), phenoData=pd)
Exp <- log(exprs(HSKJ)+1)

Old = HSKJ[,pData(HSKJ)[,1]=="Old"] # 168 samples
Young = HSKJ[,pData(HSKJ)[,1]=="Young"] # 144 samples
T2D_Old = HSKJ[,pData(HSKJ)[,1]=="T2DOld"] #96  
T1D = HSKJ[,pData(HSKJ)[,1]=="T1D"] #96  
T2D_OldC = HSKJ[,pData(HSKJ)[,1]=="T2DOldC"] #54  

Myexp = Exp[,c(colnames(Old),colnames(Young),colnames(T2D_Old),colnames(T1D))]
clu <- c(rep("2",ncol(Old)), rep("3",ncol(Young)),rep("5",ncol(T2D_Old)),rep("6",ncol(T1D)))
mT=c()
for (i in 1:nrow(Old)) {
	myexp <- Myexp[i,]
	myframe <- data.frame(myexp, clu)
	gene <- rownames(Exp)[i]
	if (sum(myexp)==0) next
	if (length( grep("OTT",gene) )==0){
		ano_results =aov(myexp~clu,data=myframe)
		pval <- summary(ano_results)[[1]][["Pr(>F)"]][[1]]
		if (pval<0.0001){
			mT = rbind(mT,c(gene,pval))
		}
	}
}
write.table(mT,"anova.txt",sep="\t",quote=F, row.names=F,col.names=F)

for (i in 1:1){
	N.cluster = 10 
	#Ano = read.table("annova.txt",row.names=1, header=F)
	sig_mT <- Myexp[mT[,1],]
	sig_d <- as.dist(1-cor(t(sig_mT)))
	sig_h <- hclust(sig_d, method="ward.D2")
	sig_dend = as.dendrogram(sig_h)
	cluster=vector(mode="character",length=nrow(sig_mT))
	rowColor=vector(mode="character",length=nrow(sig_mT))
	#rCol = vector(mode="character",length=ncol(sig_mT))
	#rCol = c(rep(col2hex("khaki"),ncol(A)), rep(col2hex("magenta"),ncol(B)))
	
#colorL=c("red","purple","blue","yellow","green","orange","brown","gray","black","coral","beige","cyan","pink","khaki","magenta")
	rCol<- c(rep(col2hex("brown"),ncol(Old)), rep(col2hex("gray"),ncol(Young)),rep(col2hex("black"),ncol(T2D_Old)),rep(col2hex("pink"),ncol(T1D)))
	#remove_gene <- paste("Destiny/Results/",title,"*",sep="")
	for (j in 1:N.cluster) {
   		sub.tf=cutree(sig_h,k=N.cluster)==j
		rowColor[sub.tf]=colorL[j]
		clustername = paste("Destiny/Results/Anova.",colorL[j],j,".txt",sep="")
	    print(clustername)
	    write.table(names(sub.tf[sub.tf==TRUE]),clustername,quote=F,col.names=F,row.names=F)
	    cluster[sub.tf]=as.character(j)
	}       
	cluster <- as.matrix(cluster)
	sig_mT <- t(scale(t(sig_mT)))
	bk <- seq(-2, 2, by=0.1)
	sig_mT[sig_mT>2]=2
	sig_mT[sig_mT< -2]= -2
	rownames(cluster)<-rownames(sig_mT)
	pngname = paste("Destiny/Anova.png",sep=".")
	png(pngname,width=1000,height=800)
	heatmap.2(as.matrix(sig_mT), Rowv=sig_dend, Colv=F, scale="none",dendrogram="row",trace="none",labCol=NULL,
	col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=1, ColSideColors=rCol,  RowSideColors=rowColor, #hline=max(breaks),
	 cexCol=0.1,cexRow=0.01,margins=c(0,20))#colsep=1:ncol(my), rowsep=1:nrow(my), sepwidth=c(0,5,0.5), sepcolor="black")# trace=T,  hline=seq(0,1,1/(ncol(my))) )
	mylegend=c("Old","Young","T2D","T1D")
	 legend("top",mylegend,cex=2.5,col=c("brown","gray","black","pink"), pch=c(15,15,15,15), bty="n",horiz=TRUE)
	 dev.off()

}


