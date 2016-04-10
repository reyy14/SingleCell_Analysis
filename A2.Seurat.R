library("Seurat")
#nbt.data = read.table("Exp.txt", sep="\t", header=T, row.names=1)
nbt.data = read.table("FPKM_OY.txt", sep="\t", header=T, row.names=1)
nbt.data=log(nbt.data+1,2)
#corner(nbt.data)
nbt=new("seurat",raw.data=nbt.data)
slotNames(nbt)
#nbt=setup(nbt,project="NBT",min.cells = 4,names.field = 1,names.delim = "_",min.genes = 1000,is.expr=1,)
nbt=setup(nbt,project="NBT",min.cells = 1,names.field = 1,names.delim = "_",min.genes = 1,is.expr=1,)
png("Seurat/ViolinPlot.Cell.png",width=1500, height=800)
vlnPlot(nbt,c("GCG","MIR1244-1","INS","SST","RPL7","PRSS1","FOXP1","LOX","CD99"))
#vlnPlot(nbt,c("MIR1244-1","INS","SST","TRY6","LOX","CD99"))
dev.off()
png("Seurat/ViolinPlot2.Cell.png",width=1500, height=800)
vlnPlot(nbt,c("GCG","NPY","SNORD83B","ATP5L","SNORD50A","PPIC","TPM1","CEACAM7"))
dev.off()
png("Seurat/ViolinPlot3.Cell.png",width=1100, height=400)
vlnPlot(nbt,c("INS","GCG","SST","CD44","PRSS1"),size.use=1)
dev.off()
png("Seurat/ViolinPlot5.Cell.png",width=1500, height=800)
vlnPlot(nbt,c("SEC11A","TPD52","NHP2","FIS1","HSPD1","FDPS","TIMM23","FOXP1","ACVR2A"))
dev.off()
#stop("d")
P.gene <- read.table("P.gene.txt")
png("Seurat/ViolinPlot4.Cell.png",width=3200, height=3000)
vlnPlot(nbt,P.gene[,1])
dev.off()

#stop("dd")
#par(mfrow=c(2,2))
#cellPlot(nbt,nbt@cell.names[1],nbt@cell.names[2],do.ident = FALSE)
#cellPlot(nbt,nbt@cell.names[3],nbt@cell.names[4],do.ident = FALSE)
# Plot two genes against each other, can do this in limited groups of cells
png("Seurat/GCGtest.png",width=800,height=800)
par(mar=c(10,10,5,5))
for (i in 1:length(unique(nbt@ident))){
	pA <- t(nbt.data[,nbt@ident==unique(nbt@ident)[i]]["GCG",])
	pB <- t(nbt.data[,nbt@ident==unique(nbt@ident)[i]]["INS",])

	if (i==1) plot(pA,pB,col=i,pch=16,xlab="GCG (log2FPKM) ",ylab="INS (log2FPKM)",xlim=c(0,20),ylim=c(0,20),cex=1.5,cex.lab=2, cex.axis=2)
	else points(pA,pB,col=i,pch=16,cex=1.5, cex.lab=2, cex.axis=2)
}
legend("topright",legend=unique(nbt@ident),col=c(1,2,3,4,5),pch=16,cex=2)
dev.off()
#stop("d")
library("scatterplot3d")
#library("rgl")
png("Seurat/GCGtest3D.png",width=800,height=800)
#for (i in 1:length(unique(nbt@ident))){
#	pA <- t(nbt.data[,nbt@ident==unique(nbt@ident)[i]]["GCG",])
#	pB <- t(nbt.data[,nbt@ident==unique(nbt@ident)[i]]["INS",])
#	pC <- t(nbt.data[,nbt@ident==unique(nbt@ident)[i]]["SST",])
#	if(i==1) scatterplot3d(pA,pB,pC, highlight.3d=F, col.axis="black", pch=16)
#	else points3d(pA,pB,pC, highlight.3d=F, col.axis="black", pch=16)
#}
#dev.off()



#stop("d")
png("Seurat/Disperse.png",width=800, height=800)
nbt=mean.var.plot(nbt,y.cutoff = 2,x.low.cutoff = 2,fxn.x = expMean,fxn.y = logVarDivMean)
dev.off()


png("Seurat/PCA.png",width=800, height=800)
nbt=pca(nbt,do.print=FALSE)
pca.plot(nbt,1,2,pt.size = 2)
dev.off()
png("Seurat/VizPCA1.png",width=1200, height=800)
viz.pca(nbt,1:2,font.size=1.2 )
dev.off()

nbt=run_tsne(nbt,dims.use = 1:11,max_iter=2000)
png("Seurat/tSNE.png",width=800, height=800)
tsne.plot(nbt,pt.size = 1.5)
dev.off()
nbt=DBclust_dimension(nbt,1,2,reduction.use = "tsne",G.use = 5,set.ident = TRUE)
png("Seurat/tSNE.cluster.png",width=800, height=800)
tsne.plot(nbt,pt.size = 1)
dev.off()
png("Seurat/tSNE.treecluster.png",width=800, height=800)
nbt=buildClusterTree(nbt,do.reorder = TRUE,reorder.numeric = TRUE,pcs.use = 1:11)
dev.off()

png("Seurat/tSNE.labelcluster.png",width=800, height=800)
tsne.plot(nbt,do.label=TRUE, label.pt.size=1, label.cex.text = 2)
dev.off()

save(nbt,file="seurat_files/A2.Robj")

stop("bawk")

my.data=fetch.data(nbt,c("ident","PC1","nGene","orig.ident"))


nbt = set.all.ident(nbt,"orig.ident")
markers.all=find_all_markers(nbt,thresh.test = 3,test.use = "roc", do.print = TRUE)
#markers.use=subset(markers.all,avg_diff>0&power>0.8)$gene
markers.use=markers.all$gene
png("Seurat/Heatmap.png",width=800, height=2000)
doHeatMap(nbt,genes.use = markers.use,slim.col.label = TRUE,remove.key = TRUE,cexRow=1)
dev.off()



nbt = set.all.ident(nbt,"ident")
markers.all=find_all_markers(nbt,thresh.test = 3,test.use = "roc", do.print = TRUE)
markers.use=subset(markers.all,avg_diff>0&power>0.8)$gene
#markers.use=markers.all$gene
png("Seurat/HeatmapCluster.png",width=800, height=2000)
doHeatMap(nbt,genes.use = markers.use,slim.col.label = TRUE,remove.key = TRUE,cexRow=1)
dev.off()

write.table(markers.all[markers.use,],"Seurat/Cluster_Gene.xls",row.names=T, col.names=T, sep="\t",quote=FALSE)

stop("=========================")
my.data=fetch.data(nbt,c("ident","PC1","nGene","orig.ident"))
P.cells=which.cells(nbt,"ESC",id = "orig.ident")
my.data=fetch.data(nbt,"ident",cells.use = ips.cells)
head(my.data,5)
ips.markers=find.markers(nbt,"ESC",thresh.use = 2)
ips.markers=find.markers(nbt,"ESC",thresh.use = 2,test.use = "roc")
png("Seurat/testExp.png",width=800, height=800)
vlnPlot(nbt,c("HPGD","FN1","NANOG","TDGF1"))
dev.off()

genes.viz = c("POU5F1","NANOG","KRT7","XIST")
png("Seurat/GeneFeature.png",width=800, height=800)
feature.plot(nbt,genes.viz,pt.size = 1)
dev.off()



stop("D")

#write.table(A,"FPKM.txt",row.names=T, col.names=T, sep="\t",quote=FALSE)
library(monocle)
library(destiny)
library(RColorBrewer)
library(gplots)
colorL=c("red","purple","blue","yellow","green","orange","brown","gray","black","coral","beige","cyan","pink","khaki","magenta")
set.seed(1)
fpkm <- read.table("../Vdat.sel.txt",header=T, row.names=1)
sample <- read.delim("../pDat.txt",row.names=1)
gene_ann <- read.delim("../attribute.txt",row.names=1)
pd <- new("AnnotatedDataFrame",data=sample)
fd <- new("AnnotatedDataFrame",data=gene_ann)
HSKJ <- newCellDataSet(as.matrix(fpkm), phenoData=pd, featureData=fd)
eSOX2 <- log(exprs(HSKJ)["SOX2",]+1,2)+0.1
eOCT4 <- log(exprs(HSKJ)["POU5F1",]+1,2)+0.1
eNANOG <- log(exprs(HSKJ)["NANOG",]+1,2)+0.1
eKRT7 <- log(exprs(HSKJ)["KRT7",]+1,2)+0.1
Exp <- exprs(HSKJ)

dm <- DiffusionMap(HSKJ)
png("De1.png",height=1000, width=1000)
plot(dm,pch=20,col.by="Diff")
dev.off()

housekeepers <- c("ATF1","ATF2")
normalizations <- colMeans(exprs(HSKJ)[housekeepers,])
normalized.ct = HSKJ 
exprs(normalized.ct) <- exprs(HSKJ)-normalizations

dm <- DiffusionMap(normalized.ct)
png("De2.png",height=1000, width=1000)
plot(dm,pch=20,col.by="Diff")
dev.off()

png("PC12.png",height=1000, width=1000)
plot(dm,1:2, cex=eSOX2,pch=20,col.by="Diff")
dev.off()

png("PC12.KRT7.png",height=1000, width=1000)
plot(dm,1:2, cex=eKRT7,pch=20,col.by="Diff")
dev.off()


png("PC23.png",height=1000, width=1000)
plot(dm,2:3, pch=20,col.by="Diff")
dev.off()

png("PC13.png",height=1000, width=1000)
plot(dm,c(1,3), pch=20,col.by="Diff")
dev.off()

my <- slot(dm,"eigenvectors")
rownames(my) <- rownames(pData(normalized.ct))
hc <- hclust(dist(cbind(my[,1],my[,2])), method='ward.D2')
png("dendrogram.png",width=1500, height=500)
plot(hc, main="dendro", cex.main=2)
N.cluster = 6
rect.hclust(hc, k=N.cluster, border='red')
dev.off()

cluster = vector(mode="character", length=ncol(normalized.ct))
for (j in 1:N.cluster) {
    sub.tf=cutree(hc,k=N.cluster)==j
    clustername = paste("Results/Cluster.",j,".txt",sep="")
    print(clustername)
    write.table(names(sub.tf[sub.tf==TRUE]),clustername,quote=F,col.names=F,row.names=F)
    cluster[sub.tf]=as.character(j)
}       
pData(normalized.ct)$Cluster <- as.factor(cluster)
# slotNames(dm)
# slot(dm,"d")

dm <- DiffusionMap(normalized.ct)
png("cluster.png",height=1000, width=1000)
plot(dm,pch=20,col.by="Cluster")
dev.off()
png("cluster.P12.png",height=700, width=700)
plot(dm,1:2,pch=20, cex=eSOX2/eSOX2*3, col.by="Cluster",xlab ="PC1", ylab="PC2", cex.lab=1.5, legend.main="Cluster")
dev.off()
png("cluster.Oct4.png",height=1000, width=1000)
plot(dm,1:2,pch=20, cex=eOCT4, col.by="Cluster")
dev.off()
png("cluster.Sox2.png",height=1000, width=1000)
plot(dm,1:2,pch=20, cex=eSOX2, col.by="Cluster")
dev.off()
png("cluster.Nanog.png",height=1000, width=1000)
plot(dm,1:2,pch=20, cex=eNANOG, col.by="Cluster")
dev.off()
png("cluster.Krt7.png",height=1000, width=1000)
plot(dm,1:2,pch=20, cex=eKRT7, col.by="Cluster")
dev.off()

png("PC12.update.png",height=700, width=700 )
plot(dm,1:2, cex=eSOX2/eSOX2*3,pch=20,col.by="Diff",xlab ="PC1", ylab="PC2", cex.lab=1.5, legend.main="Cell" )
dev.off()

#C1 <- read.table("Results/Cluster.1.txt",header=F)
C2 <- read.table("Results/Cluster.2.txt",header=F)
C3 <- read.table("Results/Cluster.3.txt",header=F)
#C4 <- read.table("Results/Cluster.4.txt",header=F)
C5 <- read.table("Results/Cluster.5.txt",header=F)
C6 <- read.table("Results/Cluster.6.txt",header=F)
CAll <- c(array(C2[,1]),array(C3[,1]),array(C5[,1]),array(C6[,1]))
CB <- c(array(C2[,1]),array(C5[,1]),array(C6[,1]))
CA <- c(array(C3[,1]))
clu <- c(rep("2",nrow(C2)), rep("3",nrow(C3)),rep("5",nrow(C5)),rep("6",nrow(C6)))
mT=c()
mX=c()
mESC=c()
for (i in 1:nrow(Exp)) {
	myexp <- Exp[i,CAll]
	myframe <- data.frame(myexp, clu)
	gene <- rownames(Exp)[i]
	if (sum(myexp)==0) next
	if (length( grep("OTT",gene) )==0){
		ano_results =aov(myexp~clu,data=myframe)
		pval <- summary(ano_results)[[1]][["Pr(>F)"]][[1]]
		if (pval<0.001){
			mT = rbind(mT,c(gene,pval))
		}
	
		tA = t.test(Exp[i,CA], Exp[i,CB])
		tESC = t.test(Exp[i,CA][1:19], Exp[i,array(C2[,1])])
		
		if (tA$estimate[2]!=0) {
			FC = tA$estimate[1] / tA$estimate[2]
			if(tA$p.value<0.01)	
					mX = rbind(mX,c(gene,pval,FC,tA$estimate[1],tA$estimate[2]))
		}
		if (tESC$estimate[1]!=0) {  # hibrid
			FC = tESC$estimate[2] / tESC$estimate[1]
			if(tESC$p.value<0.01)	
					mESC = rbind(mESC,c(gene,pval,FC,tESC$estimate[1],tESC$estimate[2]))
		}

	}
}

write.table(mT,"anova.txt",sep="\t",quote=F, row.names=F,col.names=F)
write.table(mX,"ttest.txt",sep="\t",quote=F, row.names=F,col.names=F)
write.table(mESC,"ttest_ESC.txt",sep="\t",quote=F, row.names=F,col.names=F)
#readT <- read.table("anova.txt",header=F, row.names=1)
#readT <- readT[readT[,1]<0.00001, ]
readT <- read.table("ttest.txt",header=F, row.names=1)
readT <- readT[readT[,2]>1.5, ]
sig_mT <- Exp[rownames(readT),CAll]
sig_d <- as.dist(1-cor(t(sig_mT)))
sig_h <- hclust(sig_d, method="ward.D2")
sig_dend = as.dendrogram(sig_h)
bk <- seq(-2, 2, by=0.1)
cluster=vector(mode="character",length=nrow(sig_mT))
rowColor=vector(mode="character",length=nrow(sig_mT))

#sig_mT[sig_mT>2] = 2
#sig_mT[sig_mT< -2] = -2
N2.cluster = 6
rCol <- c(rep(col2hex("red"),nrow(C2)), rep(col2hex("green") ,nrow(C3)),rep(col2hex("cyan"),nrow(C5)),rep(col2hex("violet"),nrow(C6)))
for(j in 1:N2.cluster){
	sub.tf=cutree(sig_h,k=N2.cluster)==j
	rowColor[sub.tf]=colorL[j]
	clustername = paste("Results/C2.",colorL[j],".",j,".txt",sep="")
	write.table(names(sub.tf[sub.tf==TRUE]),clustername,quote=F,col.names=F,row.names=F)
	cluster[sub.tf]=j
}
cluster <- as.matrix(cluster)
rownames(cluster)<-rownames(mT)
png("heatmap.png",width=1000,height=800)
heatmap.2(as.matrix(sig_mT), Rowv=sig_dend, Colv=F, scale="row",dendrogram="row",labRow=rownames(sig_mT),trace="none",labCol=NULL,
	col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=1, ColSideColors=rCol, # RowSideColors=rowColor, ColSideColors=rCol, #hline=max(breaks),
	 cexCol=0.1,cexRow=2.5,margins=c(0,20))#colsep=1:ncol(my), rowsep=1:nrow(my), sepwidth=c(0,5,0.5), sepcolor="black")# trace=T,  hline=seq(0,1,1/(ncol(my))) )
	 #mylegend=c(sprintf("Min: %1.1f", min(pData(BMPX)$nX)) ,sprintf("Max: %1.1f", max(pData(BMPX)$nX))  )
	# legend("top",mylegend,cex=2.5,col=c("black","black"), pch=c(0,15), bty="n",horiz=TRUE)
 dev.off()

readT <- read.table("ttest_ESC.txt",header=F, row.names=1)
readT <- readT[readT[,2]<1, ]
sig_mT <- Exp[c(rownames(readT),"POU5F1"),CAll]
sig_d <- as.dist(1-cor(t(sig_mT)))
sig_h <- hclust(sig_d, method="ward.D2")
sig_dend = as.dendrogram(sig_h)
bk <- seq(-2, 2, by=0.1)
cluster=vector(mode="character",length=nrow(sig_mT))
rowColor=vector(mode="character",length=nrow(sig_mT))

#sig_mT[sig_mT>2] = 2
#sig_mT[sig_mT< -2] = -2
N2.cluster = 6
rCol <- c(rep(col2hex("red"),nrow(C2)), rep(col2hex("green") ,nrow(C3)),rep(col2hex("cyan"),nrow(C5)),rep(col2hex("violet"),nrow(C6)))
for(j in 1:N2.cluster){
	sub.tf=cutree(sig_h,k=N2.cluster)==j
	rowColor[sub.tf]=colorL[j]
	clustername = paste("Results/C2.",colorL[j],".",j,".txt",sep="")
	write.table(names(sub.tf[sub.tf==TRUE]),clustername,quote=F,col.names=F,row.names=F)
	cluster[sub.tf]=j
}
cluster <- as.matrix(cluster)
rownames(cluster)<-rownames(mT)
png("heatmap2.png",width=1000,height=800)
heatmap.2(as.matrix(sig_mT), Rowv=sig_dend, Colv=F, scale="row",dendrogram="row",labRow=rownames(sig_mT),trace="none",labCol=NULL,
	col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=1, ColSideColors=rCol, # RowSideColors=rowColor, ColSideColors=rCol, #hline=max(breaks),
	 cexCol=0.1,cexRow=2.5,margins=c(0,20))#colsep=1:ncol(my), rowsep=1:nrow(my), sepwidth=c(0,5,0.5), sepcolor="black")# trace=T,  hline=seq(0,1,1/(ncol(my))) )
	 #mylegend=c(sprintf("Min: %1.1f", min(pData(BMPX)$nX)) ,sprintf("Max: %1.1f", max(pData(BMPX)$nX))  )
	# legend("top",mylegend,cex=2.5,col=c("black","black"), pch=c(0,15), bty="n",horiz=TRUE)
 dev.off()




