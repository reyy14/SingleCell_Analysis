library("Seurat")
#nbt.data = read.table("ExpExt.txt", sep="\t", header=T, row.names=1)
nbt.data = read.table("FPKM_X.txt", sep="\t", header=T, row.names=1)
nbt.data=log(nbt.data+1,2)
#corner(nbt.data)
nbt=new("seurat",raw.data=nbt.data)
slotNames(nbt)
#nbt=setup(nbt,project="NBT",min.cells = 4,names.field = 1,names.delim = "_",min.genes = 1000,is.expr=1,)
nbt=setup(nbt,project="NBT",min.cells = 1,names.field = 1,names.delim = "_",min.genes = 1,is.expr=1,)
png("Seurat/XViolinPlot.Cell.png",width=1500, height=800)
vlnPlot(nbt,c("GCG","MIR1244-1","INS","SST","RPL7","PRSS1","FOXP1","LOX","CD99"))
#vlnPlot(nbt,c("MIR1244-1","INS","SST","TRY6","LOX","CD99"))
dev.off()
png("Seurat/XViolinPlot2.Cell.png",width=1500, height=800)
vlnPlot(nbt,c("GCG","NPY","PRSS1","PRSS2","KLK1","TRY6","ZEB2"))
dev.off()
png("Seurat/XViolinPlot3.Cell.png",width=1500, height=800)
vlnPlot(nbt,c("INS","GCK","ARX","GCG","GHRL","NPY","PPY","SST","CD44"))
dev.off()
P.gene <- read.table("P.gene.txt")


png("Seurat/XDisperse.png",width=800, height=800)
nbt=mean.var.plot(nbt,y.cutoff = 2,x.low.cutoff = 2,fxn.x = expMean,fxn.y = logVarDivMean)
dev.off()


png("Seurat/XPCA.png",width=800, height=800)
nbt=pca(nbt,do.print=FALSE)
pca.plot(nbt,1,2,pt.size = 2)
dev.off()
png("Seurat/VizPCA1.png",width=1200, height=800)
viz.pca(nbt,1:2,font.size=1.2 )
dev.off()

nbt=run_tsne(nbt,dims.use = 1:11,max_iter=2000)
png("Seurat/XtSNE.png",width=800, height=800)
tsne.plot(nbt,pt.size = 1.5)
dev.off()
nbt=DBclust_dimension(nbt,1,2,reduction.use = "tsne",G.use = 8,set.ident = TRUE)
png("Seurat/XtSNE.cluster.png",width=800, height=800)
tsne.plot(nbt,pt.size = 1)
dev.off()
png("Seurat/XtSNE.treecluster.png",width=800, height=800)
nbt=buildClusterTree(nbt,do.reorder = TRUE,reorder.numeric = TRUE,pcs.use = 1:11)
dev.off()

png("Seurat/XtSNE.labelcluster.png",width=800, height=800)
tsne.plot(nbt,do.label=TRUE, label.pt.size=1, label.cex.text = 2)
dev.off()


my.data=fetch.data(nbt,c("ident","PC1","nGene"))#,"orig.ident"))
#nbt = set.all.ident(nbt,"orig.ident")
#markers.all=find_all_markers(nbt,thresh.test = 3,test.use = "roc", do.print = TRUE)
#markers.use=markers.all$gene
#png("Seurat/XHeatmap.png",width=800, height=2000)
#doHeatMap(nbt,genes.use = markers.use,slim.col.label = TRUE,remove.key = TRUE,cexRow=1)
#dev.off()



nbt = set.all.ident(nbt,"ident")
markers.all=find_all_markers(nbt,thresh.test = 3,test.use = "roc", do.print = TRUE)
markers.use=subset(markers.all,avg_diff>0&power>0.8)
write.table(markers.use,"Seurat/XCluster_Gene.xls",row.names=T, col.names=T, sep="\t",quote=FALSE)
#markers.use=markers.all$gene
png("Seurat/XHeatmapCluster.png",width=800, height=4000)
doHeatMap(nbt,genes.use = markers.use$gene,slim.col.label = TRUE,remove.key = TRUE,cexRow=1)
dev.off()

write.table(markers.all[markers.use$gene,],"Seurat/XCluster_Gene.xls",row.names=T, col.names=T, sep="\t",quote=FALSE)

