
mydata=read.csv("matrix.csv", header=TRUE)
library(CLAG)
M=mydata
RES=CLAG.clust(M, delta=0.05, threshold=0, analysisType=1, normalization="affine-global",rowId=row.names(M))
PCA <- prcomp(M)
clusterColors <- c("black", rainbow(RES$ncluster))
plot(PCA$x[,1], PCA$x[,2], col=clusterColors[RES$cluster+1], main=paste(RES$nclusters, "clusters"))
sink("clag-RES.txt")
lapply(RES,print)
sink()
sink("clag-clusters.txt")
lapply(RES$cluster,print)
sink()