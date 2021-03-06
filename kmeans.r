
mydata=read.csv("matrix.csv", row.names=1)
mydata = na.omit(mydata)
mydata = scale(mydata) 
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:300) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
plot(1:300, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
fit <- kmeans(mydata, 30) 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
mydata2 <- data.frame(mydata, fit$cluster)
sink("C:\\ws\\kmeans-clusters.txt")
lapply(fit$cluster,print)
sink()