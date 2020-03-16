
add.one <- function (M) {
  for (i in 1:dim(M)[1]) {
    if (M[i,i] == 0) {
      M[i,i] <- M[i,i] + 1;
    }
  }
  return (M);
}
# Inflation step of MCL
inflate <- function (M,
                     inf) {
  M <- M^(inf);
  return (M);
}
# Normalize the matrix by column
norm <- function (M) {
  colum.sum <- apply(M,2,sum)
  M <- t(M) / colum.sum
  return (t(M))
}
# MCL procedure
mcl <- function (M, 	# Matrix
                 inf, 	# Inflation value
                 iter, 	# Number of iterations
                 verbose = F
) { 
  for (i in 1:iter) {
    old.M <- M;
    M.norm <- norm(M);
    M <- M.norm%*%M.norm;
    M <- inflate(M, inf);
    M <- norm(M);
    if (sum(old.M == M) == dim(M)[1]*dim(M)[2]) {
      break;
    }
    if (verbose) {
      print (paste ("iteration", i));
    } 
  }
  return (M);
}
collect.mcl.clusters <- function (M 	# Matrix (mcl result)
) {
  M.names <- row.names(M);
  clustered.nodes <- vector(mode = "logical", length = dim(M)[1])
  for (i in 1:dim(M)[1]) {
    nodes <- M.names[which(M[i,] != 0)];
    if (length(nodes) > 0 && !clustered.nodes[which(M[i,] != 0)]) {
      print (nodes);
      clustered.nodes[which(M[i,] != 0)] = T;
    }
  }
  return (clustered.nodes);
}
mcl.data=read.csv("matrix.csv", header=TRUE)
mcl.data <- as.matrix(mcl.data)
mcl.data=mcl.data*mcl.data
inf <- 4.0
mcl.clusters <- mcl(mcl.data,inf,200, verbose = T);
x=collect.mcl.clusters(mcl.clusters)
sink("MCL.txt")
lapply(mcl.clusters,print)
sink()
write.matrix(mcl.clusters, file = "MCL-clusters.txt", sep = " ")