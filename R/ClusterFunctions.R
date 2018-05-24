## HOPKINS INDEX
Hopkins <- function(dataset, sample_size=0.1){

  dataset <- as.matrix(dataset)

  if(!is.numeric(dataset)){
    stop('dataset only takes numeric values')
  }
  if(sample_size>1|sample_size<=0){
    stop('sample size has to be greater than 0 and equal to or less than 1')
  }

  n <- (sample_size*nrow(dataset))
  dim <- ncol(dataset)

  minD <- apply(dataset, 2, min)
  maxD <- apply(dataset, 2, max)
  randomPoints <- matrix(0, ncol = dim, nrow = n)

  for(i in 1:dim) {
    randomPoints[, i] <- runif(n, min = minD[i], max = maxD[i])
  }

  randomSample <- sample(1:nrow(dataset), n, replace = FALSE)
  randomSample <- as.matrix(dataset[randomSample, ])
  colnames(randomSample) <- NULL
  rownames(randomSample) <- NULL

  oneMatrix <- rbind(randomPoints, randomSample)
  distMatrix <- as.matrix(dist(oneMatrix))

  nnRandomPoints <- NULL
  for(i in 1:n){
    sVector <- sort(distMatrix[i,(n+1):(n+n)])
    nnRandomPoints[i] <- sVector[1]
  }

  sumY <- sum(nnRandomPoints, na.rm = TRUE)

  nnRandomSample <- NULL
  for(i in (n+1):nrow(distMatrix)){
    sVector2 <- sort(distMatrix[i,(n+1):(n+n)])
    nnRandomSample[i] <- sVector2[2]
  }

  sumX <- sum(nnRandomSample, na.rm = TRUE)

  print(sumY/(sumX+sumY))

}

## EDGE CORRELATION
ClusterCor <- function(dist.obj, clusterVector, return_matrices=FALSE){

  if(class(dist.obj)!='dist'){
    stop('dist.obj needs to be of class "dist" - please use stats::dist() on your dataset')
  }

  distMatrix <- as.matrix(dist.obj)
  n <- nrow(distMatrix)

  if(length(clusterVector)!=n){
    stop('clusterVector has to be same length as observations number of observations in dataset')
  }

  incidenceMatrix <- matrix(nrow=n, ncol=n)

  for(i in 1:n){

    incidenceMatrix[i,] <- (clusterVector[i]==clusterVector)*1

  }

  correlation <- cor(c(incidenceMatrix), c(distMatrix))

  if(return_matrices==TRUE){

    return_list <- list(incidenceMatrix=incidenceMatrix, distMatrix=distMatrix, correlation=correlation)
    return(return_list)

  } else {

    return(correlation)

  }

}

## CLUSTER INDEX
ClusterIndex <- function(dist.obj, clusterVector){

  if(class(dist.obj)!='dist'){
    stop('dist.obj needs to be of class "dist" - please use stats::dist() on your dataset')
  }

  distMatrix <- as.matrix(dist.obj)
  n <- nrow(distMatrix)

  if(length(clusterVector)!=n){
    stop('clusterVector has to be same length as observations number of observations in dataset')
  }

  names(clusterVector) <- 1:n
  uniqueClust <- unique(clusterVector)

  distForeign <- NULL
  distCurrent <- NULL
  sizeForeign <- NULL
  sizeCurrent <- NULL

  for(i in 1:length(uniqueClust)){

    foreignClustObs <- as.numeric(names(clusterVector[clusterVector!=uniqueClust[i]]))
    nforeignClustObs <- length(foreignClustObs)

    currentClustObs <- as.numeric(names(clusterVector[clusterVector==uniqueClust[i]]))
    ncurrentClustObs <- length(currentClustObs)

    distForeign[i] <- sum(distMatrix[currentClustObs, foreignClustObs])
    distCurrent[i] <- sum(distMatrix[currentClustObs, currentClustObs])

    sizeForeign[i] <- nforeignClustObs*ncurrentClustObs
    sizeCurrent[i] <- (ncurrentClustObs-1)*ncurrentClustObs

  }

  avgDistForeign <- sum(distForeign)/(sum(sizeForeign))
  avgDistCurrent <- sum(distCurrent)/(sum(sizeCurrent))

  returnList <- list(avgDistCurrent=avgDistCurrent, avgDistForeign=avgDistForeign, clusterIndex=avgDistCurrent/avgDistForeign)
  return(returnList)

}

## SSE FOR HIERARCHICAL CLUSTERING
SSE <- function(dataset, clusterVector){

  dataset <- as.matrix(dataset)
  n <- nrow(dataset)
  dim <- ncol(dataset)

  if(!is.numeric(dataset)){
    stop('dataset only takes numeric values')
  }
  if(length(clusterVector)!=n){
    stop('clusterVector has to be the same length as number of observations in dataset')
  }

  sq.error <- function(x1){sum((x1 - centerGravity)^2)}

  names(clusterVector) <- 1:n
  uniqueCluster <- unique(clusterVector)
  clusterWithin <- NULL
  centroidMatrix <- matrix(NA, nrow=length(uniqueCluster), ncol = dim)

  for(i in 1:length(uniqueCluster)){

    clusterObs <- as.numeric(names(clusterVector[clusterVector==uniqueCluster[i]]))

    clusterSet <- dataset[clusterObs,]

    if(is.null(nrow(clusterSet))==TRUE){

      centerGravity <- clusterSet
      centroidMatrix[i,] <- clusterSet
      clusterWithin[i] <- 0

    } else {

      centerGravity <- colMeans(clusterSet)
      centroidMatrix[i,] <- centerGravity
      clusterWithin[i] <- sum(apply(clusterSet, 1, sq.error))

    }

  }
  return_list <- list(centroidMatrix=centroidMatrix, clusterWithin=clusterWithin, sumWithin=sum(clusterWithin))
  return(return_list)
}

## SST FOR HIERARCHICAL CLUSTERING
SST <- function(dataset){

  dataset <- as.matrix(dataset)

  if(!is.numeric(dataset)){
    stop('dataset only takes numeric values')
  }

  sq.error <- function(x1){sum((x1 - centerGravity)^2)}
  centerGravity <- colMeans(dataset)

  SST <- sum(apply(dataset, 1, sq.error))
  return(SST)

}
