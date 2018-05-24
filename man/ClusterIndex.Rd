\name{ClusterIndex}
\alias{ClusterIndex}
\title{Cluster index for cluster evaluation}
\description{Function to evaluate clustering results by calculating the ratio between intra-cluster and inter-cluster distance}
\usage{
ClusterIndex(dist.obj, clusterVector)
}
\arguments{
  \item{dist.obj}{An object of class 'dist' for dataset}
  \item{clusterVector}{A vector with integers indicating which cluster observations belong to}
}
\details{ClusterIndex computes the average distance for all observations, to observations in same cluster, and compare it with the average distance to observations in other clusters. The cluster index ranges from 0 to 1, where 0 indicates short distance to observations in same cluster, when compared to distance to observations on other clusters. 1 indicates equal distance to observations in same cluster as to observations in other clusters. The smaller the cluster index, the greater results are achieved with respect to intra- and inter-cluster distance. The theory is an extension of Tan, P.-N., Steinbach, M., Karpatne, A., & Kumar, V. (2005) idea of cohesion and separation}
\value{
\item{avgDistCurrent}{Average distance for all observations to other observations in same cluster}
\item{avgDistForeign}{Average distance for all observations to other observations in other clusters}
\item{clusterindex}{The ratio between avgDistCurrent and avgDistForeign, i.e. avgDistCurrent/avgDistForeign}
}
\author{Jacob H. Madsen}
\examples{
## Select a dataset to standardize and cluster
X <- scale(iris[,1:4])

## Cluster the dataset with a given number of clusters
cluster.obj <- kmeans(X, 3)

## Evaluate the clustering results with 'ClusterIndex'
ClusterIndex(dist(X), cluster.obj$cluster)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
