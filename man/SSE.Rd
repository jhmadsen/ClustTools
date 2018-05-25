\name{SSE}
\alias{SSE}
\title{Sum of squared error (SSE) for cluster evaluation}
\description{Function to evaluate clustering results with sum of squared error (SSE) by calculating the distance from cluster members to cluster centroids}
\usage{
SSE(dataset, clusterVector)
}
\arguments{
  \item{dataset}{The dataset for which a sum of squared error and cluster centroids are returned}
  \item{clusterVector}{A vector of with integers indicating which cluster observations belong to}
}
\details{SSE computes the sum of squared error for clustering results, given a cluster vector. The smaller the squared error, the greater clustering results are achieved with respect to intra-cluster distance. SSE also return a matrix with cluster centroids.}
\value{
\item{centroidMatrix}{A matrix of n-clusters x n-dimensions with cluster centroids}
\item{clusterWithin}{A vector of n-clusters length with sum of squared error}
\item{sumWithin}{Sum of squared error for all clusters, i.e. sum(clusterWithin)}
}
\references{Tan, P.-N., Steinbach, M., Karpatne, A., & Kumar, V. (2005). Introduction to Data Mining (Second edition). ISBN: 978-03-213-2136-7}
\author{Jacob H. Madsen}
\examples{
## Select a dataset to standardize and cluster
X <- scale(iris[,1:4])

## Cluster the dataset with a given number of clusters
cluster.obj <- hclust(dist(X), method='complete')

## Cut the hierarchical clustering tree
chosen.clusters <- cutree(cluster.obj, 3)

## Evaluate the clustering results with 'SSE' and 'SST'
clusters.SSE <- SSE(X, chosen.clusters)$sumWithin
clusters.SST <- SST(X)

## Calculate the r-squared for your cluster solution
rsq <- 1-(clusters.SSE/clusters.SST)

print(rsq)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
