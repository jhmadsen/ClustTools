\name{SST}
\alias{SST}
\title{Sum of total squares for data}
\description{Function to calculate the total sum of squares (SST) for a given dataset}
\usage{
SST(dataset)
}
\arguments{
  \item{dataset}{The dataset for which total sum of squares is returned}
}
\details{SST computes the total sum of squares, being the distance for observations to the center of gravity, using euclidean distance}
\value{The total sum of squares (SST)}
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
