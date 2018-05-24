\name{ClusterCor}
\alias{ClusterCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cluster correlation for cluster evaluation}
\description{Function to evaluate clustering results by calculating the correlation between an incidence matrix and distance matrix. Suggested by Tan, P.-N., Steinbach, M., Karpatne, A., & Kumar, V. (2005)}
\usage{
ClusterCor(dist.obj, clusterVector, return_matrices = FALSE)
}
\arguments{
  \item{dist.obj}{An object of class 'dist' for dataset}
  \item{clusterVector}{A vector with integers indicating which cluster observations belong to}
  \item{return_matrices}{Argument to return incidence- and distance matrix for observations}
}
\details{ClusterCor computes an incidence matrix for observations, given a cluster vector, by creating a n x n matrix where 1 is returned for observations in same cluster and 0 is returned for observations in different clusters. With a distance matrix as input the two matrices are vectorized and correlation is computed. A highly negative correlation indicates that observations in same cluster have small distance to each other, meaning good results with respect to minimizing intra-distance and maximizing inter-distance}
\value{
\item{incidenceMatrix}{Matrix with 0's and 1's indicating if observations belong to same cluster or not}
\item{distMatrix}{Matrix with distances between observations}
\item{correlation}{The correlation coefficient. The more negative correlation, the better results are achieved with respect to clustering objectives}
}
\references{Tan, P.-N., Steinbach, M., Karpatne, A., & Kumar, V. (2005). Introduction to Data Mining (Second edition). ISBN: 978-03-213-2136-7}
\author{Jacob H. Madsen}
\examples{
## Select a dataset to standardize and cluster
X <- scale(iris[,1:4])

## Cluster the dataset with a given number of clusters
cluster.obj <- kmeans(X, 3)

## Evaluate the clustering results with 'ClusterCor'
ClusterCor(dist(X), cluster.obj$cluster)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
