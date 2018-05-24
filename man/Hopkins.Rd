\name{Hopkins}
\alias{Hopkins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hopkins statisic for cluster tendency}
\description{Function to compute the Hopkins statistic for datasets given a certain sample size. Indicates the cluster tendency in data}
\usage{
Hopkins(dataset, sample_size=0.1)
}
\arguments{
  \item{dataset}{The dataset for which a Hopkins statistic is returned}
  \item{sample_size}{The sample size as a proportion of the total number of observations in data. The greater the sample size, the more accurate Hopkins statistic is produced. Increased sample size has exponential increased complexity}
}
\details{The Hopkins statistic is useful as a test for cluster tendency in data. By creating a uniform distribution in data space, the distance to nearest original data point is calculated. The sum of distance to original data points is compared to sum of distance between original data points. The function returns an index between 0 and 1, where 1 characterize data partitioned in clusters, 0.5 characterize random uniformly distributed data and 0 characterize random data}
\value{The Hopkins statistic}
\references{Hopkins, B. (1954). A New Method for determining the Type of Distribution of Plant Individuals. Annals of Botany. Vol. 18, pp. 213–227}
\author{Jacob H. Madsen}
\examples{
## Select some arbitrary dataset
X <- iris[,1:4]

## Run Hopkins statistic
Hopkins(X, 0.4)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
