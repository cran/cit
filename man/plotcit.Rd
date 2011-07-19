\name{plotcit}
\alias{plotcit}
\title{
  Causal Inference Test Component Condition Plots
}
\description{
   This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus and a quantitative trait. The test is applicable to data that includes genotype (discrete), possible causal mediator such as gene expression (continuous) and an outcome of interest (continuous). Boxplots and a scatterplot to visualize the component tests are generated with p-values displayed.
}
\usage{
   plotcit(L, G, T, maxit=50000)
}

\arguments{
  \item{L}{
     Vector or nxp matrix of genotypes, coded {0,1,2}.
}
  \item{G}{
     Vector or nxp matrix of candidate causal mediators (continuous variable, such as gene expression).
}
  \item{T}{
     Vector or nxp matrix of continuous traits of interest.
}
  \item{maxit}{
      Maximum number of iterations to be conducted for the conditional independence test, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit.
}
}
\details{
  Increasing maxit will increase the precision of the cit p-value, which may be useful if a very small p-value is observed and precision is desired. However, increasing maxit increases the number of permutations conducted and therefore increases run time. For each test, component p-values are evaluated after 1000 permutations have been conducted for the conditional independence test in order to increase computational efficiency. At that point, if the maximum p-value of the 4 component tests is less than .02 then more permutations are conducted. There is a reevaluation after each permutation until at least 20 permutations result in F-statistics lower than that observed or until maxit is reached, whichever comes first. Although the L, G, and T, matrices must have the same number of rows, corresponding to the sample size, they may differ in the number of columns. Residuals from linear regression were used to plot conditioned features.
}
\value{
  Four plots, corresponding to the four component CIT tests are generated, labeled by p-values.
}
\references{
 Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
}
\author{
  Joshua Millstein
}

\examples{
# Sample Size
ss = 100

# Number of variables of each type
cols = 20

# Errors
e1 = matrix(rnorm(cols * ss),ncol=cols)
e2 = matrix(rnorm(cols * ss),ncol=cols)

# Simulate genotypes, gene expression, and clinical trait matrices
L = matrix(rbinom(cols * ss,2,.5),ncol=cols)
G =  matrix(.5*L + e1,ncol=cols)
T =  matrix(.2*G + e2,ncol=cols)

trios = cbind(1:cols,1:cols,1:cols)

plotcit(L[,1], G[,1], T[,1])

}

\keyword{ htest }
\keyword{ nonparametric }
