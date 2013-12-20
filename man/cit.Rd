\name{cit}
\alias{cit}
\title{
  Causal Inference Test
}
\description{
   This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus and a quantitative trait. The test is applicable to data that includes genotype (discrete), possible causal mediator such as gene expression (continuous) and an outcome of interest (continuous).
}
\usage{
   cit(L, G, T, trios = c(1,1,1), maxit=50000)
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
  \item{trios}{
      A matrix or dataframe of three columns. Each row represents a planned test to be conducted and the number of rows is equal to the total number of tests to be conducted. The first column is an indicator for the column in L, the second is an indicator for the column in G, and the third is an indicator for the column in T. Thus, the trios matrix defines the number of tests and the variables selected for each test.
}
  \item{maxit}{
      Maximum number of iterations to be conducted for the conditional independence test, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit.
}
}
\details{
  Increasing maxit will increase the precision of the cit p-value, which may be useful if a very small p-value is observed and precision is desired. However, increasing maxit increases the number of permutations conducted and therefore increases run time. For each test, component p-values are evaluated after 1000 permutations have been conducted for the conditional independence test in order to increase computational efficiency. At that point, if the maximum p-value of the 4 component tests is less than .02 then more permutations are conducted. There is a reevaluation after each permutation until at least 20 permutations result in F-statistics lower than that observed or until maxit is reached, whichever comes first. Although the L, G, and T, matrices must have the same number of rows, corresponding to the sample size, they may differ in the number of columns.
}
\value{
  A dataframe which includes the following columns:
  \item{L_index }{column of L used in the test}
  \item{G_index }{column of G used in the test}
  \item{T_index }{column of T used in the test}
  \item{p_cit }{CIT (omnibus) p-value}
  \item{p_TassocL }{component p-value for the test association of T and L.}
  \item{p_TassocGgvnL }{component p-value for the test association of T and G|L.}
  \item{p_GassocLgvnT }{component p-value for the test association of G and L|T.}
  \item{p_LindTgvnG }{component p-value for the equivalence test of L ind T|G}
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

results = cit(L, G, T, trios)

}

\keyword{ htest }
\keyword{ nonparametric }
