\name{fdr.cit}
\alias{fdr.cit}
\title{
  Omnibus FDR Values for CIT
}
\description{
   Results including permutation results from cit.cp or cit.bp from multiple tests are used in a permutation-based approach to compute false discovery rates (q-values) for each test. q-values are computed for the omnibus test as well as all component tests. q-values are returned with confidence intervals to quantify uncertainty in the estimate. 
}
\usage{
   fdr.cit( cit.perm.list, cl=.95, c1=NA )
}

\arguments{
  \item{cit.perm.list}{
     List, where each item is a dataframe output from cit.cp or cit.bp, which must be run with n.perm > 0. 
}
  \item{cl}{
     Confidence level for the q-value confidence interval (default is .95).
}
  \item{c1}{
     Over-dispersion parameter. Setting the over-dispersion parameter to one corresponds to the assumption that all tests (referring to multiple omnibus tests conducted using different sets of variables) are independent. The default (c1=NA) is to use an empirically estimated over-dispersion parameter based on the permutation results.
}
}
\details{
  For the initial cit.cp and cit.bp test computations, we suggest setting n.perm to 100 or greater for moderate numbers of tests, say 10s or 100s. Larger numbers of permutations will generate more precise FDR estimates but require substantial computation time. If a large number of tests, say 1000s, have been conducted, then n.perm can be as small as 10 and still generate good FDR estimates and confidence intervals. FDR confidence intervals are especially useful to quantify uncertainty when FDR value > 0.05. The over-dispersion parameter should be set to one if it is known that all tests are independent, however, this is rarely the case in typical multiple testing settings. The fdr.cit function can also be used with simulated data in order to generate power calculations. In this case the over-dispersion parameter may be estimated from pilot data or other datasets and fixed to the prior estimate using the c1 argument. The omnibus CIT FDR value is generated as the maximum FDR across the four component tests, in an intesection-union type of approach analogous to the method used to generate the omnibus CIT p-value in cit.cp and cit.bp.
}
\value{
  A dataframe which includes the following columns:
  \item{ perm }{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
  \item{p_cit }{CIT (omnibus) p-value}
  \item{p_TassocL }{component p-value for the test of association between T and L.}
  \item{p_TassocGgvnL }{component p-value for the test of association between T and G|L.}
  \item{p_GassocLgvnT }{component p-value for the test of association between G and L|T.}
  \item{p_LindTgvnG }{component p-value for the equivalence test of L ind T|G}
}
\references{
 Millstein J, Volfson D. 2013. Computationally efficient permutation-based confidence interval estimation for tail-area FDR. Frontiers in Genetics | Statistical Genetics and Methodology 4(179):1-11.
}
\author{
  Joshua Millstein
}

\examples{
# Sample Size
ss = 100

n.tests = 20
myresults = vector('list', n.tests)

for( tst in 1:n.tests ){

	# Errors
	e1 = matrix(rnorm(ss),ncol=1)
	e2 = matrix(rnorm(ss),ncol=1)

	# Simulate genotypes, gene expression, and clinical traits
	L = matrix(rbinom(ss,2,.5),ncol=1)
	G =  matrix(.5*L + e1,ncol=1)
	T =  matrix(.3*G + e2,ncol=1)
	T = ifelse( T > median(T), 1, 0 )

	myresults[[ tst ]] = cit.bp(L, G, T, n.perm=20)
}
fdr.cit( myresults )

}

\keyword{ nonparametric }
