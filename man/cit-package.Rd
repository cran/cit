\name{cit-package}
\alias{cit-package}
\docType{package}
\title{
Causal Inference Test
}
\description{
This package implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus and a quantitative trait. The test is applicable to data that includes genotype (discrete), possible causal mediator such as gene expression (continuous) and an outcome of interest (continuous).
}
\details{
\tabular{ll}{
Package: \tab cit\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2013-12-19\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
}
This package implements a novel statistical framework in which existing notions of causal mediation are formalized into a hypothesis test, thus providing a standard quantitative measure of uncertainty in the form of a p-value. We treat the causal inference as a 'chain' of mathematical conditions that must be satisfied to conclude that the potential mediator is causal for the trait, where the inference is only as good as the weakest link in the chain. P-values are computed for the component conditions. The Intersection-Union Test, in which a series of statistical tests are combined to form an omnibus test, is then employed to generate the overall test result as the maximum of the component p-values. If we let L denote a locus of interest, T denote an outcome trait, and G denote a potential mediator, then the four component conditions are, 1) L and G are associated, 2) L and T are associated, 3) L is associated with G|T, and 4) L is independent of T|G. Test 4 requires an equivalence test and is implemented here using a permutation based approach (as described in Millstein et al. 2009).
}
\author{
Joshua Millstein

Maintainer: Joshua Millstein <millstei@usc.edu>
Joshua Millstein
}
\references{
Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
}

\keyword{ htest }
\keyword{ nonparametric }
