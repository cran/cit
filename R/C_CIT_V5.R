######################################################################
# Program Name: C_CIT_V5.R
# Purpose: R frontend of C++ CIT function
# Programmer: Joshua Millstein
# Date: 12/9/2011
#
# Input:
#   L: vector or nxp matrix of genotypes {0,1,2}
#   G: vector or nxp matrix of candidate causal mediators.
#   T: vector or nxp matrix of traits
#   trios: A matrix or dataframe of three columns. Each row represents a planned test to be conducted 
#          and the number of rows is equal to the total number of tests. The first column is an
#          indicator for the column in L, the second is an indicator for the column in G, and the third
#          is an indicator for the column in T.
#dyn.load(paste("cit", .Platform$dynlib.ext, sep = ""))

cit = function(L, G, T, trios = c(1,1,1), maxit=50000) {

	if(is.vector(L)) {
	   L = matrix(L,ncol=1)
	} else {
	   L = as.matrix(L)
	}
	if(is.vector(G)) {
	   G = matrix(G,ncol=1)
	} else {
	   G = as.matrix(G)
	}
	if(is.vector(T)) {
	   T = matrix(T,ncol=1)
	} else {
	   T = as.matrix(T)
	}
	 
	if(dim(L)[1] != dim(G)[1]) stop("L rows != G rows")
	if(dim(L)[1] != dim(T)[1]) stop("L rows != T rows")
	if(dim(G)[1] != dim(T)[1]) stop("G rows != T rows")

	# Recode NA's to -9999
	ms_f = function(mat) {
	   for(c_ in 1:dim(mat)[2]){
	      mat[is.na(mat[,c_]),c_] = -9999
	   }
	   return(mat);
	}
	L = ms_f(L)
	G = ms_f(G)
	T = ms_f(T)

	trios = as.matrix(trios)
	trios = na.exclude(trios)
	if(length(trios) == 3) trios = matrix(trios,nrow=1)

	pval = rep(1,dim(trios)[1])
	pval1 = rep(1,dim(trios)[1])  # output component p-values
	pval2 = rep(1,dim(trios)[1])  # output component p-values
	pval3 = rep(1,dim(trios)[1])  # output component p-values
	pval4 = rep(1,dim(trios)[1])  # output component p-values
	ntest = length(pval)

	nrow = dim(L)[1]
	ncol = dim(L)[2]

	tmp = .C("citfun", as.integer(L), as.double(G), as.double(T), as.integer(nrow), 
	   as.integer(ncol), as.integer(trios), as.double(pval), as.double(pval1), as.double(pval2), as.double(pval3), as.double(pval4),
	   as.integer(ntest), as.integer(maxit));

   rslts = as.data.frame(matrix(NA,nrow=ntest,ncol=8))
   names(rslts) = c("L_index", "G_index", "T_index", "p_cit", "p_TassocL", "p_TassocGgvnL", "p_GassocLgvnT", "p_LindTgvnG")
   rslts[,1:3] = trios
   for(i in 4:8) rslts[,i] = tmp[[i+3]]

   return(rslts)

} # End cit function

plotcit = function(L, G, T, maxit=50000) {

	L = matrix(L,ncol=1)
	G = matrix(G,ncol=1)
	T = matrix(T,ncol=1)

	if(dim(L)[1] != dim(G)[1]) stop("L elements != G elements")
	if(dim(L)[1] != dim(T)[1]) stop("L elements != T elements")
	if(dim(G)[1] != dim(T)[1]) stop("G elements != T elements")

	# Recode NA's to -9999
	ms_f = function(mat) {
	   for(c_ in 1:dim(mat)[2]){
	      mat[is.na(mat[,c_]),c_] = -9999
	   }
	   return(mat);
	}
	L = ms_f(L)
	G = ms_f(G)
	T = ms_f(T)

	trios = matrix(c(1,1,1),nrow=1)

	pval = 1
	pval1 = 1  # output component p-values
	pval2 = 1  # output component p-values
	pval3 = 1  # output component p-values
	pval4 = 1  # output component p-values
	ntest = 1

	nrow = dim(L)[1]
	ncol = dim(L)[2]

	tmp = .C("citfun", as.integer(L), as.double(G), as.double(T), as.integer(nrow), 
	   as.integer(ncol), as.integer(trios), as.double(pval), as.double(pval1), as.double(pval2), as.double(pval3), as.double(pval4),
	   as.integer(ntest), as.integer(maxit));

   pvals = rep(NA,5)
   pval_nms = c("p_CIT", "p_T~L", "p_T~G|L", "p_G~L|G", "p_LindT|G")
   for(i in 1:5) pvals[i] = tmp[[i+6]]

  fit = lm(T ~ as.factor(L),na.action="na.exclude")
  TcL = resid(fit)
  fit = lm(G ~ T,na.action="na.exclude")
  GcT = resid(fit)
  fit = lm(T ~ G,na.action="na.exclude")
  TcG = resid(fit)

  par(mfrow=c(2,2),oma=c(1,1,2,1))

  boxplot(T~L,xlab="genotype",ylab="T",main=paste(pval_nms[2],"=",round(pvals[2],3)),notch=FALSE)
  plot(TcL,G,xlab="T|L",ylab="G",main=paste(pval_nms[3],"=",round(pvals[3],3)))
  fit = lm(G ~ TcL,na.action="na.exclude")
  abline(fit$coef[1],fit$coef[2])
  boxplot(GcT~L,xlab="genotype",ylab="G|T",main=paste(pval_nms[3],"=",round(pvals[3],3)),notch=FALSE)
  boxplot(TcG~L,xlab="genotype",ylab="T|G",main=paste(pval_nms[5],"=",round(pvals[5],3)),notch=FALSE)
  mtext(paste(pval_nms[1],"=",round(pvals[1],3)),side=3,line=0,outer=TRUE)

} # End plotcit function

