#include <R.h>
#include <Rmath.h>
#include <vector>
#include <algorithm>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>

using namespace std;

/*
L: matrix of genotypes
G: matrix of candidate causal mediators
T: matrix of traits
trios: matrix with 3 columns where each row contains 3 indicies pointing to the columns of L,G,T to test.
A single test (CIT) will be conducted for each row in 'trios'.
Programmer: Joshua Millstein
*/

extern "C" {

void citfun( int *, double *, double *, int *, 
	int *, int *, double *, double *, double *, double *, double *, int *, int * );
	
double gsl_stats_tss (const double data[], size_t stride, size_t n);

int randwrapper( int n );

int randwrapper( int n ) {

	int x ;
	x = (int)(n * unif_rand() );
	return x;

}

void citfun( int *L, double *G, double *T, int *nrow, 
	int *ncol, int *triosv, double *pval, double *pval1, double *pval2, double *pval3, double *pval4, int *ntest, int *maxit )
{
	int tst, lind, gind, tind, rw, cl, alleles, i, 
		rind, cind,  df1, df2, nobs, nvars, ip, npos,nperm;
	double rss, rss1, rss2, rss3, rss4, rss5, tss, F, Fp, a1, a2,
		pv, tmp, rhs, maxp, testval;
	bool aa, bb, cc, dd;
	const int firstloop = 1000;
	const int posno = 20;
	const double alpha = .1;
	vector<vector<int> > trios;
	vector<vector<int> > LL;
	vector<vector<double> > GG;
	vector<vector<double> > TT;
	vector<double> pvec;
	vector<double> gpred;
	vector<double> gresid;
	gsl_matrix *Lm, *cov, *X;
	gsl_vector *Gm, *Tm, *Gp, *c;

  /*
	Rprintf("ntest %d\n", *ntest);
	Rprintf("nrow %d\n", *nrow);
	Rprintf("ncol %d\n", *ncol);
	Rprintf("triosv %d\n", triosv[0]);
	Rprintf("triosv %d\n", triosv[1]);
	Rprintf("triosv %d\n", triosv[2]);
	Rprintf("pval %4.2f\n", *pval);
	*/

	LL.resize( *nrow );
	GG.resize( *nrow );
	TT.resize( *nrow );
	
	GetRNGstate();
	
	for(rw = 0; rw < *nrow; rw++) {
		LL[rw].resize( *ncol );
		GG[rw].resize( *ncol );
		TT[rw].resize( *ncol );	
	}

	for(cl = 0; cl < *ncol; cl++) {
		for(rw = 0; rw < *nrow; rw++) {
			LL[rw][cl] = L[rw + *nrow * cl];
			GG[rw][cl] = G[rw + *nrow * cl];
			TT[rw][cl] = T[rw + *nrow * cl];
		}
	}

	trios.resize( *ntest );
	for(rw = 0; rw < *ntest; rw++) {
		trios[rw].resize( 3 );
	}

	for(rw = 0; rw < *ntest; rw++) {
		for(cl = 0; cl < 3; cl++) {
			trios[rw][cl] = triosv[rw + *ntest * cl];
		}
	}

// subtract 1 to convert from R indexing to C++ indexing
	for(tst = 0; tst < *ntest; tst++) {
		lind = trios[tst][0] - 1;
		gind = trios[tst][1] - 1;
		tind = trios[tst][2] - 1;
	
// create analysis vectors w/no missing data
		nobs = 0;
		alleles = 0;
		for(rw = 0; rw < *nrow; rw++) {
			aa = LL[rw][lind] != -9999;
			bb = GG[rw][gind] != -9999;
			cc = TT[rw][tind] != -9999;
			if(aa && bb && cc) {
				if(LL[rw][lind] > alleles) alleles = LL[rw][lind]; // max no. variant alleles;
				nobs++;
			}
		}
		
		Lm = gsl_matrix_alloc (nobs,alleles + 1); // + 1 for intercept
		Gm = gsl_vector_alloc (nobs);
		Tm = gsl_vector_alloc (nobs);
		rind = 0;
		for(rw = 0; rw < *nrow; rw++) {
			aa = LL[rw][lind] != -9999;
			bb = GG[rw][gind] != -9999;
			cc = TT[rw][tind] != -9999;			
			
			if(aa && bb && cc) {
				// Code SNP {0,1,2} input variable L as 2 codominant indicator variables.
				if(alleles == 2) {
					a1 = 0.;
					a2 = 0.;
					if(LL[rw][lind] == 1) a1 = 1.;
					if(LL[rw][lind] > 0) a2 = 1.;
					gsl_matrix_set(Lm, rind, 0, 1.0);  // intercept
					gsl_matrix_set(Lm, rind, 1, a1);
					gsl_matrix_set(Lm, rind, 2, a2);
				}
				if(alleles == 1) {
					gsl_matrix_set(Lm, rind, 0, 1.0);  // intercept
					gsl_matrix_set(Lm, rind, 1, LL[rw][lind]);
				}
				gsl_vector_set(Gm, rind, GG[rw][gind]);
				gsl_vector_set(Tm, rind, TT[rw][tind]);
				rind++;
			}
		}
		
		nvars = alleles;
		ip = nvars + 1;
		
//		Rprintf("nvars %d\n", alleles);
//		Rprintf("ip %d\n", ip);
		
		// fit model T ~ L1 + L2
		c = gsl_vector_alloc (ip);
		cov = gsl_matrix_alloc (ip, ip);
		gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (nobs, ip);
		gsl_multifit_linear (Lm, Tm, c, cov, &rss, work);
		gsl_multifit_linear_free (work);
		gsl_matrix_free (cov);
		gsl_vector_free (c);
	//		tss = gsl_stats_tss(Tnm, 1, nobs);
		tss = gsl_stats_tss(Tm->data, Tm->stride, Tm->size);
		df1 = nvars;
		df2 = nobs - ip;
		F = df2*(tss-rss)/(rss*df1);
		pv = gsl_cdf_fdist_Q(F, df1, df2);
		/*
		Rprintf("tst %d\n", tst);
		Rprintf("lind %d\n", lind);
		Rprintf("gind %d\n", gind);
		Rprintf("tind %d\n", tind);
		Rprintf("df1 %d\n", df1);
		Rprintf("df2 %d\n", df2);
		Rprintf("F %4.6f\n", F);
		*/
		pvec.push_back( pv ); // pval for T ~ L
//		Rprintf("pval1: %4.6f\n", log10 (pv));
		
		// fit model T ~ G + L1 + L2
		X = gsl_matrix_alloc (nobs, ip + 1);
		for(rw = 0; rw < nobs; rw++) {
			for(cind = 0; cind < ip; cind++) {
				gsl_matrix_set(X, rw, cind, gsl_matrix_get (Lm, rw, cind));
			}
			gsl_matrix_set(X, rw, ip, gsl_vector_get (Gm, rw));
		}
		c = gsl_vector_alloc (ip + 1);
		cov = gsl_matrix_alloc (ip + 1, ip + 1);
		work = gsl_multifit_linear_alloc (nobs, ip + 1);
		gsl_multifit_linear (X, Tm, c, cov, &rss1, work);
		gsl_multifit_linear_free (work);
		gsl_matrix_free (X);
		gsl_matrix_free (cov);
		gsl_vector_free (c);
		df1 = 1;
		df2 = nobs - ip -1;
		F = df2*(rss-rss1)/(rss1*df1);
		pv = gsl_cdf_fdist_Q(F, df1, df2);
		pvec.push_back( pv ); // pval for T|G ~ L

//		Rprintf("pval2: %4.6f\n", log10 (pv));
		
		// fit model G ~ T
		X = gsl_matrix_alloc (nobs,2);
		for(rw = 0; rw < nobs; rw++) {
			gsl_matrix_set(X, rw, 0, 1.0);  // intercept
			gsl_matrix_set(X, rw, 1, gsl_vector_get (Tm, rw));
		}
		c = gsl_vector_alloc (2);
		cov = gsl_matrix_alloc (2, 2);
		work = gsl_multifit_linear_alloc (nobs, 2);
		gsl_multifit_linear (X, Gm, c, cov, &rss2, work);
		gsl_multifit_linear_free (work);
		gsl_matrix_free (X);
		gsl_matrix_free (cov);
		gsl_vector_free (c);

		// fit model G ~ T + L1 + L2
		X = gsl_matrix_alloc (nobs, ip + 1);
		for(rw = 0; rw < nobs; rw++) {
			for(cind = 0; cind < ip; cind++) {
				gsl_matrix_set(X, rw, cind, gsl_matrix_get (Lm, rw, cind));
			}
			gsl_matrix_set(X, rw, ip, gsl_vector_get (Tm, rw));
		}
		c = gsl_vector_alloc (ip + 1);
		cov = gsl_matrix_alloc (ip + 1, ip + 1);
		work = gsl_multifit_linear_alloc (nobs, ip + 1);
		gsl_multifit_linear (X, Gm, c, cov, &rss3, work);
		gsl_multifit_linear_free (work);
		gsl_matrix_free (X);
		gsl_matrix_free (cov);
		gsl_vector_free (c);
		df1 = alleles;
		df2 = nobs - ip -1;
		F = df2*(rss2-rss3)/(rss3*df1);
		pv = gsl_cdf_fdist_Q(F, df1, df2);
		pvec.push_back( pv ); // pval for G|L ~ T
		
//		Rprintf("pval3: %4.6f\n", log10 (pv));
		
		// fit model T ~ G
		X = gsl_matrix_alloc (nobs,2);
		for(rw = 0; rw < nobs; rw++) {
			gsl_matrix_set(X, rw, 0, 1.0);  // intercept
			gsl_matrix_set(X, rw, 1, gsl_vector_get (Gm, rw));
		}
		c = gsl_vector_alloc (2);
		cov = gsl_matrix_alloc (2, 2);
		work = gsl_multifit_linear_alloc (nobs, 2);
		gsl_multifit_linear (X, Tm, c, cov, &rss4, work);
		gsl_multifit_linear_free (work);
		gsl_matrix_free (X);
		gsl_matrix_free (cov);
		gsl_vector_free (c);
		df1 = alleles;
		df2 = nobs - ip -1;
		F = df2*(rss4-rss1)/(rss1*df1);
		
		// fit model G ~ L1 + L2
		c = gsl_vector_alloc (ip);
		cov = gsl_matrix_alloc (ip, ip);
		work = gsl_multifit_linear_alloc (nobs, ip);
		gsl_multifit_linear (Lm, Gm, c, cov, &rss5, work);
		gsl_multifit_linear_free (work);
		gsl_matrix_free (cov);
		
		// residuals for G ~ L1 + L2
		for(rw = 0; rw < nobs; rw++) {
			rhs = 0;
			for(cl = 0; cl < ip; cl++) {
				rhs += gsl_vector_get (c, cl) * gsl_matrix_get (Lm, rw, cl);
			}
			gpred.push_back(rhs);
			tmp = gsl_vector_get (Gm, rw) - rhs;
			gresid.push_back(tmp);
		}
		gsl_vector_free (c);
		
		Gp = gsl_vector_alloc (nobs);
		npos = 0;
		for(i = 0; i < firstloop; i++){
			// randomly permute residuals
			random_shuffle( gresid.begin(), gresid.end(), randwrapper );
			
			// compute G* based on marginal L effects and permuted residuals
			for(rw = 0; rw < nobs; rw++) {
				gsl_vector_set(Gp, rw, gpred[rw] + gresid[rw] );
			}
			
			// Recompute F based on G*
			// fit model T ~ G*
			X = gsl_matrix_alloc (nobs,2);
			for(rw = 0; rw < nobs; rw++) {
				gsl_matrix_set(X, rw, 0, 1.0);  // intercept
				gsl_matrix_set(X, rw, 1, gsl_vector_get (Gp, rw));
			}
			c = gsl_vector_alloc (2);
			cov = gsl_matrix_alloc (2, 2);
			work = gsl_multifit_linear_alloc (nobs, 2);
			gsl_multifit_linear (X, Tm, c, cov, &rss4, work);
			gsl_multifit_linear_free (work);
			gsl_matrix_free (X);
			gsl_matrix_free (cov);
			gsl_vector_free (c);
			
			// fit model T ~ G* + L1 + L2
			X = gsl_matrix_alloc (nobs, ip + 1);
			for(rw = 0; rw < nobs; rw++) {
				for(cind = 0; cind < ip; cind++) {
					gsl_matrix_set(X, rw, cind, gsl_matrix_get (Lm, rw, cind));
				}
				gsl_matrix_set(X, rw, ip, gsl_vector_get (Gp, rw));
			}
			c = gsl_vector_alloc (ip + 1);
			cov = gsl_matrix_alloc (ip + 1, ip + 1);
			work = gsl_multifit_linear_alloc (nobs, ip + 1);
			gsl_multifit_linear (X, Tm, c, cov, &rss1, work);
			gsl_multifit_linear_free (work);
			gsl_matrix_free (X);
			gsl_matrix_free (cov);
			gsl_vector_free (c);
			
			df1 = alleles;
			df2 = nobs - ip - 1;
			Fp = df2*(rss4-rss1)/(rss1*df1);
			if(Fp < F) npos++;
		} // end permutation loop
		
		maxp = *max_element( pvec.begin(), pvec.end() );
		nperm = firstloop;
		aa = npos < posno;
		bb = maxp < alpha;
		cc = nperm < *maxit;
		maxp = *max_element( pvec.begin(), pvec.end());
		testval = (double) (npos + 1) / nperm ;
		dd = maxp < testval; // check that other component p-values are small
		
		if(aa && bb && cc && dd){
			while(aa && cc) {
				
				// randomly permute residuals
				random_shuffle( gresid.begin(), gresid.end(), randwrapper );
				
				// compute G* based on marginal L effects and permuted residuals
				for(rw = 0; rw < nobs; rw++) {
					gsl_vector_set(Gp, rw, gpred[rw] + gresid[rw] );
				}
				
				// Recompute F based on G*
				// fit model T ~ G*
				X = gsl_matrix_alloc (nobs,2);
				for(rw = 0; rw < nobs; rw++) {
					gsl_matrix_set(X, rw, 0, 1.0);  // intercept
					gsl_matrix_set(X, rw, 1, gsl_vector_get (Gp, rw));
				}
				c = gsl_vector_alloc (2);
				cov = gsl_matrix_alloc (2, 2);
				gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (nobs, 2);
				gsl_multifit_linear (X, Tm, c, cov, &rss4, work);
				gsl_multifit_linear_free (work);
				gsl_matrix_free (X);
				gsl_matrix_free (cov);
				gsl_vector_free (c);
				
				// fit model T ~ G* + L1 + L2
				X = gsl_matrix_alloc (nobs, ip + 1);
				for(rw = 0; rw < nobs; rw++) {
					for(cind = 0; cind < ip; cind++) {
						gsl_matrix_set(X, rw, cind, gsl_matrix_get (Lm, rw, cind));
					}
					gsl_matrix_set(X, rw, ip, gsl_vector_get (Gp, rw));
				}
				c = gsl_vector_alloc (ip + 1);
				cov = gsl_matrix_alloc (ip + 1, ip + 1);
				work = gsl_multifit_linear_alloc (nobs, ip + 1);
				gsl_multifit_linear (X, Tm, c, cov, &rss1, work);
				gsl_multifit_linear_free (work);
				gsl_matrix_free (X);
				gsl_matrix_free (cov);
				gsl_vector_free (c);
				
				df1 = alleles;
				df2 = nobs - ip - 1;
				Fp = df2*(rss4-rss1)/(rss1*df1);
				
				if(Fp < F) npos++;
				aa = npos < posno;
				cc = nperm < ( *maxit - 1 );
				nperm++;
			} // end 'while' permutation loop
		} // End if
		pv = (double) npos / nperm;
		pvec.push_back(pv); // pval for L ind T|G
		maxp = *max_element( pvec.begin(), pvec.end() );
		pval[tst] = maxp;

		pval1[tst] = pvec[0]; // pval for T ~ L
		pval2[tst] = pvec[1]; // pval for T|L ~ G
		pval3[tst] = pvec[2]; // pval for G|T ~ L
		pval4[tst] = pvec[3]; // pval for L ind T|G

		
		pvec.clear();
		gresid.clear();
		gpred.clear();
		gsl_matrix_free (Lm);
		gsl_vector_free (Gm);
		gsl_vector_free (Tm);
		gsl_vector_free (Gp);
		
	} // End tst loop
	PutRNGstate();
	trios.clear();
	LL.clear();
	GG.clear();
	TT.clear();
} // End citfun
} // End extern c wrapper
