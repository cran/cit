#include <R.h>
#include <Rmath.h>
#include <vector>
#include <algorithm>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <iostream>
#include <random> 

#include "logisticfunc.h"

using namespace std;

/*
L: matrix of continuous instrumental variables
G: matrix of candidate causal mediators
T: matrix of 0/1 variables
Programmer: Joshua Millstein
*/

extern "C" {

void citconlog2( double *, double *, double *, int &,
	int &, double &, double &, double &, double &, double &, int &);

double gsl_stats_tss (const double data[], size_t stride, size_t n);

int randwrapper1( int n );

int randwrapper1( int n )
{
	int x ;
	x = (int)(n * unif_rand() );
	return x;
}


void citconlog2( double *L, double *G, double *T, int &nrow,
	int &ncol, double &pval, double &pval1, double &pval2, double &pval3, double &pval4, int &maxit )
{
	int rw, cl, i, rind, df, df1, df2, nobs, ip, npos, nperm, nmiss, stride;
	double rss2, rss3, rss5, F, pv, pvp, tmp, rhs, maxp, testval;
	bool aa, bb, cc, dd, converged;
	const int firstloop = 1000;
	const int posno = 20;
	const double alpha = .1;
	vector<vector<double> > LL;
	vector<double> pvec;
	vector<double> gpred;
	vector<double> gresid;

	gsl_matrix *Lm, *cov, *X;
	gsl_vector *Gm, *Tm, *Gp, *c;
	unsigned se = 10;

	double *designmat = new double[ nrow * (ncol + 2) ];
	double *phenovec = new double[ nrow ];

	LL.resize( nrow );
	GetRNGstate();

	for(rw = 0; rw < nrow; rw++) {
		LL[rw].resize( ncol );
	}

	for(cl = 0; cl < ncol; cl++) {
		for(rw = 0; rw < nrow; rw++) {
			LL[rw][cl] = L[rw + nrow * cl];
		}
	}


// create analysis vectors w/no missing data
		nobs = 0;
		for(rw = 0; rw < nrow; rw++) {
		     nmiss = 0;
		     for(cl = 0; cl < ncol; cl++) {
		        if( LL[rw][cl] == -9999 ) {
					nmiss++;
			    }
		     }
			aa = nmiss == 0;
			bb = G[rw] != -9999;
			cc = T[rw] != -9999;
			if(aa && bb && cc) {
				nobs++;
			}
		}

		Lm = gsl_matrix_alloc (nobs, ncol);
		Gm = gsl_vector_alloc (nobs);
		Tm = gsl_vector_alloc (nobs);
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			nmiss = 0;
		   for(cl = 0; cl < ncol; cl++) {
		        if( LL[rw][cl] == -9999 ) {
					nmiss++;
			    }
		   }
			aa = nmiss == 0;
			bb = G[rw] != -9999;
			cc = T[rw] != -9999;

			if(aa && bb && cc) {
				for(cl = 0; cl < ncol; cl++) {
                  	gsl_matrix_set(Lm, rind, cl, LL[rw][cl]);
		      	}
				gsl_vector_set(Gm, rind, G[rw]);
				gsl_vector_set(Tm, rind, T[rw]);
				rind++;
			}
		}
		// fit model T ~ L
		ip = 1 + ncol;                               // intercept + multiple L variable
		for(rw = 0; rw < nobs; rw++) {
		   phenovec[ rw ] = gsl_vector_get(Tm, rw );
		   designmat[ rw * ip  ] = 1;      // intercept
			for(cl = 0; cl < ncol; cl++) {
                  designmat[ rw * ip + 1 + cl  ]  = gsl_matrix_get (Lm, rw, cl);
		     }
		}
		df = ncol;
		converged = logisticReg( pv, phenovec, designmat, nobs, ip, df );
		pv = ( converged ) ? pv : 9;
		pvec.push_back( pv );  // pval for T ~ L, 9 if it did not converge, p1

		// fit model T ~ L + G
		stride = ip + 1;
		for(rw = 0; rw < nobs; rw++) {
		   designmat[ rw * stride ] = 1;      // intercept
			for(cl = 0; cl < ncol; cl++) {
          	designmat[ rw * stride + 1 + cl  ]  = gsl_matrix_get (Lm, rw, cl);
		   }
		   designmat[ rw * stride + 1 + ncol  ] = gsl_vector_get(Gm, rw );
		}

		df = 1;
		converged = logisticReg( pv, phenovec, designmat, nobs, stride, df );
		pv = ( converged ) ? pv : 9;
		pvec.push_back( pv );  // pval for T ~ G|L, 9 if it did not converge, p2

		// fit model G ~ T
		X = gsl_matrix_alloc (nobs,2);
		for(rw = 0; rw < nobs; rw++) {
			gsl_matrix_set(X, rw, 0, 1.0);  // intercept
			gsl_matrix_set(X, rw, 1, gsl_vector_get (Tm, rw));
		}
		c = gsl_vector_alloc (2);
		cov = gsl_matrix_alloc (2, 2);
		gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (nobs, 2);
		gsl_multifit_linear (X, Gm, c, cov, &rss2, work);
		gsl_multifit_linear_free (work);
		gsl_matrix_free (X);
		gsl_matrix_free (cov);
		gsl_vector_free (c);

		// fit model G ~ L + T
		X = gsl_matrix_alloc (nobs, ip + 1);
		for(rw = 0; rw < nobs; rw++) {
			gsl_matrix_set(X, rw, 0, 1.0);      // intercept
			for(cl = 0; cl < ncol; cl++) {
                  gsl_matrix_set(X, rw, cl + 1, gsl_matrix_get (Lm, rw, cl));
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
		df1 = ncol;
		df2 = nobs - ip -1;
		F = df2*(rss2-rss3)/(rss3*df1);
		pv = gsl_cdf_fdist_Q(F, df1, df2);
		pvec.push_back( pv ); // pval for G ~ L|T, p3

		// fit model T ~ G + L to test L
		stride = ip + 1;
		for(rw = 0; rw < nobs; rw++) {
		   designmat[ rw * stride  ] = 1;      // intercept
		   designmat[ rw * stride + 1  ] = gsl_vector_get(Gm, rw );
			for(cl = 0; cl < ncol; cl++) {
          	designmat[ rw * stride + 2 + cl  ]  = gsl_matrix_get (Lm, rw, cl);
		   }
		}
		df = ncol;
		converged = logisticReg( pv, phenovec, designmat, nobs, stride, df );
		pv = ( converged ) ? pv : 9;    // p-value for T ~ L|G

		// fit model G ~ L
		X = gsl_matrix_alloc (nobs, ip );
		for(rw = 0; rw < nobs; rw++) {
			gsl_matrix_set(X, rw, 0, 1.0);      // intercept
			for(cl = 0; cl < ncol; cl++) {
                  gsl_matrix_set(X, rw, cl + 1, gsl_matrix_get (Lm, rw, cl));
		     }
		}
		c = gsl_vector_alloc (ip);
		cov = gsl_matrix_alloc (ip, ip);
		work = gsl_multifit_linear_alloc (nobs, ip);
		gsl_multifit_linear (X, Gm, c, cov, &rss5, work);
		gsl_multifit_linear_free (work);
		gsl_matrix_free (cov);

		// residuals for G ~ L
		for(rw = 0; rw < nobs; rw++) {
			rhs = 0;
			for(cl = 0; cl < ip; cl++) {
                  rhs += gsl_vector_get (c, cl) * gsl_matrix_get (X, rw, cl);
		     }

			gpred.push_back(rhs);
			tmp = gsl_vector_get (Gm, rw) - rhs;
			gresid.push_back(tmp);
		}
		gsl_vector_free (c);

		// Conduct an initial set of permutations

		Gp = gsl_vector_alloc (nobs);
		npos = 0;
		for(i = 0; i < firstloop; i++){
			// randomly permute residuals
			shuffle( gresid.begin(), gresid.end(), std::default_random_engine(se) );	

			// compute G* based on marginal L effects and permuted residuals
			for(rw = 0; rw < nobs; rw++) {
				gsl_vector_set(Gp, rw, gpred[rw] + gresid[rw] );
			}

			// Recompute p-value for T ~ L|G based on G*
			// fit model T ~ G* + L to test L
			stride = ip + 1;
			for(rw = 0; rw < nobs; rw++) {
		   		designmat[ rw * stride  ] = 1;      // intercept
		   		designmat[ rw * stride + 1  ] = gsl_vector_get(Gp, rw );
				for(cl = 0; cl < ncol; cl++) {
          		designmat[ rw * stride + 2 + cl  ]  = gsl_matrix_get (Lm, rw, cl);
		   		}
			}

			df = ncol;
			converged = logisticReg( pvp, phenovec, designmat, nobs, stride, df );
			pvp = ( converged ) ? pvp : 9;    // p-value for T ~ L|G*
			if( pvp > pv ) npos++;

		} // end initial permutation loop

		// Conduct additional permutations if there is some indication of statistical significance
		maxp = *max_element( pvec.begin(), pvec.end() );
		nperm = firstloop;
		aa = npos < posno;
		bb = maxp < alpha;
		cc = nperm < maxit;
		maxp = *max_element( pvec.begin(), pvec.end());
		testval = (double) (npos + 1) / nperm ;
		dd = maxp < testval; // check that other component p-values are small

		if(aa && bb && cc && dd){
			while(aa && cc) {

				// randomly permute residuals
			    shuffle( gresid.begin(), gresid.end(), std::default_random_engine(se) );	

				// compute G* based on marginal L effects and permuted residuals
				for(rw = 0; rw < nobs; rw++) {
					gsl_vector_set(Gp, rw, gpred[rw] + gresid[rw] );
				}

				// Recompute p-value for T ~ L|G based on G*
				// fit model T ~ G* + L to test L
				stride = ip + 1;
				for(rw = 0; rw < nobs; rw++) {
		   			designmat[ rw * stride  ] = 1;      // intercept
		   			designmat[ rw * stride + 1  ] = gsl_vector_get(Gp, rw );
					for(cl = 0; cl < ncol; cl++) {
          			designmat[ rw * stride + 2 + cl  ]  = gsl_matrix_get (Lm, rw, cl);
		   			}
				}

				df = ncol;
				converged = logisticReg( pvp, phenovec, designmat, nobs, stride, df );
				pvp = ( converged ) ? pvp : 9;    // p-value for T ~ L|G*
				if( pvp > pv ) npos++;

				aa = npos < posno;
				cc = nperm < ( maxit - 1 );
				nperm++;
			} // end 'while' permutation loop
		} // End if
		pv = 1.0 * npos / nperm;
		pvec.push_back(pv); // pval for L ind T|G
		maxp = *max_element( pvec.begin(), pvec.end() );
		pval = maxp;

		pval1 = pvec[0]; // pval for T ~ L
		pval2 = pvec[1]; // pval for T ~ G|L
		pval3 = pvec[2]; // pval for G ~ L|T
		pval4 = pvec[3]; // pval for L ind T|G

		pvec.clear();
		gresid.clear();
		gpred.clear();
		gsl_matrix_free (Lm);
		gsl_vector_free (Gm);
		gsl_vector_free (Tm);
		gsl_vector_free (Gp);
		gsl_matrix_free(X);

	delete [] designmat;
	delete [] phenovec;

	PutRNGstate();
	LL.clear();

} // End citconlog2
} // End extern c wrapper
