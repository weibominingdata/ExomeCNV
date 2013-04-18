#pragma once

#include <limits> // for NAN

#define BINOMIAL  1
#define POISSON   2
#define GAUSSIAN  3
#define GAMMA     4
#define NB        5
#define ZINB      6

/* Link */

#define LOGIT     1
#define LOG       2
#define IDENTITY  3
#define INVERSE   4


/* Constant  */

#define n_max (100)
#define ML_POSINF	std::numeric_limits<double>::infinity()
#define ML_NEGINF	-1*std::numeric_limits<double>::infinity()
#define ML_NAN		std::numeric_limits<double>::quiet_NaN()



class MathTools
{
public:
	MathTools(void);
	~MathTools(void);
public:
	static int	wcenter(double*, int, double*, int, double*);
	static int	wresid(double*,  int, double*, double*, double*, double*);
	static double wssq(double *y, int n, double *weights);
	static double varfun(int, double, double);
	static int muvalid(int, double);
	static double linkfun(int, double);
	static double invlink(int, double);
	static double dlink(int, double);
	static void	initialize(int, double*, double*, int, double*);


	////////////// version with prior ////////////////////////////////
	/* Fit a base model */

	static int glmFit(int* familyR, int* linkR, int* dims, int* nIter,
           double *y, double *prior, double *offset, 
           double *X, double *nTotal_binom, double *convR, 
           int *rank, double *Xb, double *fitted, double *resid, 
           double *weights, double *phi, int* trace, 
           double *scale, int *df_resid, double* beta) ;

	/*  log likelihood of Poisson */
	static double loglik_Poisson(int N, double* mu, double* y, double* w);

	/* log likelihood of negative binomial */

	static double loglik_NB(int N, double phi, double* mu, double* y, double* w);

	/* Score test for additional terms */

	static void glm_score_test(int* dims, double *Z, double *resid, 
						double *weights, double *Xb, double* scaleR,
						double *chi2, int *df);

	/* score and infor for solving MLE of phi */

	static void score_info(int N, double theta, double* mu, double* y, double* w, 
                double* score, double* info);

	/* MLE of phi */

	static int phi_ml(double* y, double* mu, int N, double* weights, 
            int limit, double eps, double* phi, int initPhi, int trace);

	/* glmNB */

	static int glmNB(int *dims, int *nIter, double *y, double *prior, 
          int *linkR, double *offset, double *X, double *convR, 
          int *rank, double *Xb, double *fitted, double *resid, 
          double *weights, double *phi, double *scale, 
          int *df_resid, int* family, double *twologlik, 
          double *scoreTestP, int *trace, double *beta);


	////////////// version without prior /////////////////////////
		static int glmFit(bool useLess, int* familyR, int* linkR, int* dims, int* nIter,
			   double *y, double *offset, double *z,
			   double *X, double *nTotal_binom, double *convR, 
			   int *rank, double *Xb, double *fitted, double *resid, 
			   double *weights, double *phi, int* trace, 
			   double *scale, int *df_resid, double* beta);

	/*  log likelihood of Poisson */
	static double loglik_Poisson(int N, double* mu, double* y);

	/* log likelihood of negative binomial */

	static double loglik_NB(int N, double phi, double* mu, double* y);


	static void score_info(int N, double theta, double* mu, double *y, 
					double* score, double* info);

	/* MLE of phi */

	static int phi_ml(double* y, double* mu, int N,  
			   int limit, double eps, double* phi, int initPhi, int trace);

	/* glmNB */

	static int glmNB(bool useLess,int *dims, int *nIter, double *y, double *z, 
			  int *linkR, double *offset, double *X, double *convR, 
			  int *rank, double *Xb, double *fitted, double *resid, 
			  double *weights, double *phi, double *scale, 
			  int *df_resid, int* family, double *twologlik, 
			  double *scoreTestP, int *trace, double *beta);

///////////////////////////////////////////////////////////////////////////////////


	static double lgammafn(double x);
	static void dpsifn(double x, int n, int kode, int m, double *ans, int *nz, int *ierr);
	static int Rf_i1mach(int i);
	static double Rf_d1mach(int i);
	static double digamma(double x);
	static double trigamma(double x);
	static double getCumulativeDensityNormal(const double x);
	static double pnorm(const double x,const double mean,const double stddev, const bool low_tail, const bool log);
	static void max(double *v, int nl, int nh, double *val, int *idx);
	static double logsumexp (double* v, const int RN);
	static double dnorm(double x, double mu, double sigma);
	static double fac(int n);
	static double dbinom(int x, int size, double p);

///////////////////////////////////////////////////////////////////////////////////// beta-binomial function
	static double betaBinomial(double n, double k, double pi, double theta);
	static double lchoose(double n, double k);
	static double lfastchoose(double n, double k);
	static double lfastchoose2(double n, double k, int *s_choose);
	static double lbeta(double a, double b);
	static double log1p(double x);
	static int chebyshev_init(double *dos, int nos, double eta);
	static double chebyshev_eval(double x, const double *a, const int n);
	static double lgammacor(double x);
	static double gammafn(double x);
	static void gammalims(double *xmin, double *xmax);
	static double stirlerr(double n);




private:
 	static const double _PI;
	static const double _LOG10_2;
	
};
