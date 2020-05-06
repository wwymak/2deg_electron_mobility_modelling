/* Routines to calculate the real and imaginary parts of the non-interacting susceptibility
(density-density response function) of a homogeneous two-dimensional electron/hole gas,
as a function of wavevector and frequency, at zero or finite temperature */


/*

************************************************************************************************

 The subroutine calculates the REAL and IMAGINARY parts of the Lindhard susceptibilty of a
 2D electron/hole gas at zero or finite temperature and zero or finite frequencies.

 To compute T=0 values, I use the formulation of Giuliani & Vignale
 To compute finite temperature values , I use the Maldague integration method over all T=0 values,
 the integral is calculated with an adaptive step-size method - that essentially treats the
 integral as a differential equation and evolves it from intial to final value of the independent
 variable.



 The range of integration is from [0.001 - 3E_f], this can be changed by changing the lowlim and
 ulim variables.

 Usage:
 double RE_chi2DEG(double q_by_kF, double w_by_wF, double hbar_wF_by_kBT, double mu_by_hbar_wF,
                   int option, int *errorcode_ptr)

 double IM_chi2DEG(double q_by_kF, double w_by_wF, double hbar_wF_by_kBT, double mu_by_hbar_wF,
                   int option, int *errorcode_ptr)

 inputs and return values :

 qF, wF etc denotes values at Fermi energy.

 q_by_kF = q/kF

 w_by_wF = omega/omegaF

 hbar_wF_by_kBT = the temperature of the system, to be entered in units of Fermi energy EF/kBT

 mu_by_hbar_wF = allows accounting for the variation of chemical potential with temperature,
                 if more accuracte values are necessary. Set this = 1, to start.

 option : DYNAMIC_ZERO_T   (finite w, T=0)
          STATIC_ZERO_T    (w=0, T=0)
          DYNAMIC_FINITE_T (finite w, finite T)
          STATIC_FINITE_T  (w=0, finite T)

 Simpler expressions exist for the static cases, including this option allows using the simpler
 formula for static cases. Using w=0, with DYNAMIC option or very large EF/kBT with FINITE_T option
 will produce almost same results as the STATIC or ZERO_T options. But this will use
 more complex formulas to do the computation and take longer.

 errorcodes : A pointer to an integer is sent to the routine. The routine will additively set the
              following errors. Check the bits of the (integer) errorcode on return. It is possible
			  that more than one error occurred.

			  TO BE WRITTEN




 The returned value is scaled by the Density of states at Fermi level.
 Return value = -\chi(q,w,T)/N(Ef)

*************************************************************************************************

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DYNAMIC_ZERO_T    0
#define STATIC_ZERO_T 100
#define DYNAMIC_FINITE_T  1
#define STATIC_FINITE_T 101

#define ACCURACY_LEVEL 1e-6
#define MAXSTEPS 1000
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON5 1.89e-4   /* the value of ERRCON5 = (5.0/SAFETY)^(1/PGROW) */
#define TINY 1.0e-15

#define CONVERGENCE_REACHED 0
#define CONVERGENCE_NOT_REACHED 1
#define NO_ERROR_OCCURRED 0
#define Q_MAGNITUDE_TOO_SMALL 2
#define W_BY_Q_LARGER_THAN_VF 4
#define UNKNOWN_CASE_ENCOUNTERED 8


double RE_chi_2DEG(double,double,double,double,int,int*);
double IM_chi_2DEG(double,double,double,double,int,int*);

void MALDAGUE_rkqs(double*,double,double*,double,double,double,double*,double*,double,double,double,double,double(*)(double,double,double,double,double));
void MALDAGUE_odeint(double*,double,double,double,double,double,int,int*,int*,double,double,double,double,double(*)(double,double,double,double,double));

double real_integrand_static(double, double, double, double, double);
double real_integrand_dynamic(double,double,double,double,double);
double imag_integrand_dynamic(double,double,double,double,double);

void nrerror(char error_text[]);


double RE_chi_2DEG(double q_by_kF, double w_by_wF, double eF_by_kT, double accuracy, int option, int *error_ptr)
{
	/* Returns the real part of the density-density response function of a non-interacting
	2D Fermi gas (i.e. the Lindhard function) */

	double nu_plus, nu_minus, nu_plus2, nu_minus2, x;
	double mu_by_kT; // chemical potential in units of thermal energy
	double upper_cutoff, lower_cutoff, largenum;
	double re_chi; // the answer!

	// variables for adaptive step size integration
	double yval,x1,x2,eps,h1,hmin;
	int nok,nbad;

	mu_by_kT = eF_by_kT; // will be approximately eF_by_kT at low enough T
	if(eF_by_kT < 700.0) {mu_by_kT = log(exp(eF_by_kT)-1.0);}
	mu_by_kT = eF_by_kT;


	switch(option)
	{

	 case STATIC_ZERO_T:

		re_chi = -1.0;
		if(q_by_kF > 2.0) {re_chi += sqrt(1.0-4.0/(q_by_kF*q_by_kF));}

		// this case should be error-proof ...
		*error_ptr = *error_ptr + NO_ERROR_OCCURRED;

		break;


	 case DYNAMIC_ZERO_T:

		nu_plus = 0.5*(w_by_wF / q_by_kF + q_by_kF);
		nu_minus = nu_plus - q_by_kF;
		nu_plus2 = nu_plus*nu_plus;
		nu_minus2 = nu_minus*nu_minus;

		re_chi = 0.0;
		if(nu_plus2 > 1.0){re_chi += sqrt(nu_plus2 - 1.0);}
		if(nu_minus2 > 1.0) { re_chi -= (nu_minus>0.0) ? sqrt(nu_minus2-1.0) : -1.0*sqrt(nu_minus2-1.0); }
		re_chi = re_chi/q_by_kF -1.0;

		*error_ptr = *error_ptr + NO_ERROR_OCCURRED;

		break;

	 case STATIC_FINITE_T:


		nu_plus = 0.5*q_by_kF;

		re_chi=0.0;

		largenum = sqrt(1.0/accuracy); // want to look for where the cosh in the denominator exceeds this value
		upper_cutoff = (2.0*acosh(largenum) + mu_by_kT)/eF_by_kT;
		upper_cutoff = sqrt(upper_cutoff); // no point integrating above this value
		lower_cutoff = (-2.0*acosh(largenum) + mu_by_kT)/eF_by_kT;
		lower_cutoff = (lower_cutoff > 0.0) ? sqrt(lower_cutoff) : 0.0;


		if(q_by_kF > lower_cutoff)
		{
			/* Brute force integration*/
			/*for(x=lower_cutoff;x<=0.5*q_by_kF;x+=(0.5*q_by_kF-lower_cutoff)/1000.0)
			{re_chi+=real_integrand_static(x,nu_plus,nu_plus,eF_by_kT,mu_by_kT);}
			re_chi = re_chi * (0.5*q_by_kF-lower_cutoff)/1000.0;
			re_chi = re_chi * eF_by_kT / q_by_kF;*/

		// Adaptive step size integration
		yval =0.0;
		x1 = lower_cutoff;
		x2 = 0.5*q_by_kF;
		eps = accuracy;
		h1 = (0.5*q_by_kF-lower_cutoff)/100.0; // suggested initial step size
		hmin = (0.5*q_by_kF-lower_cutoff)/1.0e5 + TINY;

		MALDAGUE_odeint(&yval, x1, x2,eps, h1,hmin, MAXSTEPS,&nok,&nbad,nu_minus,nu_plus,eF_by_kT,mu_by_kT,real_integrand_static);

		printf("\nodeint finished normally : y = %e : nok = %d : nbad = %d\n",yval,nok,nbad);

		re_chi = eF_by_kT * yval / q_by_kF;

		}

		re_chi -= 0.5*(1.0 + tanh(0.5*mu_by_kT));

		break;


	 case DYNAMIC_FINITE_T:

		nu_plus = 0.5*(w_by_wF / q_by_kF + q_by_kF);
		nu_minus = nu_plus - q_by_kF;
		nu_plus2 = nu_plus*nu_plus;
		nu_minus2 = nu_minus*nu_minus;

		re_chi = 0.0;

		largenum = sqrt(1.0/accuracy); // want to look for where the cosh in the denominator exceeds this value
		//largenum = 1.0e4;
		upper_cutoff = (2.0*acosh(largenum) + mu_by_kT)/eF_by_kT;
		upper_cutoff = sqrt(upper_cutoff); // no point integrating above this value
		lower_cutoff = (-2.0*acosh(largenum) + mu_by_kT)/eF_by_kT;
		lower_cutoff = (lower_cutoff > 0.0) ? sqrt(lower_cutoff) : 0.0;

		if(fabs(nu_minus) > lower_cutoff)
		{

		// Adaptive step size integration
		yval =0.0;
		x1 = lower_cutoff;
		x2 = fmin(fabs(nu_minus),upper_cutoff);
		eps = accuracy;
		h1 = (x2-x1)/100.0; // suggested initial step size
		hmin = (x2-x1)/1.0e5 + TINY;

		MALDAGUE_odeint(&yval, x1, x2,eps, h1,hmin, MAXSTEPS,&nok,&nbad,nu_minus,nu_plus,eF_by_kT,mu_by_kT,real_integrand_dynamic);

		printf("\nodeint finished normally : y = %e : nok = %d : nbad = %d\n",yval,nok,nbad);

		re_chi += 0.5 * eF_by_kT * yval / q_by_kF;

		}

		if((fabs(nu_plus) > lower_cutoff)&&(nu_minus<upper_cutoff))
		{

		// Adaptive step size integration
		yval =0.0;
		x1 = fmax(fabs(nu_minus),lower_cutoff);
		x2 = fmin(nu_plus,upper_cutoff);
		eps = accuracy;
		h1 = (x2-x1)/100.0; // suggested initial step size
		hmin = (x2-x1)/1.0e5 + TINY;

		MALDAGUE_odeint(&yval, x1, x2,eps, h1,hmin, MAXSTEPS,&nok,&nbad,nu_minus,nu_plus,eF_by_kT,mu_by_kT,real_integrand_dynamic);

		printf("\nodeint finished normally : y = %e : nok = %d : nbad = %d\n",yval,nok,nbad);

		re_chi += 0.5 * eF_by_kT * yval / q_by_kF;

		}

		re_chi -= 0.5*(1.0 + tanh(0.5*mu_by_kT));



		break;



	 default: re_chi = -1.0;

		*error_ptr = *error_ptr + UNKNOWN_CASE_ENCOUNTERED;


	}

	return (re_chi);

}

double IM_chi_2DEG(double q_by_kF, double w_by_wF, double eF_by_kT, double accuracy, int option, int *error_ptr)
{

	/* Returns the imaginary part of the density-density response function of a non-interacting
	2D Fermi gas (i.e. the Lindhard function) */

	double nu_plus, nu_minus, nu_plus2, nu_minus2;
	double mu_by_kT; // chemical potential in units of thermal energy
	double im_chi; // the answer!
	double upper_cutoff, lower_cutoff, largenum;

	// variables for adaptive step size integration
	double yval,x1,x2,eps,h1,hmin;
	int nok,nbad;

	mu_by_kT = eF_by_kT; // will be approximately eF_by_kT at low enough T
	if(eF_by_kT < 700.0) {mu_by_kT = log(exp(eF_by_kT)-1.0);}

	switch(option)
	{

	 case STATIC_ZERO_T:

		im_chi = 0.0;

		// this case should be error-proof ...
		*error_ptr = *error_ptr + NO_ERROR_OCCURRED;

		break;


	 case DYNAMIC_ZERO_T:

		nu_plus = 0.5*(w_by_wF / q_by_kF + q_by_kF);
		nu_minus = nu_plus - q_by_kF;
		nu_plus2 = nu_plus*nu_plus;
		nu_minus2 = nu_minus*nu_minus;

		im_chi = 0.0;
		if(nu_plus2 < 1.0){im_chi += sqrt(1.0 - nu_plus2);}
		if(nu_minus2 < 1.0) { im_chi -= sqrt(1.0-nu_minus2);}
		im_chi = im_chi/q_by_kF;

		*error_ptr = *error_ptr + NO_ERROR_OCCURRED;

		break;



	 case STATIC_FINITE_T:

		im_chi = 0.0;

		// this case should be error-proof ...
		*error_ptr = *error_ptr + NO_ERROR_OCCURRED;

		break;



	case DYNAMIC_FINITE_T:

		nu_plus = 0.5*(w_by_wF / q_by_kF + q_by_kF);
		nu_minus = nu_plus - q_by_kF;
		nu_plus2 = nu_plus*nu_plus;
		nu_minus2 = nu_minus*nu_minus;

		im_chi = 0.0;

		largenum = sqrt(1.0/accuracy); // want to look for where the cosh in the denominator exceeds this value
		largenum = 1.0e4;
		upper_cutoff = (2.0*acosh(largenum) + mu_by_kT)/eF_by_kT;
		upper_cutoff = sqrt(upper_cutoff); // no point integrating above this value
		lower_cutoff = (-2.0*acosh(largenum) + mu_by_kT)/eF_by_kT;
		lower_cutoff = (lower_cutoff > 0.0) ? sqrt(lower_cutoff) : 0.0;


		if(fabs(nu_plus) < upper_cutoff)
		{

		// Adaptive step size integration
		yval =0.0;
		x1 = fmax(lower_cutoff,nu_plus);
		x2 = upper_cutoff;
		eps = accuracy;
		h1 = (x2-x1)/100.0; // suggested initial step size
		hmin = (x2-x1)/1.0e6 + TINY;

		MALDAGUE_odeint(&yval, x1, x2,eps, h1,hmin, MAXSTEPS,&nok,&nbad,nu_minus,nu_plus,eF_by_kT,mu_by_kT,imag_integrand_dynamic);

		printf("\nodeint finished normally : y = %e : nok = %d : nbad = %d\n",yval,nok,nbad);

		im_chi -= 0.5 * eF_by_kT * yval / q_by_kF;

		}

		if((fabs(nu_plus) > lower_cutoff)&&(nu_minus<upper_cutoff))
		{

		// Adaptive step size integration
		yval =0.0;
		x2 = fmax(fabs(nu_minus),lower_cutoff);
		x1 = fmin(nu_plus,upper_cutoff);
		eps = accuracy;
		h1 = (x2-x1)/100.0; // suggested initial step size
		hmin = (x2-x1)/1.0e6 + TINY;

		MALDAGUE_odeint(&yval, x1, x2,eps, h1,hmin, MAXSTEPS,&nok,&nbad,nu_minus,nu_plus,eF_by_kT,mu_by_kT,imag_integrand_dynamic);

		printf("\nodeint finished normally : y = %e : nok = %d : nbad = %d\n",yval,nok,nbad);

		im_chi += 0.5 * eF_by_kT * yval / q_by_kF;

		}

		break;



	 default: im_chi = 0.0;

		*error_ptr = *error_ptr + UNKNOWN_CASE_ENCOUNTERED;


	}

	return (im_chi);


}



void MALDAGUE_rkqs(double *y, double dydx, double *x, double htry, double eps, double yscal, double *hdid, double *hnext, double nu_minus, double nu_plus, double eF_by_kT, double mu_by_kT, double (*MALDAGUE_INTEGRAND)(double,double,double,double,double))
{
	/* Takes one integration step, starting at x, with current integrated total y and integrand dydx at x;
	tries to step by htry, with the change in y calculated using fifth order RK, and its error using diff between
	fourth and fifth order RK methods;
	if the error is too big the step h is reduced, until the error is small enough */

	static double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
                c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
		        dc5=-277.0/14336.0;
	double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
         dc4=c4-13525.0/55296.0,dc6=c6-0.25;    /*  Runge-Kutta Cash-Carp paramters
		                                            from Numerical recipes...
		                                         */
	double ak2, ak3, ak4, ak5, ak6;
	int i;
	double errmax, h, htemp, xnew, yerr, ytemp;

	h = htry;
	for (;;) // infinite loop of decreasing step size, exitted when error small enough
	{
		// Fifth-order Runge-Kutta needs five evaluations of the integrand in the interval x to x+h
		ak2 = (*MALDAGUE_INTEGRAND)((*x+a2*h),nu_minus,nu_plus,eF_by_kT,mu_by_kT);
		//printf("ak2 = %le\n",ak2);
		ak3 = (*MALDAGUE_INTEGRAND)((*x+a3*h),nu_minus,nu_plus,eF_by_kT,mu_by_kT);
		ak4 = (*MALDAGUE_INTEGRAND)((*x+a4*h),nu_minus,nu_plus,eF_by_kT,mu_by_kT);
		ak5 = (*MALDAGUE_INTEGRAND)((*x+a5*h),nu_minus,nu_plus,eF_by_kT,mu_by_kT);
		ak6 = (*MALDAGUE_INTEGRAND)((*x+a6*h),nu_minus,nu_plus,eF_by_kT,mu_by_kT);

		// Estimate new y
		ytemp = *y + h*(c1*dydx + c3*ak3 + c4*ak4 + c6*ak6); /* Accumulate increments with proper weights */
		//printf("estimated ystep is %le\n",ytemp-(*y));

		// Estimate the error in the change in y using diff between fourth and fifth order methods
		yerr=h*(dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6);
		// Compare error with yscal and see if relative error is too big
		errmax = fabs(yerr/yscal); errmax /= eps;
		/* printf("rkqs : errmax = %e\n",errmax); */

		if(errmax <= 1.0) break;  /* step succeeded, compute size of next step
	                               after coming out of the infinite loop */

	    htemp=SAFETY*h*pow(errmax,PSHRNK); /* need to reduce stepsize, according to NR prescription */

	            /* h can be +ve or -ve depending on whether the integration
	               runs from a low value to high value or not. But in one
				   step don't reduce by more than a factor of 10 */

	   if(h>0.0){ h = (htemp < 0.1*h)? 0.1*h : htemp;}
	   if(h<0.0){ h = (htemp > 0.1*h)? 0.1*h : htemp;}

	   xnew = (*x) + h; // check if the new step size is too small
	   if(xnew == *x) nrerror("stepsize underflow in rkqs");

	}

	// Now have made the step so can return results to calling routine and suggest increased step
	*hnext = (errmax > ERRCON5)? SAFETY*h*pow(errmax,PGROW) : 1.5*h; /* not more than a
	                                                                factor of 5 increase */

    (*x) += (*hdid =h);
    (*y) = ytemp;

	return;
}



void MALDAGUE_odeint(double *yval,double x1, double x2, double eps, double h1, double hmin, int maxsteps, int *nok, int *nbad, double nu_minus, double nu_plus, double eF_by_kT, double mu_by_kT, double(*MALDAGUE_INTEGRAND)(double,double,double,double,double))
{
	/* Integrates the function MALDAGUE_INTEGRAND(x,nu_minus,nu_plus,eF/kT,mu/kT) between x1 and x2,
	using and adaptive step size routine, with relative error eps;
	yval is given as the inital value and will be changed to add on the integrated total */

	int i, nstep;
	double x,hnext,hdid,h;
	double yscal,dydx,y;

	x = x1;
	h=(x2 > x1) ? fabs(h1) : -fabs(h1); // make sure direction of step is right
	*nok = (*nbad) = 0;
	y = *yval;

	for(nstep=1; nstep<=maxsteps; nstep++)
	{
		dydx = (*MALDAGUE_INTEGRAND)(x,nu_minus,nu_plus,eF_by_kT,mu_by_kT);
		yscal = fabs(y) + fabs(dydx*h) + TINY; // erros will be compared with this

		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x; /* check if this step takes us out of the integration region;
												if the stepsize overshoots the upper limit decrease it */

		MALDAGUE_rkqs(&y,dydx,&x,h,eps,yscal,&hdid,&hnext,nu_minus,nu_plus,eF_by_kT,mu_by_kT,MALDAGUE_INTEGRAND); // call the routine to make the step
        *yval=y;

		if (hdid == h) ++(*nok); else ++(*nbad);

		if ((x-x2)*(x2-x1) >= 0.0)  /* upper limit reached ? */
			{
			   return;  /* normal exit to calling program */
		     }

		 printf("[%d]",nstep);
		 if (fabs(hnext) <= hmin) nrerror("Step size too small in ODEINT");
		 h=hnext;
	}

	return;
}

double real_integrand_static(double x, double nu_minus, double nu_plus, double eF_by_kT, double mu_by_kT)
{
	double cosh_factor;
	double numerator=0.0;
	double x2, nu_plus2;

	x2 = x*x; nu_plus2 = nu_plus*nu_plus;

	cosh_factor = 0.5*(x2*eF_by_kT - mu_by_kT);
	cosh_factor = cosh(cosh_factor);
	cosh_factor = cosh_factor*cosh_factor;

	if (nu_plus2 > x2)
	{numerator += sqrt(nu_plus2 - x2);}

	return(x*numerator / cosh_factor);
}



double real_integrand_dynamic(double x, double nu_minus, double nu_plus, double eF_by_kT, double mu_by_kT)
{
	double cosh_factor;
	double numerator=0.0;
	double x2, nu_plus2, nu_minus2;

	x2 = x*x;
	nu_plus2 = nu_plus*nu_plus; nu_minus2 = nu_minus*nu_minus;

	cosh_factor = 0.5*(x2*eF_by_kT - mu_by_kT);
	cosh_factor = cosh(cosh_factor);
	cosh_factor = cosh_factor*cosh_factor;

	if (nu_plus2 > x2)
	{numerator += sqrt(nu_plus2 - x2);}

	if (nu_minus2 > x2)
	{numerator -= ((nu_minus>0.0) ? sqrt(nu_minus2 - x2) : -1.0*sqrt(nu_minus2 - x2));}

	return(x*numerator / cosh_factor);
}


double imag_integrand_dynamic(double x, double nu_minus, double nu_plus, double eF_by_kT, double mu_by_kT)
{
	double cosh_factor;
	double numerator=0.0;
	double x2, nu_plus2, nu_minus2;

	x2 = x*x;
	nu_plus2 = nu_plus*nu_plus; nu_minus2 = nu_minus*nu_minus;

	cosh_factor = 0.5*(x2*eF_by_kT - mu_by_kT);
	cosh_factor = cosh(cosh_factor);
	cosh_factor = cosh_factor*cosh_factor;

	if (x2 > nu_plus2)
	{numerator -= sqrt(x2 - nu_plus2);}

	if (x2 > nu_minus2)
	{numerator += sqrt(x2 - nu_minus2);}

	return(x*numerator / cosh_factor);
}


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

