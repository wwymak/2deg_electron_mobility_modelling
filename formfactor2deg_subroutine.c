#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define ZERO_THICKNESS_WITH_ITSELF 0
#define ZERO_THICKNESS_WITH_POINT_CHARGE 1
#define COSINE_WAVEFN_WITH_ITSELF  2
#define COSINE_WAVEFN_WITH_POINT_CHARGE 3
#define FANG_HOWARD_WAVEFN_WITH_ITSELF 4
#define FANG_HOWARD_WAVEFN_WITH_POINT_CHARGE 5
#define POINT_CHARGE_AND_WAVEFN_NUMERIC 6
#define WAVEFN_WITH_ITSELF_NUMERIC 7

#define PI 3.14159265



double  FormFactor2DEG(int, double, double, double, double, double  **, int *);




double  FormFactor2DEG(int WAVEFN_TYPE,
                       double  q,
					   double  param1,
					   double  param2,
					   double  param3,
					   double  **wvfn_val,
					   int *errcode_ptr
					   )

		  /* For zero thickness plane

		     if WAVEFN_TYPE = ZERO_THICKNESS_WITH_ITESELF

			 param1 is not used
		     param2 is not used
		     param3 is not used

		     if WAVEFN_TYPE = ZERO_THICKNESS_WITH_POINT_CHARGE

			 param1 is not used
		     param2 is not used
		     param3 is the distance of the charge from the plane



		     For heterointerface...

		     if WAVEFN_TYPE= FANG_HOWARD_WITH_ITSELF

		     param1 is the Fang Howard b paramater
		     param2 is the ratio k_insulator/k_semiconductor
		     param3 is not used

		     if WAVEFN_TYPE= FANG_HOWARD_WITH_POINT_CHARGE

		     param1 is the Fang Howard b paramater
		     param2 is the ratio k_insulator/k_semiconductor
			 param3 is the z co-ordinate of the point charge (0,0,z)



   		    For quantum wells...

		    if WAVEFN_TYPE = COSINE_WAVEFN_WITH_ITSELF

		    param1 is the QW width
		    param2 is not used
		    param3 is not used

		    if WAVEFN_TYPE = COSINE_WAVEFN_WITH_POINT_CHARGE

			param1 is the QW width
			param2 is not used
		    param3 is the co-ordinate of the point charge, assuming  z=0 is the center of the well

		  */




{ double  fq, fq1,fq2, qL, qbyb, dielectric_ratio, fourPIsq;
  double  a0,a1,a2,bz,four_pi_sq_by_l_sq;


 /* calling  program has to ensure that the units of q and lengths are inverses of each other
    This is possible because the calculation can be done in terms of qL, qz0 etc., which are
	dimensionless. The return value itself is dimensionless.

 */

  q=fabs(q); /* only the magnitude of q should matter - scattering  is isotropic */

  switch(WAVEFN_TYPE)
   { case ZERO_THICKNESS_WITH_ITSELF:
          fq= 1.0;
		  break;

	 case ZERO_THICKNESS_WITH_POINT_CHARGE:
          /* param3 is treated as the distance of the point charge from the plane */
          fq=exp(-1.0*fabs(param3)*q);
          break;


	 case COSINE_WAVEFN_WITH_POINT_CHARGE :
		   /* We need two parameters,
		      param1 is the quantum well width
		      param2 is not used
			  param3 is the co-ordinate of the point charge, assume that
			  z=0 is the center of the well.

		      Consider two cases, fabs(z) > L/2 & fabs(z) < L/2
		   */
              param1=fabs(param1);

			  qL=q*param1;
			  fourPIsq=4.0*PI*PI;


           if(fabs(param3) > 0.5*param1)  /* The charge is outside the QW */
			   {
				 fq=exp(-q*fabs(param3))*fourPIsq/(fourPIsq + qL*qL);
				 qL *=0.5;
				 fq = (qL < 1e-6) ? fq : fq*sinh(qL)/qL;

				}

		   if(fabs(param3) <= 0.5*param1)  /* The charge is inside the QW */
			   {

                 fq=exp(-0.5*qL)*sinh(q*param3)*2.0*qL/(fourPIsq + qL*qL);
                 qL *=0.5;
				 fq = (qL < 1e-6) ?  1.0 - fq : (1- exp(-1.0*qL)*cosh(q*param3))/qL - fq;


				}

              break;

	 case COSINE_WAVEFN_WITH_ITSELF:
		  /* param1 is the quantum well width

		   */
		  qL = q*param1;
          fourPIsq=4.0*PI*PI;

		  /* The analytic expression has the correct limit as q-> 0, but it has a zero/zero form.
		   Part of the expression need to be expanded to get rid of this non-analytic point, so that
		   the numerical calculation does not encounter a problem. The analytic expression can be found
		   in Jauho and Smith, PRB 47, 4420 (1993)
		   */

		  if(qL < 1e-4)
			{fq= 1.0 + qL/(fourPIsq + qL*qL) - 0.23201214969*qL + 0.032672741512*qL*qL;
			 }else{fq = 2/qL + qL/(fourPIsq + qL*qL) - (2/qL)*pow(1 + qL*qL/fourPIsq, -2)*(1 - exp(-1.0*qL))/qL;
				  }

		  break;


	 case FANG_HOWARD_WAVEFN_WITH_ITSELF:
		  /* param1 is the Fang Howard b paramater
		     param2 is the ratio k_insulator/k_semiconductor
		  */
		  qbyb=q/param1;
		  dielectric_ratio=param2;

          fq = (1+dielectric_ratio)*pow(1+qbyb,-3.0)*(0.5+ 0.5625*qbyb + 0.1865*qbyb*qbyb) + 0.5*(1-dielectric_ratio)*pow(1+qbyb,-6.0);

          /* printf("q= %e \t b= %e \t FANG HOWARD fq = %e\n",q,param1,fq); */
          break;



     case FANG_HOWARD_WAVEFN_WITH_POINT_CHARGE:
		  /* param1 is the Fang Howard b paramater
		     param2 is the ratio k_insulator/k_semiconductor
			 param3 is the z co-ordinate of the point charge (0,0,z)

		     Using calculation given in Ando, Fowler and Stern RMP 54,437 (1982) p448,499
		  */

		  qbyb=q/param1;
		  dielectric_ratio=param2;
		  bz=param1*param3;

          if(param3 < 0.0){fq=pow(1 + qbyb,-3)*exp(q*param3); /* z_point_charge < 0 */
			              } else{ /* z_point_charge > 0 */
							     if(fabs(1-qbyb) < 1e-3){ /* the q=b case */
							                               fq1=0.0625*(1+dielectric_ratio)*(1 + 2*bz + 2*bz*bz + 1.33333*bz*bz*bz)*exp(-1.0*bz);
							                               fq2=0.5*(1-dielectric_ratio)*pow(1+qbyb,-3)*exp(-1.0*q*param3);
														   fq=fq1+fq2;

						                                  }else{ /* q neq b */
															    a0=2.0*qbyb*(3+qbyb*qbyb)*pow(1+qbyb,-3);
																a1=4.0*q*(1-qbyb)*pow(1+qbyb,-2);
																a2=param1*q*pow(1-qbyb,2)/(1+qbyb);
																fq1=0.5*(1+dielectric_ratio)*pow(1-qbyb,-3)*(exp(-1.0*q*param3) - (a0+a1*param3+a2*param3*param3)*exp(-1.0*bz));
																fq2=0.5*(1-dielectric_ratio)*pow(1+qbyb,-3)*exp(-1.0*q*param3);
																fq=fq1+fq2;
																}


							     }



          break;

     case POINT_CHARGE_AND_WAVEFN_NUMERIC:  /* TO BE WRITTEN */
		  /* param1 is the z co-ordinate of the point charge */
		  fq=1.0;
          break;

	 case WAVEFN_WITH_ITSELF_NUMERIC:     /* TO BE WRITTEN */
          fq=1.0;
	      break;

     default:                      /* !!!!! SHOULD NOT GET HERE !!!!  */
		   fq=1.0;
          break;


	}


 return fq;

}
