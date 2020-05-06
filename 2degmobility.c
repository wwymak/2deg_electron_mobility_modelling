#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


#define ZERO_THICKNESS_WITH_ITSELF 0
#define ZERO_THICKNESS_WITH_POINT_CHARGE 1
#define COSINE_WAVEFN_WITH_ITSELF  2
#define COSINE_WAVEFN_WITH_POINT_CHARGE 3
#define FANG_HOWARD_WAVEFN_WITH_ITSELF 4
#define FANG_HOWARD_WAVEFN_WITH_POINT_CHARGE 5
#define POINT_CHARGE_AND_WAVEFN_NUMERIC 6
#define WAVEFN_WITH_ITSELF_NUMERIC 7

#define PI 3.14159265


#define DYNAMIC_ZERO_T    0
#define STATIC_ZERO_T 100
#define DYNAMIC_FINITE_T  1
#define STATIC_FINITE_T 101

#define HETEROINTERFACE 0
#define QUANTUM_WELL    1


#define YES 1
#define NO 0
#define WRITE_VTHETA_TO_FILE 1
#define DO_NOT_WRITE_VTHETA_TO_FILE 0


extern double RE_chi_2DEG(double,double,double,double,int,int*);
extern double IM_chi_2DEG(double,double,double,double,int,int*);

extern double  FormFactor2DEG(int, double, double, double, double, double  **, int *);



extern int getline_from_file(FILE *, char *,char *);
extern int remove_comments_from_linestring(char *);
extern int getsection_from_file(FILE *, char *, char *, int *, char *);
extern int get_number_of_numbers_in_a_string(char *);
extern int search_string_and_update_value(char *, char *, char *);



double one_by_tau_for_a_slab_of_dopant(int, double,double , double, double, double, double, double, double, int,int,int,char *);
double one_by_tau_for_a_delta_doped_layer(int, double,double, double, double, double, double,double,int,int,int, char *);
double one_by_tau_for_interface_roughness(int, double,double, double, double, double, double, double, int, int,int,char *);

/* one_by_tau_for_alloy_scattering()  to be written */




double hbar,epsilon0,m0,electron_charge;  /* Values are set in main */



int main (void)
{
  double  q, qwell_width=20,fq1,fq2,fang_howard_b_in_nanometer;  /* unit of length is  nanometer */
  double  **wavefn, param1=0.0,param2=0.0,param3=0.0;
  double  k_Fermi, k_Fermi_in_nanometer, epsilon_qw_2DEG, q_epsilon_qw_2DEG,qtf,qtf_in_nanometer,
          qtf_by_2k_Fermi, density_of_states_at_Fermi_level, bare_coulomb_prefactor;

  double  q_by_kF, w_by_wF, hbar_wF_by_kBT,  mu_by_hbar_wF;

  double areal_density_of_delta_doping,z_co_ord_of_delta_doping_in_nanometer,one_by_tau_delta_doping;

  double modulation_doping_start_in_nanometer, modulation_doping_end_in_nanometer, modulation_doping_density_in_per_m3,
         slab_width_in_nanometer,one_by_tau_modulation_doping;

  double background_doping_start_in_nanometer, background_doping_end_in_nanometer, background_doping_density_in_per_m3,
         one_by_tau_background_doping;

  double roughness_amplitude_in_nanometer, roughness_correl_length_in_nanometer, one_by_tau_interface_roughness,
         N_depletion_charge_density;

  double one_by_tau_total;

  double ss2,sum,theta,dtheta,V_theta,V_theta_squared_sum,mobility,tmp;

  double  N2d, epsilonr, background_in_GaAs, background_in_AlGaAs, GaAs_cap_width_in_nanometer, distance_of_interface_from_surface_in_nanometer,
          surface_charge_density;

  double quantum_well_width_in_nanometer;

 //  double datapts[]={1.2e14,1.75e14,2.5e14,4e14,5.5e14,7e14,8.5e14,1e15,1.25e15,1.5e15,1.75e15,2e15,2.25e15,2.5e15,2.75e15,3.0e15,3.5e15,4.0e15,5.0e15};

  double *datapts;

  FILE *fpout, *fpin;

  double *ImpurityProfileArray, *DeltaDopingArray, *ModulationDopingArray;
  int    NLayersForImpurityProfile=0, NLayersForDeltaDoping=0, NLayersForModulationDoping=0;

  int i,ii,j,ndatapts,errcode=0.0;

  int OUTPUT_IN_SCALED_1E10CM2_AND_1E6CM2VS_UNITS ; /* if NO   output will be in SI units */


  int HETEROINTERFACE_OR_QUANTUM_WELL,
      INCLUDE_SCATTERING_FROM_SURFACE_STATES,
	  INCLUDE_SCATTERING_FROM_DELTA_DOPED_LAYER,
	  INCLUDE_SCATTERING_FROM_BACKGR,
	  INCLUDE_SCATTERING_FROM_MODULATION_DOPED_LAYER,
	  INCLUDE_SCATTERING_FROM_INTERFACE_ROUGHNESS;


  char datastring[65536], errmsg[2048], returnvaluestring[256], outputfilename[2048]; /* just make these big enough for the purpose.....don't bother input length check and realloc... */
  char *pch;
  int errorcode;

  hbar=1.0557e-34; m0=9.1e-31;electron_charge=1.601e-19;epsilon0=8.854e-12; /* All in SI units */

  //ndatapts=sizeof(datapts)/sizeof(double);




  /************************************
   The parameters of the structure
   ************************************/




/*******************************************BEGIN INPUT FILE READ SECTION *********************************/

 fpin=fopen("input_for_mobility.txt","r");

 /********************************************************************************************************/

 /*section for OuputSpecifications */

 OUTPUT_IN_SCALED_1E10CM2_AND_1E6CM2VS_UNITS = YES; /* if NO   output will be in SI units */
 outputfilename[0]='\0';


 datastring[0]='\0';
 errorcode=0;
 getsection_from_file(fpin, "OutputSpecifications", datastring, &errorcode, errmsg); printf("%s\n",errorcode? "NOT FOUND":"FOUND");
 if(!errorcode){

			    if(search_string_and_update_value("outputfile",datastring,returnvaluestring))
				  {
					 strcpy(outputfilename,returnvaluestring);
                  }


 			    if(search_string_and_update_value("outputunits",datastring,returnvaluestring))
				  {ii=strlen(returnvaluestring);
			       for(i=0;i<ii;i++){returnvaluestring[i]=toupper(returnvaluestring[i]);}
				   if( strcmp(returnvaluestring,"SI")==0)
					  {OUTPUT_IN_SCALED_1E10CM2_AND_1E6CM2VS_UNITS = NO;
				      }
                  }



			  }



 /*******************************************************************************************************/

 /* section for DeviceType */
 /* default values */

  HETEROINTERFACE_OR_QUANTUM_WELL = QUANTUM_WELL;
  quantum_well_width_in_nanometer=20.00;
  distance_of_interface_from_surface_in_nanometer = -310;  /* this is negative in our convention  */
                                                           /* z=0 is the first heterointerface    */
                                                           /* it may be the top surface of the QW */

 datastring[0]='\0';
 errorcode=0;
 getsection_from_file(fpin, "DeviceType", datastring, &errorcode, errmsg); printf("%s\n",errorcode? "NOT FOUND":"FOUND");
 if(!errorcode){

				if(search_string_and_update_value("quantum_well_width_in_nanometer",datastring,returnvaluestring))
				  {sscanf(returnvaluestring,"%lf",&quantum_well_width_in_nanometer);
				  }


				if(search_string_and_update_value("distance_of_interface_from_surface_in_nanometer",datastring,returnvaluestring))
                  {sscanf(returnvaluestring,"%lf",&distance_of_interface_from_surface_in_nanometer);
				  }


                if(search_string_and_update_value("devicetype",datastring,returnvaluestring))
				  {ii=strlen(returnvaluestring);
			       for(i=0;i<ii;i++){returnvaluestring[i]=toupper(returnvaluestring[i]);}
				   if( (strcmp(returnvaluestring,"HETEROINTERFACE")==0) || (strcmp(returnvaluestring,"HETEROSTRUCTURE")==0))
					  {HETEROINTERFACE_OR_QUANTUM_WELL = HETEROINTERFACE;
				      }else {HETEROINTERFACE_OR_QUANTUM_WELL = QUANTUM_WELL;}
				  }



               }

 /*******************************************************************************************************************/
 /* section for CarrierDensity */

 /* default - just the point 1e15/m^2 */

  ndatapts = 1;
  datapts=(double *)malloc(ndatapts*sizeof(double));
  datapts[0]=1e15;

 datastring[0]='\0';
 errorcode=0;
 getsection_from_file(fpin, "CarrierDensity", datastring, &errorcode, errmsg); printf("%s\n",errorcode? "NOT FOUND":"FOUND");

 if(!errorcode){
	             ndatapts=get_number_of_numbers_in_a_string(datastring);

				 if(ndatapts>0)
				   {
				     free(datapts);
				     datapts=(double *)malloc(ndatapts*sizeof(double));

					 	   i=0;             /* Hoping this is a safe way to read the string of values to the array */
						   pch = strtok (datastring," ,\t");
                           while (pch != NULL)
                                {
                                  sscanf(pch,"%lf",&datapts[i++]);
                                  pch = strtok (NULL, " ,\t");
                                 }

                 }

			   }


 /********************************************************************************************************************/

  /* section  for SurfaceCharge */
  /* default values */


  INCLUDE_SCATTERING_FROM_SURFACE_STATES = NO;

  surface_charge_density = 2.5e14;


 datastring[0]='\0';
 errorcode=0;
 getsection_from_file(fpin, "SurfaceCharge", datastring, &errorcode, errmsg); printf("%s\n",errorcode? "NOT FOUND":"FOUND");


 if(!errorcode){INCLUDE_SCATTERING_FROM_SURFACE_STATES = YES;
                 /* errorcode !=0 means the section is either commented out or has errors */

				if(search_string_and_update_value("surface_charge_density",datastring,returnvaluestring))
				  {sscanf(returnvaluestring,"%lf",&surface_charge_density);
				  }


			    /* distance of the interface from surface is already known so there
				   is nothing more needed

				*/

			   }


 /***************************************************************************************************/


  /* section  for ImpurityProfile or the background */
  /* default values */

  INCLUDE_SCATTERING_FROM_BACKGR = NO;


  datastring[0]='\0';
  errorcode=0;
  getsection_from_file(fpin, "ImpurityProfile", datastring, &errorcode, errmsg); printf("%s\n",errorcode? "NOT FOUND":"FOUND");

  if(!errorcode){

                 NLayersForImpurityProfile =  get_number_of_numbers_in_a_string(datastring);

				 if( ((NLayersForImpurityProfile%3) == 0) && (NLayersForImpurityProfile > 0) )
					{
					  INCLUDE_SCATTERING_FROM_BACKGR = YES;
					  /* The array must occur as start stop amount of impurity */

				      ImpurityProfileArray = (double *) malloc (NLayersForImpurityProfile*sizeof(double));
                      NLayersForImpurityProfile = NLayersForImpurityProfile/3 ;

					 	   i=0;             /* Hoping this is a safe way to read the string of values to the array */
						   pch = strtok (datastring," ,\t");
                           while (pch != NULL)
                                {
                                  sscanf(pch,"%lf",&ImpurityProfileArray[i++]);
                                  pch = strtok (NULL, " ,\t");
                                 }
					 }


				}



  /**********************************************************************************************/

  /* section for ModulationDoping is similar to background doping */


  INCLUDE_SCATTERING_FROM_MODULATION_DOPED_LAYER = NO;

  datastring[0]='\0';
  errorcode=0;
  getsection_from_file(fpin,"ModulationDoping", datastring, &errorcode, errmsg); printf("%s\n",errorcode? "NOT FOUND":"FOUND");


  if(!errorcode){

                 NLayersForModulationDoping =  get_number_of_numbers_in_a_string(datastring);

				 if( ((NLayersForModulationDoping%3) == 0)  && (NLayersForModulationDoping > 0))
					{
					  INCLUDE_SCATTERING_FROM_MODULATION_DOPED_LAYER = YES;
					  /* The array must occur as start : stop : amount of impurity */

				      ModulationDopingArray = (double *) malloc (NLayersForModulationDoping*sizeof(double));
				      NLayersForModulationDoping = NLayersForModulationDoping/3 ;




                      puts(datastring);

					 	   i=0;
						   pch = strtok (datastring," ,\t");
                           while (pch != NULL)
                                {
                                  sscanf(pch,"%lf",&ModulationDopingArray[i++]);
                                  pch = strtok (NULL, " ,\t");
                                 }



					 }


				}


  /***************************************************************************************************/

  INCLUDE_SCATTERING_FROM_DELTA_DOPED_LAYER=NO;

  datastring[0]='\0';
  errorcode=0;
  getsection_from_file(fpin, "DeltaDoping", datastring, &errorcode, errmsg);  printf("%s\n",errorcode? "NOT FOUND":"FOUND");

  if(!errorcode){

                 INCLUDE_SCATTERING_FROM_DELTA_DOPED_LAYER=YES;

                 NLayersForDeltaDoping =  get_number_of_numbers_in_a_string(datastring);

				 if(((NLayersForDeltaDoping%2) == 0)  && (NLayersForDeltaDoping>0))
					{
					  INCLUDE_SCATTERING_FROM_DELTA_DOPED_LAYER = YES;
					  /* The array must occur as location : impurity conc*/

				      DeltaDopingArray = (double *) malloc (NLayersForDeltaDoping*sizeof(double));
				      NLayersForDeltaDoping = NLayersForDeltaDoping/2 ;



                      puts(datastring);

					 	   i=0;
						   pch = strtok (datastring," ,\t");
                           while (pch != NULL)
                                {
                                  sscanf(pch,"%lf",&DeltaDopingArray[i++]);
                                  pch = strtok (NULL, " ,\t");
                                 }


					 }


				}



 /***************************************************************************************************/
  INCLUDE_SCATTERING_FROM_INTERFACE_ROUGHNESS = NO;

  roughness_amplitude_in_nanometer = 0.2;
  roughness_correl_length_in_nanometer=29;
  N_depletion_charge_density=1e13;        /* in SI : this should not make a difference compared to N2d */


  datastring[0]='\0';
  errorcode=0;
  getsection_from_file(fpin, "InterfaceRoughness", datastring, &errorcode, errmsg);  printf("%s\n",errorcode? "NOT FOUND":"FOUND");


  if(!errorcode){ INCLUDE_SCATTERING_FROM_INTERFACE_ROUGHNESS = YES;

				if(search_string_and_update_value("roughness_amplitude_in_nanometer",datastring,returnvaluestring))
				  {sscanf(returnvaluestring,"%lf",&roughness_amplitude_in_nanometer);
				  }


				if(search_string_and_update_value("roughness_correl_length_in_nanometer",datastring,returnvaluestring))
                  {sscanf(returnvaluestring,"%lf",&roughness_correl_length_in_nanometer);
				  }


				if(search_string_and_update_value("N_depletion_charge_density",datastring,returnvaluestring))
                  {sscanf(returnvaluestring,"%lf",&N_depletion_charge_density);
				  }


               }


  fclose(fpin);

 /**************************** READING INPUT FILE DONE ****************************************************/


 /******************************* Verify that things have been read correctly *****************************/

 printf("DeviceType = %s\n",HETEROINTERFACE_OR_QUANTUM_WELL==QUANTUM_WELL? "Quantum Well" : "Heterointerface");
 printf("distance_of_interface_from_surface_in_nanometer =%e\n",distance_of_interface_from_surface_in_nanometer);
 if(HETEROINTERFACE_OR_QUANTUM_WELL==QUANTUM_WELL)
   {printf("quantum_well_width_in_nanometer=%e\n",quantum_well_width_in_nanometer);
   }


 if(INCLUDE_SCATTERING_FROM_SURFACE_STATES)
   { printf("Surface at %.2f nm away: Charge = %e\n", distance_of_interface_from_surface_in_nanometer, surface_charge_density);
   }



  if(INCLUDE_SCATTERING_FROM_BACKGR)
	{

       printf("NLayersForImpurityProfile=%d\n",NLayersForImpurityProfile);

       for(i=0; i< NLayersForImpurityProfile; i++)
	      { printf("%e\t",ImpurityProfileArray[3*i]);
	        printf("%e\t",ImpurityProfileArray[3*i +1]);
	        printf("%e\n",ImpurityProfileArray[3*i +2]);
	      }

    }




  if(INCLUDE_SCATTERING_FROM_DELTA_DOPED_LAYER)
    {

       printf("NLayersForDeltaDoping=%d\n",NLayersForDeltaDoping);

       for(i=0; i< NLayersForDeltaDoping; i++)
	      { printf("%e\t",DeltaDopingArray[2*i]);
	        printf("%e\n",DeltaDopingArray[2*i +1]);
	      }

    }



  if(INCLUDE_SCATTERING_FROM_MODULATION_DOPED_LAYER)
	{

       printf("NLayersForModulationDoping=%d\n",NLayersForModulationDoping);

       for(i=0; i< NLayersForModulationDoping; i++)
	      { printf("%e\t",ModulationDopingArray[3*i]);
	        printf("%e\t",ModulationDopingArray[3*i +1]);
	        printf("%e\n",ModulationDopingArray[3*i +2]);
	      }

    }



  if(INCLUDE_SCATTERING_FROM_INTERFACE_ROUGHNESS)
    {

	printf("roughness_amplitude_in_nanometer=%e\n",roughness_amplitude_in_nanometer);
    printf("roughness_correl_length_in_nanometer=%e\n",roughness_correl_length_in_nanometer);
    printf("N_depletion_charge_density=%e\n",N_depletion_charge_density);


	}


  /********************************** Verification of parameters done..... **************************************/

 /****************************************************************************************************************/

  if(strlen(outputfilename))
	{fpout = fopen(outputfilename,"w");
    }else{fpout = fopen("output.dat","w");
	     }


  fprintf(fpout,"N2d%s\t",OUTPUT_IN_SCALED_1E10CM2_AND_1E6CM2VS_UNITS == YES? "(1e10percmsq)":"(permtsq)");

  fprintf(fpout,"%s",INCLUDE_SCATTERING_FROM_SURFACE_STATES == YES? "musurface\t" : "");

  fprintf(fpout,"%s",INCLUDE_SCATTERING_FROM_MODULATION_DOPED_LAYER == YES? "mumodulation\t" : "");

  fprintf(fpout,"%s",INCLUDE_SCATTERING_FROM_DELTA_DOPED_LAYER == YES? "mudelta\t" : "");

  fprintf(fpout,"%s",INCLUDE_SCATTERING_FROM_BACKGR == YES ? "mubg\t" : "");

  fprintf(fpout,"%s",INCLUDE_SCATTERING_FROM_INTERFACE_ROUGHNESS == YES? "muir\t" : "");

  fprintf(fpout,"mutotal%s\n",OUTPUT_IN_SCALED_1E10CM2_AND_1E6CM2VS_UNITS == YES? "(1e6cm^2/Vs)":"(m^2/Vs)");

  /*******************************************************************************/



  /* The calculation loops over each density from the datapts array */

  for(ii=0; ii < ndatapts ; ii++)
  {

   N2d=datapts[ii]; printf("[N2d = %e m^-2 ]\n",N2d);



   fprintf(fpout,"%e\t",OUTPUT_IN_SCALED_1E10CM2_AND_1E6CM2VS_UNITS == YES?  N2d/1e14 : N2d);




  epsilonr=12.8;                                                                   /* SI */
  density_of_states_at_Fermi_level = m0*0.067/(PI*hbar*hbar);                      /* SI */
  bare_coulomb_prefactor = electron_charge*electron_charge/(2*epsilon0*epsilonr);  /* SI */
  qtf= density_of_states_at_Fermi_level*bare_coulomb_prefactor;                    /* SI */
  qtf_in_nanometer= 1e-9*density_of_states_at_Fermi_level*bare_coulomb_prefactor;
  k_Fermi= sqrt(2*PI*N2d);                                                         /* SI */
  k_Fermi_in_nanometer = k_Fermi*1e-9;
  qtf_by_2k_Fermi=0.5*qtf/k_Fermi;




  one_by_tau_total=0.0;




  /**********************************************************************************************************************
   TREATED AS ONE  DELTA DOPED LAYER. COULD BE SURFACE STATES ALSO: All units in SI unless there is a suffix stating
   some other unit. Can have as many layers as needed.
   **********************************************************************************************************************/

  if(INCLUDE_SCATTERING_FROM_SURFACE_STATES == YES)
	{

      areal_density_of_delta_doping = surface_charge_density;
      z_co_ord_of_delta_doping_in_nanometer =distance_of_interface_from_surface_in_nanometer;


	 switch(HETEROINTERFACE_OR_QUANTUM_WELL)
	     {
	        case  HETEROINTERFACE :
            one_by_tau_delta_doping = one_by_tau_for_a_delta_doped_layer(HETEROINTERFACE_OR_QUANTUM_WELL,
			                                                             quantum_well_width_in_nanometer,
			                                                             z_co_ord_of_delta_doping_in_nanometer,
			                                                             areal_density_of_delta_doping,
                                                                         N2d, N_depletion_charge_density, m0*0.067, epsilonr,
																		 FANG_HOWARD_WAVEFN_WITH_POINT_CHARGE,
					               										 FANG_HOWARD_WAVEFN_WITH_ITSELF,
																		 DO_NOT_WRITE_VTHETA_TO_FILE,
																		 "vth_for_delta_doping.dat");

	               break;

	        case QUANTUM_WELL :

            one_by_tau_delta_doping = one_by_tau_for_a_delta_doped_layer(HETEROINTERFACE_OR_QUANTUM_WELL,
			                                                             quantum_well_width_in_nanometer,
			                                                             z_co_ord_of_delta_doping_in_nanometer,
			                                                             areal_density_of_delta_doping,
                                                                         N2d, N_depletion_charge_density, m0*0.067, epsilonr,
																		 COSINE_WAVEFN_WITH_POINT_CHARGE,
					               										 COSINE_WAVEFN_WITH_ITSELF,
																		 DO_NOT_WRITE_VTHETA_TO_FILE,
																		 "vth_for_delta_doping.dat");

	                break;


		    default:  /* generate some error or exit... */
		            break;

		 }

	  mobility = electron_charge/(m0*0.067*one_by_tau_delta_doping);
      /* printf("mu (surface charge) = %e cm^2V^-1s^-1\n", mobility*1e4);
	  */

      fprintf(fpout,"%e\t",OUTPUT_IN_SCALED_1E10CM2_AND_1E6CM2VS_UNITS == YES? mobility/100 : mobility);
      one_by_tau_total += one_by_tau_delta_doping;

    }




  /******************************************************************************************
   MODULATION DOPING LAYER  All units in SI unless there is a suffix stating some other unit
   ******************************************************************************************/

  if(INCLUDE_SCATTERING_FROM_MODULATION_DOPED_LAYER == YES)
	{

	 one_by_tau_modulation_doping = 0.00;

	 for(i=0; i <NLayersForModulationDoping; i++)

	  {
	   modulation_doping_start_in_nanometer = ModulationDopingArray[3*i];
       modulation_doping_end_in_nanometer   = ModulationDopingArray[3*i +1] ;
       modulation_doping_density_in_per_m3=   ModulationDopingArray[3*i +2];


	  switch(HETEROINTERFACE_OR_QUANTUM_WELL)
	        { case HETEROINTERFACE:

	          one_by_tau_modulation_doping += one_by_tau_for_a_slab_of_dopant(HETEROINTERFACE_OR_QUANTUM_WELL,
			                                                               quantum_well_width_in_nanometer,
	                                                                       modulation_doping_start_in_nanometer,
	                                                                       modulation_doping_end_in_nanometer,
                                                                           modulation_doping_density_in_per_m3,
																           N2d, N_depletion_charge_density, m0*0.067, epsilonr,
										                                   FANG_HOWARD_WAVEFN_WITH_POINT_CHARGE,
																           FANG_HOWARD_WAVEFN_WITH_ITSELF,
										                                   DO_NOT_WRITE_VTHETA_TO_FILE,
																		   "vth_for_modulation_doping.dat");


	               break;


		     case QUANTUM_WELL:

			 one_by_tau_modulation_doping += one_by_tau_for_a_slab_of_dopant(HETEROINTERFACE_OR_QUANTUM_WELL,
			                                                              quantum_well_width_in_nanometer,
	                                                                      modulation_doping_start_in_nanometer,
	                                                                      modulation_doping_end_in_nanometer,
                                                                          modulation_doping_density_in_per_m3,
																          N2d, N_depletion_charge_density, m0*0.067, epsilonr,
										                                  FANG_HOWARD_WAVEFN_WITH_POINT_CHARGE,
																          FANG_HOWARD_WAVEFN_WITH_ITSELF,
										                                  DO_NOT_WRITE_VTHETA_TO_FILE,
																		  "vth_for_modulation_doping.dat");



	              break;

			 default:  /* generate error or exit... */

				  break;

	         }


	     }


	 mobility = electron_charge/(m0*0.067*one_by_tau_modulation_doping);
     printf("mu (modulation doping) = %e cm^2V^-1s^-1 \n", mobility*1e4);


     fprintf(fpout,"%e\t",OUTPUT_IN_SCALED_1E10CM2_AND_1E6CM2VS_UNITS == YES? mobility/100 : mobility);
     one_by_tau_total += one_by_tau_modulation_doping;


   }


  /************************************************************************************************************/

   if(INCLUDE_SCATTERING_FROM_DELTA_DOPED_LAYER == YES)
     {

	    one_by_tau_delta_doping=0.00;

       for(i=0; i< NLayersForDeltaDoping; i++)
          {
	         z_co_ord_of_delta_doping_in_nanometer =DeltaDopingArray[2*i];
			 areal_density_of_delta_doping = DeltaDopingArray[2*i + 1];



	 switch(HETEROINTERFACE_OR_QUANTUM_WELL)
	      {
	        case  HETEROINTERFACE :
            one_by_tau_delta_doping += one_by_tau_for_a_delta_doped_layer(HETEROINTERFACE_OR_QUANTUM_WELL,
			                                                             quantum_well_width_in_nanometer,
			                                                             z_co_ord_of_delta_doping_in_nanometer,
			                                                             areal_density_of_delta_doping,
                                                                         N2d, N_depletion_charge_density, m0*0.067, epsilonr,
																		 FANG_HOWARD_WAVEFN_WITH_POINT_CHARGE,
					               										 FANG_HOWARD_WAVEFN_WITH_ITSELF,
																		 DO_NOT_WRITE_VTHETA_TO_FILE,
																		 "vth_for_delta_doping.dat");

	               break;

	        case QUANTUM_WELL :

            one_by_tau_delta_doping += one_by_tau_for_a_delta_doped_layer(HETEROINTERFACE_OR_QUANTUM_WELL,
			                                                             quantum_well_width_in_nanometer,
			                                                             z_co_ord_of_delta_doping_in_nanometer,
			                                                             areal_density_of_delta_doping,
                                                                         N2d, N_depletion_charge_density, m0*0.067, epsilonr,
																		 COSINE_WAVEFN_WITH_POINT_CHARGE,
					               										 COSINE_WAVEFN_WITH_ITSELF,
																		 DO_NOT_WRITE_VTHETA_TO_FILE,
																		 "vth_for_delta_doping.dat");

	                break;


		    default:  /* generate some error or exit... */
		            break;

		  }
	    }


   /******** Now the total contribution */
	  mobility = electron_charge/(m0*0.067*one_by_tau_delta_doping);
  //    printf("mu (delta doping) = %e cm^2V^-1s^-1\n", mobility*1e4);

      fprintf(fpout,"%e\t",OUTPUT_IN_SCALED_1E10CM2_AND_1E6CM2VS_UNITS == YES? mobility/100 : mobility);
      one_by_tau_total += one_by_tau_delta_doping;

   }

 /*************************************************************************************************************/





  /*****************************************************************
   BACKGROUND DOPING : treat it as a very low modulation doping :
   All units in SI unless there is a suffix stating some other unit
   Can be divided into as many slabs as needed to model different
   background in GaAs and AlGaAs etc.
   *****************************************************************/


   if(INCLUDE_SCATTERING_FROM_BACKGR == YES)
     {

	  one_by_tau_background_doping=0.00;

   for(i=0; i< NLayersForImpurityProfile; i++)
     {
       background_doping_start_in_nanometer = ImpurityProfileArray[3*i];
       background_doping_end_in_nanometer   = ImpurityProfileArray[3*i +1];
       background_doping_density_in_per_m3  = ImpurityProfileArray[3*i +2];




    switch(HETEROINTERFACE_OR_QUANTUM_WELL)
      {
		case HETEROINTERFACE:

	    one_by_tau_background_doping += one_by_tau_for_a_slab_of_dopant(HETEROINTERFACE_OR_QUANTUM_WELL,
                                                                     quantum_well_width_in_nanometer,
                                                                     background_doping_start_in_nanometer,
                                                                     background_doping_end_in_nanometer,
                                                                     background_doping_density_in_per_m3,
						        									 N2d, N_depletion_charge_density, m0*0.067, epsilonr,
										                             FANG_HOWARD_WAVEFN_WITH_POINT_CHARGE,
																     FANG_HOWARD_WAVEFN_WITH_ITSELF,
										                             DO_NOT_WRITE_VTHETA_TO_FILE,
																     "vth_for_bg1_doping.dat");


	    break;

	    case QUANTUM_WELL:

	    one_by_tau_background_doping += one_by_tau_for_a_slab_of_dopant(HETEROINTERFACE_OR_QUANTUM_WELL,
                                                                     quantum_well_width_in_nanometer,
                                                                     background_doping_start_in_nanometer,
                                                                     background_doping_end_in_nanometer,
                                                                     background_doping_density_in_per_m3,
						        									 N2d, N_depletion_charge_density, m0*0.067, epsilonr,
										                             COSINE_WAVEFN_WITH_POINT_CHARGE,
																     COSINE_WAVEFN_WITH_ITSELF,
										                             DO_NOT_WRITE_VTHETA_TO_FILE,
																     "vth_for_bg1_doping.dat");



	     break;


	     default :  /*  generate error or exit.... */


		 break;


	      }
        }

   /******** Now the total contribution */

   mobility = electron_charge/(m0*0.067*one_by_tau_background_doping);
   //printf("mu (background)= %e cm^2V^-1s^-1 \n", mobility*1e4);
   fprintf(fpout,"%e\t",OUTPUT_IN_SCALED_1E10CM2_AND_1E6CM2VS_UNITS == YES? mobility/100 : mobility);

   one_by_tau_total += one_by_tau_background_doping;


	  }




 /***************************************************************
  INTERFACE ROUGHNESS : All units in SI unless there is a suffix
  stating some other unit
  ***************************************************************/


  if(INCLUDE_SCATTERING_FROM_INTERFACE_ROUGHNESS == YES)

   {

    switch(HETEROINTERFACE_OR_QUANTUM_WELL)

       {
         case HETEROINTERFACE:

         one_by_tau_interface_roughness = one_by_tau_for_interface_roughness(HETEROINTERFACE_OR_QUANTUM_WELL,
                                                                             quantum_well_width_in_nanometer,
                                                                             roughness_amplitude_in_nanometer,
		    															     roughness_correl_length_in_nanometer,
                                                                             N2d,N_depletion_charge_density,
																	         m0*0.067, epsilonr,
																	         FANG_HOWARD_WAVEFN_WITH_POINT_CHARGE,
																	         FANG_HOWARD_WAVEFN_WITH_ITSELF,
										                                     DO_NOT_WRITE_VTHETA_TO_FILE,
																	         "vth_for_int_rough.dat");

         break;


		case QUANTUM_WELL:

        one_by_tau_interface_roughness = one_by_tau_for_interface_roughness(HETEROINTERFACE_OR_QUANTUM_WELL,
                                                                            quantum_well_width_in_nanometer,
                                                                            roughness_amplitude_in_nanometer,
		    															    roughness_correl_length_in_nanometer,
                                                                            N2d,N_depletion_charge_density,
																	        m0*0.067, epsilonr,
																	        COSINE_WAVEFN_WITH_POINT_CHARGE,
																	        COSINE_WAVEFN_WITH_ITSELF,
										                                    DO_NOT_WRITE_VTHETA_TO_FILE,
																	        "vth_for_int_rough.dat");


          break;


       default : /* generate error or exit ..... */

		  break;

      }

   mobility = electron_charge/(m0*0.067*one_by_tau_interface_roughness);
   fprintf(fpout,"%e\t",OUTPUT_IN_SCALED_1E10CM2_AND_1E6CM2VS_UNITS == YES? mobility/100 : mobility);


   one_by_tau_total += one_by_tau_interface_roughness;


	}


   mobility = electron_charge/(m0*0.067*one_by_tau_total);
   fprintf(fpout,"%e\n",OUTPUT_IN_SCALED_1E10CM2_AND_1E6CM2VS_UNITS == YES? mobility/100 : mobility);

  }

  fclose(fpout);
  return 0;
}





 double one_by_tau_for_a_slab_of_dopant(int heterointerface_or_quantum_well,
                                        double quantum_well_width_in_nanometer, /* ignored if heterointerface */
										double doping_start_in_nanometer,
										double doping_end_in_nanometer,
                                        double doping_density_in_per_m3,
										double N2d,
                                        double N_depletion_charge_density,
										double m_eff,
										double epsilonr,
										int option1,
										int option2,
										int write_vtheta_to_file,
										char *outfilename
										)

  {


  double  density_of_states_at_Fermi_level, bare_coulomb_prefactor,
          qtf, qtf_in_nanometer, k_Fermi,k_Fermi_in_nanometer,qtf_by_2k_Fermi;

  double tmp,q,ss2,theta,dtheta,V_theta,fq1,fq2, sum, slab_width_in_nanometer,one_by_tau;

  double  q_by_kF, w_by_wF, hbar_wF_by_kBT,  mu_by_hbar_wF;  /* for calling the 2DEG_chi routine */

  double fang_howard_b_in_nanometer,param1,param2,param3;    /* for calling the form factor routine */

  double **wavefn;

  int errcode;

  FILE *fpout;

  if(write_vtheta_to_file){fpout=fopen(outfilename,"w");  fprintf(fpout,"theta \t fq1 \t fq2 \t Vtheta\n");}

  density_of_states_at_Fermi_level = m_eff/(PI*hbar*hbar);                         /* SI */
  bare_coulomb_prefactor = electron_charge*electron_charge/(2*epsilon0*epsilonr);  /* SI */
  qtf= density_of_states_at_Fermi_level*bare_coulomb_prefactor;                    /* SI */
  qtf_in_nanometer= 1e-9*density_of_states_at_Fermi_level*bare_coulomb_prefactor;
  k_Fermi= sqrt(2*PI*N2d);                                                         /* SI */
  k_Fermi_in_nanometer = k_Fermi*1e-9;
  qtf_by_2k_Fermi=0.5*qtf/k_Fermi;

  fang_howard_b_in_nanometer = 1.0e-9*pow(((33*m_eff*electron_charge*electron_charge*N2d)/(8*hbar*hbar*epsilon0*epsilonr))+((12*m_eff*electron_charge*electron_charge*N_depletion_charge_density)/(epsilon0*epsilonr*hbar*hbar)),1.0/3.0);

   if(doping_start_in_nanometer > doping_end_in_nanometer)
	  {tmp=doping_start_in_nanometer;
       doping_start_in_nanometer = doping_end_in_nanometer;
	   doping_end_in_nanometer=tmp;
	  }


   slab_width_in_nanometer=(doping_end_in_nanometer-doping_start_in_nanometer)/1000.0;
   if(slab_width_in_nanometer < 0.5)slab_width_in_nanometer=0.5; /* slab_width < (one lattice layer) is  unphysical */


  one_by_tau=0.0;
  for(param3=doping_start_in_nanometer+0.5*slab_width_in_nanometer; param3<=doping_end_in_nanometer; param3 +=slab_width_in_nanometer)
     {

	   dtheta=PI/1000; sum=0.0;
       for(theta=1e-4; theta < PI; theta+= dtheta)
	      {
		   ss2=sin(theta/2);

	       q_by_kF=2*ss2;
           w_by_wF=0.0;             /* we need static case only here */
		   hbar_wF_by_kBT = 1000.0; /* would not be used for static zero T calculation */
		   mu_by_hbar_wF = 1.0;

           q=2*k_Fermi_in_nanometer*ss2;



           /* param3 is the co-ordinate (0,0,z) of the charged layer in nanometer */

		   switch(heterointerface_or_quantum_well)
		          {
		            case HETEROINTERFACE:

						 param1=fang_howard_b_in_nanometer;   /* param1 is the Fang Howard b paramater, has dimension of inverse length   */
	                     param2=1.00;                         /* param2 is the ratio k_insulator/k_semiconductor */

						 fq1 = FormFactor2DEG(option1, q, param1,param2, param3, wavefn, &errcode);

						 /* param3 is irrelevant for fq2, no need to reset it */
						 fq2 = FormFactor2DEG(option2, q, param1,param2, param3, wavefn, &errcode);
						 break;

					case QUANTUM_WELL:

						 param1=quantum_well_width_in_nanometer;
	                     param2=1.00;    /* not used */

						 fq1 = FormFactor2DEG(option1, q, param1,param2, param3-0.5*quantum_well_width_in_nanometer, wavefn, &errcode);

						 /* param3 is irrelevant for fq2, no need to reset it */
						 fq2=FormFactor2DEG(option2, q, param1,param2, param3, wavefn, &errcode);


						 break;

					default : /* generate error or exit...shouldn't come here */
			             break;
					/* this will be multiplied with the areal_density (slab_width_in_nanometer*1e-9*doping_density_in_per_m3)
		               at the end
					 */
                   }




           V_theta = fq1/(ss2 - qtf_by_2k_Fermi*fq2*RE_chi_2DEG(q_by_kF, w_by_wF, hbar_wF_by_kBT, mu_by_hbar_wF,STATIC_ZERO_T, &errcode));

		   /* The value returned by RE_chi_2DEG is scaled by the 2D density of states - but this is absorbed in the
		      definition of qtf - the Thomas Fermi wavevector. The e^2/2epsilon and a factor of k_Fermi has been pulled out
		   */

         if(write_vtheta_to_file){fprintf(fpout,"%e\t %e \t %e \t %e \n",theta,fq1,fq2,V_theta);}
		    sum+= V_theta*V_theta*(1-cos(theta));

           }


           one_by_tau += (slab_width_in_nanometer*1e-9*doping_density_in_per_m3)*(m_eff/(PI*hbar*hbar*hbar))*pow(0.5*bare_coulomb_prefactor/k_Fermi,2)*sum*dtheta;

	  }

  if(write_vtheta_to_file){fclose(fpout);}

  return one_by_tau;

  }


 double one_by_tau_for_a_delta_doped_layer(int heterointerface_or_quantum_well,
                                           double quantum_well_width_in_nanometer,  /* ignored if QW */
										   double doping_location_in_nanometer,
										   double doping_density_in_per_m2,
										   double N2d,
                                           double N_depletion_charge_density,
                                           double m_eff,
										   double epsilonr,
										   int option1,
										   int option2,
										   int write_vtheta_to_file,
										   char *outfilename
										   )

  {


  double  density_of_states_at_Fermi_level, bare_coulomb_prefactor, qtf, qtf_in_nanometer, k_Fermi,k_Fermi_in_nanometer,
          qtf_by_2k_Fermi;
  double tmp,q,ss2,theta,dtheta,V_theta,fq1,fq2, sum,one_by_tau;
  double  q_by_kF, w_by_wF, hbar_wF_by_kBT,  mu_by_hbar_wF;  /* for calling the 2DEG_chi routine */
  double fang_howard_b_in_nanometer,param1,param2,param3;    /* for calling the form factor routine */
  double **wavefn;

  int errcode;

  FILE *fpout;

  if(write_vtheta_to_file){fpout=fopen(outfilename,"w");  fprintf(fpout,"theta \t fq1 \t fq2 \t Vtheta\n");}


  density_of_states_at_Fermi_level = m_eff/(PI*hbar*hbar);                         /* SI */
  bare_coulomb_prefactor = electron_charge*electron_charge/(2*epsilon0*epsilonr);  /* SI */
  qtf= density_of_states_at_Fermi_level*bare_coulomb_prefactor;                    /* SI */
  qtf_in_nanometer= 1e-9*density_of_states_at_Fermi_level*bare_coulomb_prefactor;
  k_Fermi= sqrt(2*PI*N2d);                                                         /* SI */
  k_Fermi_in_nanometer = k_Fermi*1e-9;
  qtf_by_2k_Fermi=0.5*qtf/k_Fermi;


  fang_howard_b_in_nanometer =  1.0e-9*pow(((33*m_eff*electron_charge*electron_charge*N2d)/(8*hbar*hbar*epsilon0*epsilonr))+((12*m_eff*electron_charge*electron_charge*N_depletion_charge_density)/(epsilon0*epsilonr*hbar*hbar)),1.0/3.0);



       one_by_tau=0.0;
	   dtheta=PI/1000; sum=0.0;
       for(theta=1e-4; theta < PI; theta+= dtheta)
	      {
		   ss2=sin(theta/2);

	       q_by_kF=2*ss2;
           w_by_wF=0.0;             /* we need static case only here */
		   hbar_wF_by_kBT = 1000.0; /* would not be used for static zero T calculation */
		   mu_by_hbar_wF = 1.0;

           q=2*k_Fermi_in_nanometer*ss2;


           param3=doping_location_in_nanometer; /* co-ordinate (0,0,z) of the charged layer in nanometer */


		   switch(heterointerface_or_quantum_well)
		          {
		            case HETEROINTERFACE:
						 param1=fang_howard_b_in_nanometer;   /* param1 is the Fang Howard b paramater, has dimension of inverse length   */
	                     param2=1.00;                         /* param2 is the ratio k_insulator/k_semiconductor */

						 fq1 = FormFactor2DEG(option1, q, param1,param2, param3, wavefn, &errcode);
                         /* param3 is irrelevant for fq2, no need to reset it */
						 fq2=FormFactor2DEG(option2, q, param1,param2, param3, wavefn, &errcode);
						 break;

					case QUANTUM_WELL:

						 param1=quantum_well_width_in_nanometer;
	                     param2=1.00;    /* not used */

						 fq1 = FormFactor2DEG(option1, q, param1,param2, param3, wavefn, &errcode);
						 /* param3 is irrelevant for fq2, no need to reset it */
						 fq2=FormFactor2DEG(option2, q, param1,param2, param3, wavefn, &errcode);


						 break;

					default : /* generate error or exit...shouldn't come here */
			             break;
					/* this will be multiplied with the areal_density (slab_width_in_nanometer*1e-9*doping_density_in_per_m3)
		               at the end
					 */
                   }


           V_theta = fq1/(ss2 - qtf_by_2k_Fermi*fq2*RE_chi_2DEG(q_by_kF, w_by_wF, hbar_wF_by_kBT, mu_by_hbar_wF,STATIC_ZERO_T, &errcode));

		   /* The value returned by RE_chi_2DEG is scaled by the 2D density of states - but this is absorbed in the
		      definition of qtf - the Thomas Fermi wavevector. The e^2/2epsilon and a factor of k_Fermi has been pulled out
		   */
      if(write_vtheta_to_file){fprintf(fpout,"%e\t %e \t %e \t %e \n",theta,fq1,fq2,V_theta);}

		    sum+= V_theta*V_theta*(1-cos(theta));

           }



           one_by_tau = doping_density_in_per_m2*(m_eff/(PI*hbar*hbar*hbar))*pow(0.5*bare_coulomb_prefactor/k_Fermi,2)*sum*dtheta;

       if(write_vtheta_to_file){fclose(fpout);}


  return one_by_tau;

  }



 double one_by_tau_for_interface_roughness(int heterointerface_or_quantum_well,
                                           double quantum_well_width_in_nanometer, /* ignored if heterointerface */
                                           double roughness_amplitude_in_nanometer,
										   double roughness_correl_length_in_nanometer,
                                           double N2d,
										   double N_depletion_charge_density,
										   double m_eff,
										   double epsilonr,
										   int option1,
										   int option2,
										   int write_vtheta_to_file,
										   char *outfilename
										   )
 {


  /* Following  T. Ando, Journal of the Physical Society of Japan, vol 43, no 5 p 1616 (1977)
     But converted to SI for Coulomb interaction etc.

    For QW: H. Sakaki et al APL 51, 1934 (1987) But this does not treat the shape of the wavefn in
    the QW explicitly........a formulation with F(q) is needed.

  */

  double  density_of_states_at_Fermi_level, bare_coulomb_prefactor, qtf, qtf_in_nanometer, k_Fermi,k_Fermi_in_nanometer,
          qtf_by_2k_Fermi;
  double tmp,q,ss2,theta,dtheta,V_theta,fq1,fq2, sum,one_by_tau;
  double  q_by_kF, w_by_wF, hbar_wF_by_kBT,  mu_by_hbar_wF;  /* for calling the 2DEG_chi routine */
  double fang_howard_b_in_nanometer,param1,param2,param3,gammaq,k_Fermi_mult_roughness_correl_length_sq;    /* for calling the form factor routine */
  double **wavefn;

  int errcode;

  FILE *fpout;

  if(write_vtheta_to_file){fpout=fopen(outfilename,"w");  fprintf(fpout,"theta \t fq1 \t fq2 \t Vtheta\n");}

  density_of_states_at_Fermi_level = m_eff/(PI*hbar*hbar);                         /* SI */
  bare_coulomb_prefactor = electron_charge*electron_charge/(2*epsilon0*epsilonr);  /* SI */
  qtf= density_of_states_at_Fermi_level*bare_coulomb_prefactor;                    /* SI */
  qtf_in_nanometer= 1e-9*qtf;
  k_Fermi= sqrt(2*PI*N2d);                                                         /* SI */
  k_Fermi_in_nanometer = k_Fermi*1e-9;
  qtf_by_2k_Fermi=0.5*qtf/k_Fermi;
  k_Fermi_mult_roughness_correl_length_sq=k_Fermi_in_nanometer*roughness_correl_length_in_nanometer;
  k_Fermi_mult_roughness_correl_length_sq *= k_Fermi_mult_roughness_correl_length_sq;


  fang_howard_b_in_nanometer = 1.0e-9*pow(((33*m_eff*electron_charge*electron_charge*N2d)/(8*hbar*hbar*epsilon0*epsilonr))+((12*m_eff*electron_charge*electron_charge*N_depletion_charge_density)/(epsilon0*epsilonr*hbar*hbar)),1.0/3.0);
  gammaq=bare_coulomb_prefactor*(N2d+2*N_depletion_charge_density);    /* SI */
                                       /* since this is independent of q it is pulled out of integration  */
                                       /* If dielectric ratio is not equal to 1, then this is not correct */

       one_by_tau=0.0;
	   dtheta=PI/1000; sum=0.0;
       for(theta=1e-4; theta < PI; theta+= dtheta)
	      {
		   ss2=sin(theta/2);

	       q_by_kF=ss2+ss2;
           w_by_wF=0.0;             /* we need static case only here */
		   hbar_wF_by_kBT = 1000.0; /* would not be used for static zero T calculation */
		   mu_by_hbar_wF = 1.0;

           q=2*k_Fermi_in_nanometer*ss2;
	       param1=fang_howard_b_in_nanometer;   /* param1 is the Fang Howard b paramater, has dimension of inverse length   */
	       param2=1.00;                         /* param2 is the ratio k_insulator/k_semiconductor */
           param3=1.00;                         /* co-ordinate (0,0,z) of the charged layer in nanometer, irrelevant here */



		    /* fq1 = exp(-0.125*k_Fermi_mult_roughness_correl_length_sq*ss2*ss2);*/

		       fq1 = exp(-0.5*k_Fermi_mult_roughness_correl_length_sq*ss2*ss2);

		    /* param3 is irrelevant for fq2, no need to reset it */
		       fq2=FormFactor2DEG(option2, q, param1,param2, param3, wavefn, &errcode);






           V_theta = fq1*ss2/(ss2 - qtf_by_2k_Fermi*fq2*RE_chi_2DEG(q_by_kF, w_by_wF, hbar_wF_by_kBT, mu_by_hbar_wF,STATIC_ZERO_T, &errcode));

		   /* The value returned by RE_chi_2DEG is scaled by the 2D density of states - but this is absorbed in the
		      definition of qtf - the Thomas Fermi wavevector. The e^2/2epsilon and a factor of k_Fermi has been pulled out
		   */
            if(write_vtheta_to_file){fprintf(fpout,"%e\t %e \t %e \t %e \n",theta,fq1,fq2,V_theta);}

		    sum+= V_theta*V_theta*(1-cos(theta));

           }




		switch(heterointerface_or_quantum_well)
		     {  case HETEROINTERFACE:

                     one_by_tau = (m_eff/(hbar*hbar*hbar))*(1e-36)*pow(gammaq*roughness_amplitude_in_nanometer*roughness_correl_length_in_nanometer,2)*sum*dtheta;

					 break;

			    case QUANTUM_WELL:

					 one_by_tau =  2.0*(hbar/m_eff)*(1e18)*pow(quantum_well_width_in_nanometer,-6)*pow(PI*PI*roughness_amplitude_in_nanometer*roughness_correl_length_in_nanometer,2)*sum*dtheta;

					 break;

			    default :   /* generate error or exit....shouldn't come here. */
				     break;


			}



            if(write_vtheta_to_file){fclose(fpout);}


  return one_by_tau;


  }



