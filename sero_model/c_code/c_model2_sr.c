 /* Set up likelihood for vector of infection history */

#include <math.h>
#include <stdio.h>

/* Inputs needed:
	- test strain data vector
	- infection history vector
	- d.ij vector
 */

void c_model2_sr(int *nin, int *itot, int *nsin, double *x, double *x1, double *titre, 
                  double *titrepred, double *dd, int *ntheta, 
                  double *theta, int *inputtestyr)
{
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Calculate lambda */
	
	int n = nin[0];
	int total_inf = itot[0];
	int nsamp = nsin[0];
	double T_1 = theta[1];
	double T_2 = theta[2];
	double wane = theta[3];
	double mu = theta[0];
	
	// This to be made an argument of the function -- gives test year
	int t_sample = inputtestyr[0]; 
  
  double yrTitre[n*nsamp];
  int maskedInfectionHistory[n];
  double distanceFromTest[n];
  double cumInfectionHistory[n];
  
	/* Add for loop over k*/
	
	int k;
	int i;
	int j;
	int m;
	
	double xx2; 
	
	for (k=0; k<nsamp; k++){
		
	  for (j=0; j<n; j++) {
	
		xx2=0;

	  // Make a masked infection history
	  for (m=0;m<n;m++) {
	    if (m <= j) {
	      maskedInfectionHistory[m]=x[m];
	    } else {
	      maskedInfectionHistory[m]=0;
	    }
	  }
	  
	  // Make an index for waning
	  for (m=0;m<n;m++) {
		  distanceFromTest[m]=exp(-wane * (t_sample-m));
	  }

	  // Make a cumulative infection history
	  cumInfectionHistory[0] = maskedInfectionHistory[0];
	  for (m=1;m<n;m++) {
	    cumInfectionHistory[m] = cumInfectionHistory[m-1] + 
	      maskedInfectionHistory[m];
	  }
	  	    
		/* Calculate expected titre	- note k indexed from 0 */

		for (i=0; i<n; i++){
			x1[i] = mu *
			  dd[k*n+i] * 
			  maskedInfectionHistory[i] *
			  exp(-1.0 * T_2 * ( cumInfectionHistory[i]  - 1.0)) *
			  pow(1+T_1 , (total_inf - cumInfectionHistory[i])) *
			  distanceFromTest[i];
		}
	
		for (i=0; i<n; i++){
			xx2 =  xx2 + x1[i];
		}
	
	  yrTitre[k+j*nsamp]=xx2;
	
	  } // end test year loop (j)
	
	} // end sample loop (k)
	
	for (k=0;k<nsamp;k++) {
	  titrepred[k]=yrTitre[k+(t_sample-1)*nsamp]; //Need (t-1) as index from 0
	}
	
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Compare to observed titre*/

}