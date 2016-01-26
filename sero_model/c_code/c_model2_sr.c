 /* Set up likelihood for vector of infection history */

#include <math.h>
#include <stdio.h>

/* Inputs needed:
	- test strain data vector
	- infection history vector
	- d.ij vector
 */

void c_model2_sr(int *nin, int *nsin, double *x, double *x1, double *titre, 
                  double *titrepred, double *dd, double *mu)
{
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Calculate lambda */
	
	int n = nin[0];
	int nsamp = nsin[0];
	
	// This to be made an argument of the function
	int t_sample = n; 
  
  double yrTitre[n*nsamp]; 	

	/* Add for loop over k*/
	
	int k;
	int i;
	int j;
	
	double xx2; 
	
	for (k=0; k<nsamp; k++){
		
	  for (j=0; j<n; j++) {
	
		xx2=0;

		/* Calculate expected titre	- note k indexed from 0 */

		for (i=0; i<n; i++){
			x1[i] =  mu[0] * dd[k*n+i] * x[i];
		}
	
		for (i=0; i<n; i++){
			xx2 =  xx2 + x1[i];
		}
	
	  yrTitre[k+j*nsamp]=xx2;
	
	  }
	
	}
	
	for (k=0;k<nsamp;k++) {
	  titrepred[k]=yrTitre[k+(t_sample-1)*nsamp];
	}
	
	/* End loop over k*/
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Compare to observed titre*/

}