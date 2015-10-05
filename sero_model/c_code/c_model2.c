/* Set up likelihood for vector of infection history */

#include <math.h>
#include <stdio.h>

/* Inputs needed:
	- test strain data vector
	- infection history vector
	- d.ij vector
 */

void c_model2(int *nin, int *nsin, double *x, double *x1, double *titre, double *titrepred, double *dd, double *mu)
{
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Calculate lambda */
	
	int n = nin[0];
	int nsamp = nsin[0];
	

	/* Add for loop over k*/
	
	int k;
	int i;
	
	double xx2; 
	
	for (k=0; k<nsamp; k++){
		
		xx2=0;

		/* Calculate expected titre	*/

		for (i=0; i<n; i++){
			x1[i] =  mu[0] * dd[k*n+i] * x[i];
		}
	
		for (i=0; i<n; i++){
			xx2 =  xx2 + x1[i];
		}
	
		titrepred[k]=xx2;
	
	}
	
	/* End loop over k*/
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Compare to observed titre*/

	
	
}