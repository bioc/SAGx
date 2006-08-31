/*samrocNboot.c is replacing samrocN and is linked with samrocNboot.r */

#include <R.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

void samrocNboot(double *indata,int *NR,int *NC,double *xhat,int *NRP ,int *NCP,double *xmat,double *contrast,int *smooth,int *B, double *npar,double *scalek,double *ss,double *diffs,double *sses){  /*NR =ncol, *NC=nrow*/
   int i,j,k,I;
   int *result = (int *)malloc(*NCP * sizeof(int));
   double *P = (double *)malloc(*NRP*(*NCP) * sizeof(double));
   double *beta = (double *)malloc(*NRP*(*NR) * sizeof(double));
   double *yhat = (double *)malloc(*NR*(*NC) * sizeof(double));
   double *resid = (double *)malloc(*NR*(*NC) * sizeof(double));
   double *var = (double *)malloc(*NR * sizeof(double));

/*Here is where the bootstrap loop should begin*/
  for(I=0;I<*B;I=I+1){							

/*Random number generator no replacement*/							
 /* srand((long)time(NULL)); *//*thisline need to be included for randomness between tries*/							
  	for(i=0;i<*NCP;i=i+1){						
  		resample:					
 		result[i] = (rand () %(*NCP));					
  		for(j=0;j<=i-1;j=j+1){ 					
			if(result[i]==result[j]){				
				goto resample;			
							}
					    }		
                         }							


/*Rprintf(""\n\nResampled P* \n\n"");*/							
/*Make resampled P*=(xhat) */	
 for(i=0;i<(*NRP); i=i+1){   	
    	for(j=0;j<(*NCP); j=j+1){	
    		P[i+j*(*NRP)]=xhat[i+result[j]*(*NRP)]; 	
     		/*Rprintf("%8.2f",P[i+j*(*NRP)]);*/
					}
     /*Rprintf("\n");*/ 
				}	

/*beta<-crossprod(t(P*),t(indata))=P*%*%t(indata)*/	
/*Rprintf(""\n\nbeta \n\n"");*/	

for(i=0;i<*NRP; i=i+1){	
   for(j=0;j<*NR;j=j+1){	
      beta[i+j*(*NRP)]=0;	
      for(k=0;k<*NC;k=k+1){	
	beta[i+j*(*NRP)]=beta[i+j*(*NRP)]+P[i+k*(*NRP)]*indata[j+k*(*NR)];   
                        }	

   /*Rprintf("%8.2f",beta[i+j*(*NRP)]);*/}
     /*Rprintf("\n");*/ }

   	
/*crossprod (contrast,beta)*/	
/*Rprintf(""\n\ndiffs \n\n"");*/	

    for(j=0;j<*NR;j=j+1){	
      diffs[j+I*(*NR)]=0;	
      for(k=0;k<*NRP;k=k+1){	
	diffs[j+I*(*NR)]=diffs[j+I*(*NR)]+contrast[k]*beta[k+j*(*NRP)];
}   	
                         /* Rprintf("%8.2f",diffs[j+I*(*NR)]);*/}	
    	



/*yhat<-t(beta)%*%t(xmat)*/	
/*Rprintf(""\n\nyhat \n\n"");*/
for(i=0;i<*NR; i=i+1){	
   for(j=0;j<*NC;j=j+1){	
      yhat[i+j*(*NR)]=0;	
      for(k=0;k<*NRP;k=k+1){	
	yhat[i+j*(*NR)]=yhat[i+j*(*NR)]+beta[k+i*(*NRP)]*xmat[j+k*(*NCP)];   
                        }	
      resid[i+j*(*NR)]=indata[i+j*(*NR)]-yhat[i+j*(*NR)];	
      	
     /*Rprintf("%8.2f",yhat[i+j*(*NR)]);*/}
    /*Rprintf("\n");*/ }	
       	
      	

/*var*/
 /*Rprintf(""\n\nvar \n\n"");*/
 for(i=0;i<*NR;i=i+1){
     var[i]=0;
       for(j=0;j<*NC;j=j+1){ var[i]=var[i]+(resid[i+j*(*NR)])*(resid[i+j*(*NR)]);} 
       var[i]=var[i]/((*NC)-(*npar));
/*Rprintf(""%10.2f"",var[i]);*/}
         


/*Calculate sses*/
/*Rprintf(""\n\n sses \n\n"");*/
if(*smooth==1){  /*1<-> TRUE, 0 <-> FALSE */
  for(i=0;i<*NR;i=i+1){ sses[i+I*(*NR)]=ss[i];
 /*Rprintf(""%8.2f"",sses[i+I*(*NR)]);*/}}

 else{  for(i=0;i<*NR;i=i+1){ sses[i+I*(*NR)]=sqrt(var[i]/(*scalek));	
 /*Rprintf(""%8.2f"",sses[i+I*(*NR)]);*/}}	
 	
  	
}	
free(result); free(P); free(beta); free(yhat); free(resid); free(var);

/*End of bootstrap loop*/	
	}	
