/*MiniGSEA replaces the bootstrap loop in gsea, and is linked with gseawcan.r */
#include <R.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

void MiniGSEA(double *indata,int*B,int *N1,int *N2, int *NR ,int *NC,double *s0,double *d1){  /*NR =ncol, *NC=nrow*/
   int i,j,k,I;
   int *result = (int *)malloc(*NC * sizeof(int));
   double *data1= (double *)malloc(*NR *(*N1) * sizeof(double));
   double *data2= (double *)malloc(*NR *(*N2) * sizeof(double));
   double *m1= (double *)malloc(*NR * sizeof(double));
   double *m2=  (double *)malloc(*NR * sizeof(double));       /*Vector of rowmeans of data1 and data2 resp.*/
   double *var1= (double *)malloc(*NR * sizeof(double));
   double *var2= (double *)malloc(*NR * sizeof(double));
   double *SE= (double *)malloc(*NR * sizeof(double));
   double N1d,N2d;

  

/*Here is where the bootstrap loop should begin*/
  for(I=0;I<*B;I=I+1){

/*Random number generator no replacement*/
 /* srand((long)time(NULL));*/ /*thisline need to be included for randomness between tries*/
  for(i=0;i<*NC;i=i+1){
  resample:
  /*result[i]=i+1;*/
  result[i] = (rand () %(*NC));
  for(j=0;j<=i-1;j=j+1){ if(result[i]==result[j]){goto resample;}}
  /*Rprintf("result=%d\n",result[i]);*/
}

/*Make data1 and data2*/
 /*Rprintf("\n\ndata1\n\n");*/
for(i=0;i<(*NR); i=i+1){   
    for(j=0;j<(*N1); j=j+1){
    data1[i+j*(*NR)]=indata[i+result[j]*(*NR)]; 
     /* Rprintf("%6.2f",data1[i][j]);*/}
      /*Rprintf("\n"); */}

/*Rprintf("\n\ndata2\n\n"); */

for(i=0;i<(*NR); i=i+1){   
    for(j=0;j<(*N2); j=j+1){
    data2[i+j*(*NR)]=indata[i+result[j]*(*NR)]; 
      /*Rprintf("%6.2f",data2[i][j]);*/}
    /*  Rprintf("h\n")*/;
       }


/*Calculate rowmeans*/
for(i=0;i<*NR;i=i+1){
     m1[i]=0; m2[i]=0;
       for(j=0;j<*N1;j=j+1){ m1[i]=m1[i]+data1[i+j*(*NR)];}
     m1[i]=m1[i]/(*N1);
       for(j=0;j<*N2;j=j+1){ m2[i]=m2[i]+data2[i+j*(*NR)];}
     m2[i]=m2[i]/(*N2);
   }
      /* for(i=0;i<*NR;i=i+1){Rprintf("%1.2f\n",m2[i]);} Rprintf("\n");*/
       /*Variances*/
     for(i=0;i<*NR;i=i+1){
     var1[i]=0; var2[i]=0;
       for(j=0;j<*N1;j=j+1){ var1[i]=var1[i]+(data1[i+j*(*NR)]-m1[i])*(data1[i+j*(*NR)]-m1[i]);} 
       for(j=0;j<*N2;j=j+1){ var2[i]=var2[i]+(data2[i+j*(*NR)]-m2[i])*(data2[i+j*(*NR)]-m2[i]);}
     }
     /*  for(i=0;i<*NR;i=i+1){Rprintf("%1.2f\n",var2[i]);} Rprintf("\n");*/
N1d=*N1; N2d=*N2;      

     for(i=0;i<*NR;i=i+1){
       SE[i]=sqrt((var1[i]+var2[i])/(N1d+N2d-2))*sqrt(1/(N1d)+1/(N2d));
       d1[i+I*(*NR)]=(m1[i]-m2[i])/(SE[i]+(*s0));
       k=I*(*NR)+i;
       /*Rprintf("%6.2f %6.2f d1[%d]\n",SE[i],d1[I*(*NR)+i], k);*/}

    }
/*End of bootstrap loop*/
free(result); free(data1); free(data2); free(m1); free(m2); free(var1); free(var2); free(SE);


}
