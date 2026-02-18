#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <stdio.h>
#include "utils.h"
#include "distribution.h"

float *  invertLT( float * L , int n)
{
   float *iLT = (float*)calloc(n * n, sizeof(float));
   

   for ( int i = 0; i < n; i++ ) 
   {
      if ( L[i*n+i] == 0.0 )
      {
       //  cout << "*** Singular matrix ***\n";
       //  return M;
      }
      iLT[i*n+i] = 1.0 / L[i*n+i];
      for ( int j = 0; j < i; j++ )
      {
         for ( int k = j; k < i; k++ ) iLT[i*n+j] += L[i*n+k] * iLT[k*n+j];
         iLT[i*n+j] = -iLT[i*n+j] / L[i*n+i];
      }
   }
   return iLT;
}


float *cholesky(float *A, int n) {
  // A=LL^T; cholesky decomposition
    float *L = (float*)calloc(n * n, sizeof(float));
    if (L == NULL)
        exit(EXIT_FAILURE);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < (i+1); j++) {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += L[i * n + k] * L[j * n + k];
            L[i * n + j] = (i == j) ?
                           sqrt(A[i * n + i] - s) :
                           (1.0 / L[j * n + j] * (A[i * n + j] - s));
        }

    return L;
}

void ForwardSubstitution (float * L,int n,float *x,float *b){
for (int i=0;i<n;i++){
float xx=0;
int j;
for (j=0;j<i;j++)   xx+=L[i * n + j]*x[j];
x[i]=(b[i]-xx)/L[i * n + i];}
}

  void BackSubstitution (float * L, int n,float *x,float *b){
   //System L^Tx=b; L=lower triangular matrix
   int i,j;
   for (i=0;i<n;i++){
   float xx=0;
   for (j=0;j<i;j++)
       xx+=L[(n-1-i)  + (n-1-j)*n]*x[n-1-j];
   x[n-1-i]=(b[n-1-i]-xx)/L[(n-1-i) * n + (n-1-i)];
   }
   }

void LinearSystem(float *L,int n, float *x, float *b){
    //Resolve Ax=b; with A=LL^T; cholesky decomposition
 ForwardSubstitution (L, n,x,b);
 BackSubstitution (L,n,x,x);

}

float max(int n, float * x){
float xmax=x[0];
int i;
for (i=0;i<n;i++){
if (x[i]>xmax)
xmax=x[i];
}
return xmax;
}

int maxIdx(int n, double * x){
    double xmax=x[0];int j=0;
int i;
for (i=0;i<n;i++){
    if (x[i]>xmax){
        xmax=x[i];j=i;
    }
}
return j;
}

// Find indexes that are different from a

void findc(int n,int R[n],int a,int * IDX, int *nx)
{
  int ii_data[n];
  int idx = 0;
  int  ii = 0;
  _Bool  exitg2 = 0;
  _Bool guard2=0;
  while ((exitg2 == 0) && (ii < n)) {
    guard2 = 0;
    if (R[ii]!= a) {
      idx++;
      ii_data[idx - 1] = ii;
      if (idx >= n) {
        exitg2 = 1;
      } else {
        guard2 = 1;
      }
    } else {
      guard2 = 1;
    }
    
    if (guard2 == 1) {
      ii++;
    }
  }
  
  int loop_ub=idx;
  for (idx = 0; idx < loop_ub; idx++) {
    IDX[idx] = ii_data[idx];
  }
  *nx=loop_ub;
  
}

float Mean(int n, float *x){
    float mean=0;
    for (int i=0;i<n;i++){
        mean+=x[i];
    }
    return mean/n;
}

