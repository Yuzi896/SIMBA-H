#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
void MCMCfunction(int *Method, int *OutcomeType, int *M1, int *p1, int *nn, double *Data1, double *Cov1, int * pc,
double *Youtcome, int *KK, double *truncatedvalue, double *rhosigma,
double *gamsigma, double *alphaP, double *delt, double *sigma20Mu1,
int *Nsample1, int *burnin1, int *seed1, double *probZMean1,
double *probkappajk1, double *TransMatMean1, double *MuMean1,
double *CoefYsample, double *sig2yMean, double *logpostsample,
double * Betasample, double *DIC,int *nnTest,double *DataTest1,
double *CovTest1, double * YpredTest, int *Ztest);
void logposterior(int n, int M, int *p, int pc, float ***dat, double *Y, float **Cov,
                  float *Beta, float *CoefY, float sig2y, float ****Mu, float ***Muk,
                  float ***sigma2Muk, float ****sigma2, float **sigma2j3, int *Zk,
                  int ***kappajk, float ***TransMat, float *Omeg, double *alphaP,
                  float truncval, float alpha, float rho, float gam, double sigma20Mu,
                  double *delt, int K, double *logpost1, double *loglik1);
float LGammaFF(float x);
float max(int n, float *x);
void findc(int n, int R[n], int a, int *IDX, int *nx);
int maxIdx(int n, double *x);
float Mean(int n, float *x);
void SampleCoefYSig2y(double *Y, float *CoefY, int *Zk, float *sig2y, float sig2a, int K, int n, float aa, float ba, int OutcomType);
void SampleZ(int * Zk,int M, int n, int K,int* p,int ***kappajk,float ****sigma2,float ***dat,float ****Mu,float **sigma2j3, double *Y,
float *CoefY,  float sig2y, float *Omeg, int pred, double *probSample);
void SampleKappa(int K, int p, int n, int *Zk, int **kappajk, float ***Mu, float ***sigma2, float **TransMat, float **dat, float *sigma2j3, double *probSample);
void ClusterProba(float *Omeg, int K, int n, double *alphaP,int *Zk);
void SampleMuSigma(int K, int p, int n, int *Zk, int **kappajk, float ***Mu,  float **Muk, float ***sigma2, float **sigma2Muk,
float **TransMat, float **dat, float *sigma2j3, double rho, double gam, double alpha, float sigma20Mu, float truncval);
void SampleTransitionMatrice(int K, int p, int Method1, float **TransMat, int **kappajk, double *delt);
float *cholesky(float *A, int n);
void LinearSystem(float *L,int n, float *x, float *b);
void BackSubstitution (float * L, int n,float *x,float *b);
void ForwardSubstitution (float * L,int n,float *x,float *b);
void SampleFixedEffectsY(double *Y, int *Zk, int K, int n, float *CoefY, float *Beta, float **Cov, int pr);
double *ResiXy(double *Y, int n, float *Beta, float **Cov, int pc);
float *  invertLT( float * L , int n);
void SampleLatentYBinary(int pc, int K, int n, double *ylatent, double *y, float **Cov, float *CoefY, int *Zk, float *Beta);
void SampleLatentYRightCensored(int pc, int K, int n, double *ylatent, double *Time, double *Cens,
float **Cov, float *CoefY, int *Zk, float *Beta, float sig2y);
void YPredict(float *Ypred, int *Zk, int pc, float **Cov, int n, int K, float *Beta, float *CoefY);


void MCMCfunctionWithProbSamples(int *Method, int *OutcomeType, int *M1, int *p1, int *nn, double *Data1,
                  double *Cov1, int *pc, double *Youtcome, int *KK, double *truncatedvalue,
                  double *rhosigma, double *gamsigma, double *alphaP, double *delt,
                  double *sigma20Mu1, int *Nsample1, int *burnin1, int *seed1, double *probZMean1,
                  double *probkappajk1, double *TransMatMean1, double *MuMean1,
                  double *CoefYsample, double *sig2yMean, double *logpostsample,
                  double *Betasample, double *DIC, int *nnTest, double *DataTest1,
                  double *CovTest1, double *YpredTest, int *Ztest, double *MuSample, double *probZSample, double *probKappaSample);

void MCMCfunctionWithProbAndZSamples(int *Method, int *OutcomeType, int *M1, int *p1, int *nn, double *Data1,
                  double *Cov1, int *pc, double *Youtcome, int *KK, double *truncatedvalue,
                  double *rhosigma, double *gamsigma, double *alphaP, double *delt,
                  double *sigma20Mu1, int *Nsample1, int *burnin1, int *seed1, double *probZMean1,
                  double *probkappajk1, double *TransMatMean1, double *MuMean1,
                  double *CoefYsample, double *sig2yMean, double *logpostsample,
                  double *Betasample, double *DIC, int *nnTest, double *DataTest1,
                  double *CovTest1, double *YpredTest, int *Ztest, double *MuSample, double *probZSample, double *probKappaSample, int *Zsample);
