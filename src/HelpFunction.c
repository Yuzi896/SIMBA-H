#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "distribution.h"
#include "GammaFunction.hpp"
#include "utils.h"

void SampleLatentYBinary(int pc, int K, int n, double *ylatent, double *y, float **Cov, float *CoefY, int *Zk, float *Beta)
{
  for (int i = 0; i < n; i++)
  {
    float Mean = 0;
    for (int k = 0; k < K; k++)
    {
      if (Zk[i] == (k + 1))
      {
        Mean += CoefY[k];
      }
    }
    for (int l = 0; l < pc; l++)
    {
      Mean += Cov[i][l] * Beta[l];
    }
    if (y[i] == 1)
    {
      ylatent[i] = r_lefttruncnorm(0, Mean, 1.0);
    }
    else
    {
      ylatent[i] = r_righttruncnorm(0, Mean, 1.0);
      // printf(" %.2lf ", ylatent[i]);
    }
  }
}

void SampleLatentYRightCensored(int pc, int K, int n, double *ylatent, double *Time, double *Cens,
                                float **Cov, float *CoefY, int *Zk, float *Beta, float sig2y)
{
  for (int i = 0; i < n; i++)
  {
    float Mean = 0;
    for (int k = 0; k < K; k++)
    {
      if (Zk[i] == (k + 1))
      {
        Mean += CoefY[k];
      }
    }
    for (int l = 0; l < pc; l++)
    {
      Mean += Cov[i][l] * Beta[l];
    }
    if (Cens[i] == 1)
      ylatent[i] = r_lefttruncnorm(Time[i], Mean, sqrt(sig2y));
  }
}

void SampleTransitionMatrice(int K, int p, int Method1, float **TransMat, int **kappajk, double *delt)
{
  // if Method1==1, then it's HMM, otherwise, we don't account for HMM
  if (Method1 != 1)
  {
    int *Nr = (int *)calloc(3, sizeof(int));
    // Nr=rep(0,3)
    float *alphaD = (float *)calloc(3, sizeof(float));
    for (int s = 0; s < 3; s++)
    {
      for (int k = 0; k < K; k++)
        for (int j = 0; j < p; j++)
          Nr[s] = Nr[s] + (kappajk[j][k] == (s + 1));
      // Nr[s]=sum(kappajk==s)
      alphaD[s] = Nr[s] + delt[s];
    }
    // ma=rdirichlet(1,as.vector(delt+Nr))
    float *theta = calloc(3, sizeof(float));
    gsl_ran_dirichlet(3, alphaD, theta);
    // TransMat=matrix(c(rep(ma,3)),3,byrow = T)
    for (int l = 0; l < 3; l++)
      for (int l1 = 0; l1 < 3; l1++)
        TransMat[l][l1] = theta[l1];
    free(Nr);
    free(alphaD);
    free(theta);
  }
  else
  {
    for (int r = 0; r < 3; r++)
    {
      // Nr=rep(0,3)
      int *Nr = (int *)calloc(3, sizeof(int));
      // Nr=rep(0,3)
      float *alphaD = (float *)calloc(3, sizeof(float));
      for (int s = 0; s < 3; s++)
      {
        for (int j = 0; j < (p - 1); j++)
        {
          for (int k = 0; k < K; k++)
            Nr[s] = Nr[s] + ((kappajk[j][k] == (r + 1)) && (kappajk[j + 1][k] == (s + 1)));
          //            Nr[s]=Nr[s]+sum((kappajk[j,]==r) & (kappajk[j+1,]==s))
        }
        alphaD[s] = Nr[s] + delt[s];
      }
      gsl_ran_dirichlet(3, alphaD, TransMat[r]);
      // TransMat[r,]=rdirichlet(1,as.vector(delt+Nr))
      free(Nr);
      free(alphaD);
    }
  }
}

void SampleMuSigma(int K, int p, int n, int *Zk, int **kappajk, float ***Mu, float **Muk, float ***sigma2, float **sigma2Muk,
                   float **TransMat, float **dat, float *sigma2j3, double rho, double gam, double alpha, float sigma20Mu, float truncval)
{
  for (int j = 0; j < p; j++)
  {
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < K; k++)
      {
        int Nkl = 0;
        float Skl = 0;
        // Nkl=(kappajk[j,k]==l)*sum(Zk==k)
        for (int i = 0; i < n; i++)
        {
          Nkl = Nkl + (kappajk[j][k] == (l + 1)) * (Zk[i] == (k + 1));
          // Skl=(kappajk[j,k]==l)*sum((Zk==k)*dat[,j]);
          Skl = Skl + (kappajk[j][k] == (l + 1)) * (Zk[i] == (k + 1)) * dat[i][j];
        }
        if (l < 2)
        {
          float Mean = (Muk[k][l] * sigma2[j][k][l] / sigma2Muk[k][l] + Skl) / (sigma2[j][k][l] / sigma2Muk[k][l] + Nkl);
          float SD2 = 1 / (Nkl / sigma2[j][k][l] + 1 / sigma2Muk[k][l]);

          Mu[j][k][l] = Mean + GenerateNormal(sqrt(SD2));
        }

        if (l < 2)
        {
          float S2kl = 0;
          for (int i = 0; i < n; i++)
            S2kl = S2kl + (kappajk[j][k] == (l + 1)) * (Zk[i] == (k + 1)) * pow(dat[i][j] - Mu[j][k][l], 2);
          // S2kl=(kappajk[j,k]==l)*sum((Zk==k)*(dat[,j]-Mu[j,k,l])^2)

          sigma2[j][k][l] = 1 / ran_gamma_knuth(rho + Nkl / 2.0, 1 / (gam + S2kl / 2));
          // sigma2[j,k,l]=rinvgamma(1, rho+Nkl/2, gam+S2kl/2)
        }
      }
    }
    int N3 = 0;
    float S32 = 0;
    for (int i = 0; i < n; i++)
    {
      N3 = N3 + (kappajk[j][Zk[i] - 1] == 3);
      // N3=sum(kappajk[j,Zk]==3)
      // S32=sum((kappajk[j,Zk]==3)*(dat[,j])^2)
      S32 = S32 + (kappajk[j][Zk[i] - 1] == 3) * pow(dat[i][j], 2);
    }
    sigma2j3[j] = 1 / ran_gamma_knuth(rho + N3 / 2.0, 1 / (gam + S32 / 2));
    // sigma2j3[j]=rinvgamma(1, rho+N3/2, gam+S32/2)
  }

  for (int l = 0; l < 2; l++)
  {
    for (int k = 0; k < K; k++)
    {
      int Nkl = 0;
      float Skl = 0;
      for (int j = 0; j < p; j++)
      {
        Nkl = Nkl + (kappajk[j][k] == (l + 1));
        Skl = Skl + (kappajk[j][k] == (l + 1)) * Mu[j][k][l];
      }
      float Mean = (alpha * sigma2Muk[k][l] / sigma20Mu + Skl) / (sigma2Muk[k][l] / sigma20Mu + Nkl);
      float SD2 = 1 / (Nkl / sigma2Muk[k][l] + 1 / sigma20Mu);
      if (l == 0)
      {
        Muk[k][l] = r_lefttruncnorm(truncval, Mean, sqrt(SD2));
        // rtruncnorm(1, a=t, b=Inf, mean =Mean , sd = sqrt(SD2))
      }
      else if (l == 1)
      {
        Muk[k][l] = r_righttruncnorm(-truncval, Mean, sqrt(SD2));
        // rtruncnorm(1, a=-Inf, b=-t, mean =Mean , sd = sqrt(SD2))
      }
      float S2kl = 0;
      for (int j = 0; j < p; j++)
        S2kl = S2kl + (kappajk[j][k] == (l + 1)) * pow(Mu[j][k][l] - Muk[k][l], 2);
      sigma2Muk[k][l] = 1 / ran_gamma_knuth(rho + Nkl / 2.0, 1 / (gam + S2kl / 2));
      //       printf("%lf ",Muk[k][l]);
    }
  }
}

// ## Sample kappa
void SampleKappa(int K, int p, int n, int *Zk, int **kappajk, float ***Mu, float ***sigma2, float **TransMat, float **dat, float *sigma2j3)
{
  for (int k = 0; k < K; k++)
  {
    for (int j = 0; j < p; j++)
    {
      // logprobkappa=rep(0,3);
      float *logprobkappa = calloc(3, sizeof(float));
      for (int l = 0; l < 3; l++)
      {
        float MN = Mu[j][k][l];
        float SD;
        if (l < 2)
        {
          SD = sqrt(sigma2[j][k][l]);
        }
        else if (l == 2)
        {
          SD = sqrt(sigma2j3[j]);
        }
        if (j == 0)
        {
          logprobkappa[l] = log(TransMat[l][kappajk[1][k] - 1]);
          for (int i = 0; i < n; i++)
            logprobkappa[l] = logprobkappa[l] + ((Zk[i] == (k + 1)) * (-log(SD) - 0.5 * log(2 * M_PI) - pow(dat[i][j] - MN, 2) / (2 * pow(SD, 2))));
          // logprobkappa[l]=log(TransMat[l,kappajk[2,k]])+
          //  sum((Zk==k)*dnorm(dat[,j],mean =MN ,sd=SD,log = T))
        }
        else if (j < p - 1)
        {
          logprobkappa[l] = log(TransMat[kappajk[j - 1][k] - 1][l]) + log(TransMat[l][kappajk[j + 1][k] - 1]);
          for (int i = 0; i < n; i++)
            logprobkappa[l] = logprobkappa[l] + ((Zk[i] == (k + 1)) * (-log(SD) - 0.5 * log(2 * M_PI) - pow(dat[i][j] - MN, 2) / (2 * pow(SD, 2))));
          //            logprobkappa[l]=log(TransMat[kappajk[j-1,k],l])+log(TransMat[l,kappajk[j+1,k]])+
          //            sum((Zk==k)*dnorm(dat[,j],mean = MN,sd=SD,log = T))
        }
        else if (j == (p - 1))
        {
          logprobkappa[l] = log(TransMat[kappajk[j - 1][k] - 1][l]);
          for (int i = 0; i < n; i++)
            logprobkappa[l] = logprobkappa[l] + ((Zk[i] == (k + 1)) * (-log(SD) - 0.5 * log(2 * M_PI) - pow(dat[i][j] - MN, 2) / (2 * pow(SD, 2))));
          // logprobkappa[l]=log(TransMat[kappajk[j-1,k],l])+
          //  sum((Zk==k)*dnorm(dat[,j],mean = MN,sd=SD,log = T))
        }
      }
      double maxpr = max(3, logprobkappa);
      for (int l = 0; l < 3; l++)
      {
        logprobkappa[l] = exp(logprobkappa[l] - maxpr);
      }
      int *rho1 = (int *)calloc(3, sizeof(int));
      ran_multinomial(3, 1, logprobkappa, rho1);
      int kl;
      int *IDX = calloc(3, sizeof(int));
      findc(3, rho1, 0, IDX, &kl);
      kappajk[j][k] = 1 + IDX[0];
      free(IDX);
      free(rho1);
      //        kappajk[j,k]=sample(1:3,size=1,prob=exp(logprobkappa-max(logprobkappa)))
    }
  }
}

// #### Sample cluster probabilities
void ClusterProba(float *Omeg, int K, int n, double *alphaP, int *Zk)
{
  float *alphaPP = (float *)calloc(K, sizeof(float));
  for (int k = 0; k < K; k++)
  {
    for (int i = 0; i < n; i++)
      alphaPP[k] = alphaPP[k] + (Zk[i] == (k + 1));
    alphaPP[k] = alphaPP[k] + alphaP[k];
  }
  gsl_ran_dirichlet(K, alphaPP, Omeg);
  free(alphaPP);
}

void YPredict(float *Ypred, int *Zk, int pc, float **Cov, int n, int K, float *Beta, float *CoefY)
{
  for (int i = 0; i < n; i++)
  {
    float xb = 0;
    for (int l = 0; l < pc; l++)
    {
      xb += Cov[i][l] * Beta[l];
    }
    Ypred[i] = xb;
    for (int k = 0; k < K; k++)
    {
      if (Zk[i] == (k + 1))
      {
        Ypred[i] += CoefY[k];
      }
    }
  }
}

void SampleZ(int *Zk, int M, int n, int K, int *p, int ***kappajk, float ****sigma2, float ***dat, float ****Mu,
             float **sigma2j3, double *Y, float *CoefY, float sig2y, float *Omeg, int pred)
{
  // ## Sample Z

  for (int i = 0; i < n; i++)
  {
    // logprobZ=rep(0,K);
    float *logprobZ = (float *)calloc(K, sizeof(float));
    for (int k = 0; k < K; k++)
    {
      for (int m = 0; m < M; m++)
        for (int j = 0; j < p[m]; j++)
          logprobZ[k] = logprobZ[k] + ((kappajk[m][j][k] == 1) * (-0.5 * log(2 * M_PI) - log(sqrt(sigma2[m][j][k][0])) - pow(dat[m][i][j] - Mu[m][j][k][0], 2) / (2 * sigma2[m][j][k][0]))) +
                        ((kappajk[m][j][k] == 2) * (-0.5 * log(2 * M_PI) - log(sqrt(sigma2[m][j][k][1])) - pow(dat[m][i][j] - Mu[m][j][k][1], 2) / (2 * sigma2[m][j][k][1]))) +
                        ((kappajk[m][j][k] == 3) * (-0.5 * log(2 * M_PI) - log(sqrt(sigma2j3[m][j])) - pow(dat[m][i][j], 2) / (2 * sigma2j3[m][j])));
      logprobZ[k] += log(Omeg[k]);
      if (pred == 0)
        logprobZ[k] += -0.5 * pow(Y[i] - CoefY[k], 2) / sig2y;

      // logprobZ[k]=sum((kappajk[,k]==1)*dnorm(dat[i,],mean = Mu[,k,1],sd=sqrt(sigma2[,k,1]),log = T))+
      // sum((kappajk[,k]==2)*dnorm(dat[i,],mean = Mu[,k,2],sd=sqrt(sigma2[,k,2]),log=T))+
      // sum((kappajk[,k]==3)*dnorm(dat[i,],mean = 0,sd=sqrt(sigma2j3),log=T))+log(Omeg[k])
    }
    float maxZ = max(K, logprobZ);
    // float suml=0;
    for (int k = 0; k < K; k++)
    {
      logprobZ[k] = exp(logprobZ[k] - maxZ);
      // suml=suml+exp(logprobZ[k]-maxZ);
    }
    int *rho1 = (int *)calloc(K, sizeof(int));
    ran_multinomial(K, 1, logprobZ, rho1);
    int kl;
    int *IDX = calloc(K, sizeof(int));
    findc(K, rho1, 0, IDX, &kl);
    Zk[i] = 1 + IDX[0];
    free(IDX);
    free(rho1);
    free(logprobZ);
  } // End sample Z
}
void SampleFixedEffectsY(double *Y, int *Zk, int K, int n, float *CoefY, float *Beta, float **Cov, int pr)
{
  float *Yr = (float *)calloc(n, sizeof(float));

  for (int i = 0; i < n; i++)
  {
    for (int k = 0; k < K; k++)
    {
      if (Zk[i] == (k + 1))
      {
        Yr[i] = Y[i] - CoefY[k];
        // printf("\nYr=%d \n",i);
      }
    }
  }
  float *Matr = (float *)calloc(pr * pr, sizeof(float));
  float *No = (float *)calloc(pr, sizeof(float));

  for (int l = 0; l < pr; l++)
  {
    No[l] = GenerateNormal(1);
    for (int l1 = 0; l1 < l + 1; l1++)
    {
      for (int i = 0; i < n; i++)
      {
        Matr[l * pr + l1] += Cov[i][l] * Cov[i][l1];
      }
      if (l == l1)
      {
        Matr[l * pr + l1] += 1.0 / 100;
      }
      Matr[l1 * pr + l] = Matr[l * pr + l1];
      // printf("\nMatCov=%lf, %d, %d \n",Matr[l * pr + l1],l,l1);
      // printf("MatCov=%lf ",Matr[l * pr + l1]);
    }
    // printf("\n");
  }
  float *xb = (float *)calloc(pr, sizeof(float));
  for (int l = 0; l < pr; l++)
  {
    for (int i = 0; i < n; i++)
    {
      xb[l] += Cov[i][l] * Yr[i];
    }
  }
  float *L = cholesky(Matr, pr);
  float *BetaMean = (float *)calloc(pr, sizeof(float));
  LinearSystem(L, pr, BetaMean, xb);
  float *iLT = invertLT(L, pr);

  float *No1 = (float *)calloc(pr, sizeof(float));
  for (int l = 0; l < pr; l++)
  {
    for (int l1 = 0; l1 < pr; l1++)
      No1[l] += iLT[l * pr + l1] * No[l1];
    Beta[l] = BetaMean[l] + No1[l];
    // printf("\nBet11=%lf ",BetaMean[l]);
  }

  free(No);
  free(No1);
  free(L);
  free(iLT);
  free(Yr);
  free(BetaMean);
  free(xb);
}

double *ResiXy(double *Y, int n, float *Beta, float **Cov, int pc)
{
  double *Y1 = (double *)calloc(n, sizeof(double));
  for (int i = 0; i < n; i++)
  {
    float xb = 0;
    for (int l = 0; l < pc; l++)
    {
      xb += Cov[i][l] * Beta[l];
    }
    Y1[i] = Y[i] - xb;
  }
  return Y1;
}
void SampleCoefYSig2y(double *Y, float *CoefY, int *Zk, float *sig2y, float sig2a, int K, int n, float aa, float ba, int OutcomType)
{
  if (OutcomType == 2)
  {
    *sig2y = 1;
  }
  float sumX2 = 0;
  for (int k = 0; k < K; k++)
  {
    float meany = 0;
    int sizek = 0;
    for (int i = 0; i < n; i++)
    {
      if (Zk[i] == (k + 1))
      {
        meany += Y[i];
        sizek += 1;
      }
    }
    float sig2 = 1 / (1 / sig2a + sizek / (*sig2y));
    CoefY[k] = sig2 * meany / (*sig2y) + GenerateNormal(sqrt(sig2));
    if (sizek == 0)
    {
      CoefY[k] = 0;
    } // that is if the size of cluster k is empty
    for (int i = 0; i < n; i++)
    {
      if (Zk[i] == (k + 1))
        sumX2 += pow(Y[i] - CoefY[k], 2);
    }
  }
  if (OutcomType == 2)
  {
    *sig2y = 1;
  }
  else
  {
    *sig2y = 1 / ran_gamma_knuth(0.5 * n + aa, 1 / (0.5 * sumX2 + ba));
  }
  // printf("Sigma2=%f ",*sig2y);
}

void logposterior(int n, int M, int *p, int pc, float ***dat, double *Y, float **Cov,
                  float *Beta, float *CoefY, float sig2y, float ****Mu, float ***Muk,
                  float ***sigma2Muk, float ****sigma2, float **sigma2j3, int *Zk,
                  int ***kappajk, float ***TransMat, float *Omeg, double *alphaP,
                  float truncval, float alpha, float rho, float gam, double sigma20Mu,
                  double *delt, int K, double *logpost1, double *loglik1)
{
  float loglik = 0;
  float logpost = 0;
  // n=nrow(dat);p=ncol(dat);
  for (int i = 0; i < n; i++)
  {
    float xb = 0;
    for (int l = 0; l < pc; l++)
    {
      xb += Cov[i][l] * Beta[l];
    }
    for (int k = 0; k < K; k++)
    {
      if (Zk[i] == (k + 1))
      {
        for (int m = 0; m < M; m++)
          for (int j = 0; j < p[m]; j++)
            loglik = loglik + ((kappajk[m][j][k] == 1) * (-0.5 * log(2 * M_PI) - log(sqrt(sigma2[m][j][k][0])) - pow(dat[m][i][j] - Mu[m][j][k][0], 2) / (2 * sigma2[m][j][k][0]))) +
                     ((kappajk[m][j][k] == 2) * (-0.5 * log(2 * M_PI) - log(sqrt(sigma2[m][j][k][1])) - pow(dat[m][i][j] - Mu[m][j][k][1], 2) / (2 * sigma2[m][j][k][1]))) +
                     ((kappajk[m][j][k] == 3) * (-0.5 * log(2 * M_PI) - log(sqrt(sigma2j3[m][j])) - pow(dat[m][i][j], 2) / (2 * sigma2j3[m][j])));
        // loglik=loglik+sum((kappajk[,k]==1)*dnorm(dat[i,],mean = Mu[,k,1],sd=sqrt(sigma2[,k,1]),log = T))+
        // sum((kappajk[,k]==2)*dnorm(dat[i,],mean = Mu[,k,2],sd=sqrt(sigma2[,k,2]),log=T))+
        // sum((kappajk[,k]==3)*dnorm(dat[i,],mean = Mu[,k,3],sd=sqrt(sigma2j3),log=T))
        loglik = loglik - 0.5 * log(2 * M_PI) - log(sqrt(sig2y)) - 0.5 * pow(Y[i] - xb - CoefY[k], 2) / sig2y;
      }
    }
  }
  //  #logprior for Zk
  for (int k = 0; k < K; k++)
  {
    for (int i = 0; i < n; i++)
    {
      logpost = logpost + (Zk[i] == (k + 1)) * log(Omeg[k]);
      // printf("logpost=%lf ",logpost);
    }
    for (int m = 0; m < M; m++)
    {

      for (int j = 0; j < p[m]; j++)
      {
        for (int l = 0; l < 3; l++)
        {
          if (j == 0)
          {
            for (int m = 0; m < M; m++)
              logpost = logpost + (kappajk[m][0][k] == (l + 1)) * log(1.0 / 3);
            // printf("logpost=%lf ",(kappajk[0][k]==(l+1))*log(1/3));
          }
          else
          {

            logpost = logpost + (kappajk[m][j][k] == (l + 1)) * log(TransMat[m][kappajk[m][j - 1][k] - 1][l]);
          }
        }
      }
    }
  }
  // logpost=logpost+ddirichlet(Omeg,as.vector(alphaP))
  /*  printf("\n");
  for (int k=0;k<K;k++){
  printf("%lf,", Omeg[k]);
  }
   printf("\nAlphaP\n");
  for (int k=0;k<K;k++){
  printf("%lf,", alphaP[k]);
  }
  */
  // printf("logpost=%lf ",gsl_ran_dirichlet_lnpdf (K, alphaP, Omeg));
  logpost = logpost + gsl_ran_dirichlet_lnpdf(K, alphaP, Omeg);
  for (int r = 0; r < 3; r++)
  {
    for (int m = 0; m < M; m++)
      logpost = logpost + gsl_ran_dirichlet_lnpdf(3, delt, TransMat[m][r]);
    // logpost=logpost+ddirichlet(TransMat[r,],as.vector(delt))
  }

  for (int k = 0; k < K; k++)
  {
    for (int l = 0; l < 2; l++)
    {
      double z1 = (truncval - alpha) / sqrt(sigma20Mu);
      double z2 = (-truncval - alpha) / sqrt(sigma20Mu);
      for (int m = 0; m < M; m++)
      {
        double z = (Muk[m][k][l] - alpha) / sqrt(sigma20Mu);
        if (l == 0)
        {
          logpost = logpost + (-0.5 * log(2 * M_PI) - log(sqrt(sigma20Mu)) - z * z / 2) - log(.5 * erfc(z1 / sqrt(2)));

          // logpost=logpost+dtruncnorm(Muk[k,l], a=t, b=Inf, mean =alpha , sd = sqrt(sigma2Muk[k,l]))
        }
        else if (l == 1)
        {

          logpost = logpost + (-0.5 * log(2 * M_PI) - log(sqrt(sigma20Mu)) - z * z / 2) - log(.5 * (1 + erf(z2 / sqrt(2))));
          // logpost=logpost+dtruncnorm(Muk[k,l], a=-Inf, b=-t, mean =alpha , sd = sqrt(sigma2Muk[k,l]))
        }
        for (int j = 0; j < p[m]; j++)
        {
          logpost = logpost - 0.5 * log(2 * M_PI) - .5 * log(sigma2Muk[m][k][l]) - pow(Mu[m][j][k][l] - Muk[m][k][l], 2) / (2 * sigma2Muk[m][k][l]);
          // logpost=logpost+sum(dnorm(Mu[,k,l],mean = Muk[k,l],sd=sqrt(sigma2Muk[k,l]),log = T))
          logpost = logpost - (rho + 1) * log(sigma2[m][j][k][l]) - gam / sigma2[m][j][k][l];
        }
        // logpost=logpost+sum(dinvgamma(sigma2[,k,l], rho, gam))
        logpost = logpost - (rho + 1) * log(sigma2Muk[m][k][l]) - gam / sigma2Muk[m][k][l];
        // logpost=logpost+dinvgamma(sigma2Muk[k,l], rho, gam)
      }
    }
  }
  for (int m = 0; m < M; m++)
    for (int j = 0; j < p[m]; j++)
      logpost = logpost - (rho + 1) * log(sigma2j3[m][j]) - gam / sigma2j3[m][j];
  for (int j = 0; j < pc; j++)
    logpost = logpost - 0.5 * log(2 * M_PI) - .5 * log(100) - pow(Beta[j], 2) / (2 * 100);

  logpost = logpost - (rho + 1) * log(sig2y) - gam / sig2y;
  // logpost=logpost+sum(dinvgamma(sigma2j3, rho, gam))

  logpost = logpost + loglik;
  // printf("LL= %.2lf ", loglik); printf("LogPost= %.2lf ", logpost);
  *loglik1 = loglik;
  *logpost1 = logpost;

  // list(loglik=loglik,logpost=logpost+loglik)
}
