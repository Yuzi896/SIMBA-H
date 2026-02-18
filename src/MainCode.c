#include <string.h>
#include <time.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "utils.h"
#include "distribution.h"
#include "GammaFunction.hpp"
#include <R_ext/Print.h>

// #include "functions.h"
// To format a selection: CMD+K, CMD+F
// #include "MCMC.h"
void MCMCfunction(int *Method, int *OutcomeType, int *M1, int *p1, int *nn, double *Data1,
                  double *Cov1, int *pc, double *Youtcome, int *KK, double *truncatedvalue,
                  double *rhosigma, double *gamsigma, double *alphaP, double *delt,
                  double *sigma20Mu1, int *Nsample1, int *burnin1, int *seed1, double *probZMean1,
                  double *probkappajk1, double *TransMatMean1, double *MuMean1,
                  double *CoefYsample, double *sig2yMean, double *logpostsample,
                  double *Betasample, double *DIC, int *nnTest, double *DataTest1,
                  double *CovTest1, double *YpredTest, int *Ztest)
// MCMCfunction<-function(Method="HMM",Data=Data,K=K,t=0.01,alpha=0,rho=.1,gam=.1,alphaP=rep(1,K),
//                      delt=rep(1,3),sigma20Mu=1000,
//                    Nsample=1000, burnin=300,seed=1)
// M1=number of platforms
{
  //  set.seed(seed)
  // Manage some data
  int seed = seed1[0];
  int M = M1[0];
  srand(seed);
  int n = nn[0];
  int nTest = nnTest[0];
  int K = KK[0];
  int pc1 = pc[0]; // # of covariates
  double *Y = (double *)calloc(n, sizeof(double));
  float *Ypred = (float *)calloc(nTest, sizeof(float));
  double *Cens = (double *)calloc(n, sizeof(double));
  if (OutcomeType[0] == 1)
  { // Continous outcome
    for (int i = 0; i < n; i++)
      Y[i] = Youtcome[i];
  }
  else if (OutcomeType[0] == 2)
  { // Binary
    for (int i = 0; i < n; i++)
    {
      if (Youtcome[i] == 1)
        Y[i] = r_lefttruncnorm(0, 0, 1.0);
      else
        Y[i] = r_righttruncnorm(0, 0, 1.0);
    }
  }
  else if (OutcomeType[0] == 3)
  { // Right censored survival
    for (int i = 0; i < n; i++)
    {
      Y[i] = Youtcome[i];
      Cens[i] = Youtcome[i + n];
    }
  }

  int l = 0;
  int ltt = 0;
  float *Beta = (float *)calloc(pc1, sizeof(float));
  float **CovTest = (float **)calloc(nTest, sizeof(float *));
  float **Cov = (float **)calloc(n, sizeof(float *));
  for (int i = 0; i < n; i++)
  {
    Cov[i] = (float *)calloc(pc1, sizeof(float));
    for (int j = 0; j < pc1; j++)
    {
      Cov[i][j] = Cov1[l];
      l++;
    }
  }

  for (int i = 0; i < nTest; i++)
  {
    CovTest[i] = (float *)calloc(pc1, sizeof(float));
    for (int j = 0; j < pc1; j++)
    {
      CovTest[i][j] = CovTest1[ltt];
      ltt++;
    }
  }
  l = 0;
  ltt = 0;
  float ***dat = (float ***)calloc(M, sizeof(float **));
  float ***datTest = (float ***)calloc(M, sizeof(float **));
  for (int m = 0; m < M; m++)
  {
    int p = p1[m];
    dat[m] = (float **)calloc(n, sizeof(float *));
    // float **dat = (float **)calloc(n, sizeof(float *));
    for (int i = 0; i < n; i++)
    {
      dat[m][i] = (float *)calloc(p, sizeof(float));
      for (int j = 0; j < p; j++)
      {
        dat[m][i][j] = Data1[l];
        l++;
      }
    }
    datTest[m] = (float **)calloc(nTest, sizeof(float *));
    for (int i = 0; i < nTest; i++)
    {
      datTest[m][i] = (float *)calloc(p, sizeof(float));
      for (int j = 0; j < p; j++)
      {
        datTest[m][i][j] = DataTest1[ltt];
        ltt++;
      }
    }
  }
  double rho = rhosigma[0];
  double gam = gamsigma[0];
  int Method1 = Method[0];
  float sigma20Mu = sigma20Mu1[0];
  float truncval = truncatedvalue[0];
  //  ### Initialization
  float ***TransMat = (float ***)calloc(M, sizeof(float **));
  for (int m = 0; m < M; m++)
  {
    TransMat[m] = (float **)calloc(3, sizeof(float *));
    for (int l = 0; l < 3; l++)
    {
      TransMat[m][l] = (float *)calloc(3, sizeof(float));
      for (int l1 = 0; l1 < 3; l1++)
      {
        TransMat[m][l][l1] = 1.0 / 3;
      }
    }
  }
  // TransMat=matrix(rep(c(1/3,1/3,1/3),3),3,byrow = T)
  // Zk=rep(0,n)
  int *Zktest = (int *)calloc(nTest, sizeof(int));
  int *Zk = (int *)calloc(n, sizeof(int));
  double alpha = 0; // Prior mean of Mukl
                    // dat=Data

  for (int i = 0; i < n; i++)
  {
    float *probrho = (float *)calloc(K, sizeof(float));
    for (int k = 0; k < K; k++)
    {
      probrho[k] = 1.0 / K;
    }
    int *rho1 = (int *)calloc(K, sizeof(int));
    ran_multinomial(K, 1, probrho, rho1);
    int kl;
    int *IDX = calloc(K, sizeof(int));
    findc(K, rho1, 0, IDX, &kl);
    Zk[i] = 1 + IDX[0];
    // Rprintf("%d ",Zk[i]);
    free(rho1);
    free(probrho);
  }

  // for (i in 1:n){
  //   Zk[i]=sample(1:K,size=1)}
  //  kappajk=matrix(0,p,K)
  float *prob = (float *)calloc(3, sizeof(float));
  for (int l = 0; l < 3; l++)
  {
    prob[l] = 1.0 / 3;
  }

  int ***kappajk = (int ***)calloc(M, sizeof(int **));
  for (int m = 0; m < M; m++)
  {
    int p = p1[m];
    kappajk[m] = (int **)calloc(p, sizeof(int *));
    for (int j = 0; j < p; j++)
    {
      kappajk[m][j] = (int *)calloc(K, sizeof(int));
      for (int k = 0; k < K; k++)
      {
        int *kap = (int *)calloc(3, sizeof(int));
        ran_multinomial(3, 1, prob, kap);
        int kl;
        int *IDX = calloc(3, sizeof(int));
        // Rprintf("kap \n");
        // for (int l=0;l<3;l++)
        // Rprintf("%d ",kap[l]);

        findc(3, kap, 0, IDX, &kl);
        kappajk[m][j][k] = 1 + IDX[0];
        // Rprintf("kappaInit \n");
        // Rprintf("%d, ",kappajk[j][k]);
        free(kap);
        free(IDX);
      }
    }
  }
  free(prob);
  // for (k in 1:K){
  // for (j in 1:p){
  //   kappajk [j,k]=sample(1:3,size=1,prob=c(1,1,1))
  //  }
  //  }
  // Muk=matrix(0,K,3)
  float ***Muk = (float ***)calloc(M, sizeof(float **));
  float ***sigma2Muk = (float ***)calloc(M, sizeof(float **));
  for (int m = 0; m < M; m++)
  {
    Muk[m] = (float **)calloc(K, sizeof(float *));
    sigma2Muk[m] = (float **)calloc(K, sizeof(float *));
    for (int k = 0; k < K; k++)
    {
      Muk[m][k] = (float *)calloc(3, sizeof(float));
      sigma2Muk[m][k] = (float *)calloc(2, sizeof(float));
      for (int l = 0; l < 2; l++)
        sigma2Muk[m][k][l] = 1.0;
    }
  }
  float sig2y = 1;
  // sigma2=array(1,dim = c(p,K,2))
  float *CoefY = (float *)calloc(K, sizeof(float));
  float *Omeg = (float *)calloc(K, sizeof(float));
  float ****sigma2 = (float ****)calloc(M, sizeof(float ***));
  float ****Mu = (float ****)calloc(M, sizeof(float ***));
  float ****probkappajk = (float ****)calloc(M, sizeof(float ***));

  for (int m = 0; m < M; m++)
  {
    int p = p1[m];
    sigma2[m] = (float ***)calloc(p, sizeof(float **));
    Mu[m] = (float ***)calloc(p, sizeof(float **));
    probkappajk[m] = (float ***)calloc(p, sizeof(float **));

    for (int j = 0; j < p; j++)
    {
      sigma2[m][j] = (float **)calloc(K, sizeof(float *));
      Mu[m][j] = (float **)calloc(K, sizeof(float *));
      probkappajk[m][j] = (float **)calloc(K, sizeof(float *));
      for (int k = 0; k < K; k++)
      {
        sigma2[m][j][k] = (float *)calloc(2, sizeof(float));
        Mu[m][j][k] = (float *)calloc(3, sizeof(float));
        probkappajk[m][j][k] = (float *)calloc(3, sizeof(float));
        for (int l = 0; l < 2; l++)
        {
          sigma2[m][j][k][l] = 1.0;
        }
      }
    }
  }
  // sigma2j3=rep(1,p);
  float **sigma2j3 = (float **)calloc(M, sizeof(float *));
  for (int m = 0; m < M; m++)
  {
    int p = p1[m];
    sigma2j3[m] = (float *)calloc(p, sizeof(float));
    for (int j = 0; j < p; j++)
      sigma2j3[m][j] = 1.0;
  }
  // Mu=array(0,dim=c(p,K,3))
  // sigma2Muk=matrix(1,K,2)

  // We set output
  // probZMean=matrix(0,n,K);
  float **probZMean = (float **)calloc(n, sizeof(float *));
  for (int i = 0; i < n; i++)
  {
    probZMean[i] = (float *)calloc(K, sizeof(float));
  }

  // probkappajk=array(0, dim=c(p,K,3));
  //  TransMatMean=matrix(0,3,3)

  float ***TransMatMean = (float ***)calloc(M, sizeof(float **));
  for (int m = 0; m < M; m++)
  {
    TransMatMean[m] = (float **)calloc(3, sizeof(float *));
    for (int l = 0; l < 3; l++)
    {
      TransMatMean[m][l] = (float *)calloc(3, sizeof(float));
    }
  }

  // MuMean=matrix(0,K,3)
  double ***MuMean = (double ***)calloc(M, sizeof(double **));
  double ***SigmaMuMean = (double ***)calloc(M, sizeof(double **));
  for (int m = 0; m < M; m++)
  {
    MuMean[m] = (double **)calloc(K, sizeof(double *));
    SigmaMuMean[m] = (double **)calloc(K, sizeof(double *));
    for (int k = 0; k < K; k++)
    {
      MuMean[m][k] = (double *)calloc(3, sizeof(double));
      SigmaMuMean[m][k] = (double *)calloc(2, sizeof(double));
    }
  }
  // SigmaMuMean=matrix(0,K,2);

  // logpostsample=NULL;
  int Nsample = Nsample1[0];
  int burnin = burnin1[0];
  float *loglik = (float *)calloc(Nsample + burnin, sizeof(float));

  // loglik=NULL
  // Rprintf("\n");
  int ll = 0;
  int lll = 0;
  int lt = 0;
  int lp = 0;
  for (int tt = 0; tt < (Nsample + burnin); tt++)
  {
    // Rprintf("Mcmc iteration: %d\n", tt);
    int printNumber = (Nsample + burnin) / 10;
    if ((tt) % printNumber == 1)
    {
      // #print(paste("Number of mcmc iterations is ", tt))
      Rprintf("Mcmc iteration: %d\n", tt);
      // Rprintf("\rMcmc iteration: %d", tt);
      // cat(paste("\rmcmc iteration", tt))
    }

    if (OutcomeType[0] == 2)
    { // Binary
      SampleLatentYBinary(pc1, K, n, Y, Youtcome, Cov, CoefY, Zk, Beta);
    }
    else if (OutcomeType[0] == 3)
    { // Right censored survival
      SampleLatentYRightCensored(pc1, K, n, Y, Youtcome, Cens,
                                 Cov, CoefY, Zk, Beta, sig2y);
    }
    if (pc1 != 0)
    {
      // Sample efects of covariates
      SampleFixedEffectsY(Y, Zk, K, n, CoefY, Beta, Cov, pc1);
      if (tt >= burnin)
      {
        for (int l = 0; l < pc1; l++)
        {
          Betasample[ll] = Beta[l];
          // Rprintf("\nBet=%lf ",Beta[l]);
          ll++;
        }
      }
    }

    double *Y1 = ResiXy(Y, n, Beta, Cov, pc1);
    // Sample parameterers response
    double sig2a = 10; // variance of coefY
    SampleCoefYSig2y(Y1, CoefY, Zk, &sig2y, sig2a, K, n, rho, gam, OutcomeType[0]);
    if (tt > burnin)
    {
      for (int k = 0; k < K; k++)
      {
        CoefYsample[lll] = CoefY[k];
        lll++;
      }

      sig2yMean[0] = sig2yMean[0] + sig2y / Nsample;
    }
    // if (Method!="HMM"){
    for (int m = 0; m < M; m++)
    {
      int p = p1[m];
      // Sample Transition matrix
      SampleTransitionMatrice(K, p, Method1, TransMat[m], kappajk[m], delt);
      if (tt >= burnin)
      {
        // TransMatMean=TransMatMean+TransMat/Nsample
        for (int l = 0; l < 3; l++)
          for (int l1 = 0; l1 < 3; l1++)
            TransMatMean[m][l][l1] = TransMatMean[m][l][l1] + TransMat[m][l][l1] / Nsample;
      }
      //    ### Sample MU and Sigma
      SampleMuSigma(K, p, n, Zk, kappajk[m], Mu[m], Muk[m], sigma2[m], sigma2Muk[m], TransMat[m], dat[m], sigma2j3[m], rho, gam, alpha, sigma20Mu, truncval);

      if (tt > burnin)
      {
        for (int l = 0; l < 2; l++)
        {
          for (int k = 0; k < K; k++)
          {
            MuMean[m][k][l] = MuMean[m][k][l] + Muk[m][k][l] / Nsample;
            SigmaMuMean[m][k][l] = SigmaMuMean[m][k][l] + sigma2Muk[m][k][l] / Nsample;
          }
        }
        // MuMean=MuMean+Muk/Nsample;SigmaMuMean=SigmaMuMean+sigma2Muk/Nsample
      }

      // ### Sample sig2 variance of the noise

      // ### Sample latent kappa
      // ## Sample kappa

      SampleKappa(K, p, n, Zk, kappajk[m], Mu[m], sigma2[m], TransMat[m], dat[m], sigma2j3[m]);
      for (int k = 0; k < K; k++)
      {
        for (int j = 0; j < p; j++)
        {
          //        kappajk[j,k]=sample(1:3,size=1,prob=exp(logprobkappa-max(logprobkappa)))
          if (tt >= burnin)
          {
            for (int l = 0; l < 3; l++)
            {
              probkappajk[m][j][k][l] = probkappajk[m][j][k][l] + (kappajk[m][j][k] == (l + 1));
            }
          }
        }
      }

    } // end loop on m

    // #### Sample cluster probabilities
    ClusterProba(Omeg, K, n, alphaP, Zk);
    // Omeg=rdirichlet(1,as.vector(alphaP+sapply(1:K,function(x) sum(Zk==x))))

    // ## Sample Z
    SampleZ(Zk, M, n, K, p1, kappajk, sigma2, dat, Mu, sigma2j3, Y1, CoefY, sig2y, Omeg, 0);

    // Sample Z on the test data
    SampleZ(Zktest, M, nTest, K, p1, kappajk, sigma2, datTest, Mu, sigma2j3, Y1, CoefY, sig2y, Omeg, 1);
    // Predict Y on test data

    YPredict(Ypred, Zktest, pc1, CovTest, nTest, K, Beta, CoefY);

    free(Y1);

    for (int i = 0; i < nTest; i++)
    {

      if (tt >= burnin)
      {
        Ztest[lt] = Zktest[i];
        lt++;
        YpredTest[lp] = Ypred[i];
        lp++;
      }
    }

    for (int i = 0; i < n; i++)
    {

      if (tt >= burnin)
      {
        for (int k = 0; k < K; k++)
        {
          probZMean[i][k] = probZMean[i][k] + (Zk[i] == (k + 1));
        }
      }
    }
    //    LogP=logpost(dat,Mu,Muk,sigma2Muk,sigma2,sigma2j3,Zk,kappajk,TransMat,Omeg,alphaP,t,alpha,rho,gam,delt,K);

    //  logpostsample[tt]=LogP$logpost
    // loglik[tt]=LogP$loglik
    double logs = 0;
    double lk = 0;

    logposterior(n, M, p1, pc1, dat, Y, Cov, Beta, CoefY, sig2y, Mu, Muk,
    sigma2Muk, sigma2, sigma2j3, Zk, kappajk, TransMat, Omeg, alphaP, truncval,
    alpha, rho, gam, sigma20Mu, delt, K, &logs, &lk);
    logpostsample[tt] = logs;
    loglik[tt] = lk;
    // Rprintf("logpostsample= %.2lf ", logpostsample[tt]);
    // void logposterior(int n, int p,float ** dat,float ***Mu,float **Muk,float **sigma2Muk, float *** sigma2,float *sigma2j3,int *Zk, int ** kappajk,float ** TransMat,float* Omeg,float *alphaP,float truncval,float alpha,float rho,float gam,double sigma20Mu,float* delt,int K, float * logpost1,float *loglik1)
  } // end of loop tt

  // We get outputs here as vectors...

  ll = 0;
  for (int i = 0; i < n; i++)
  {
    // Rprintf("\n%d,", i);
    for (int k = 0; k < K; k++)
    {
      // Rprintf(" %.2lf ", probZMean[i][k] / Nsample);
      probZMean1[ll] = probZMean[i][k] / Nsample;
      ll++;
    }
  }

  ll = 0;
  for (int m = 0; m < M; m++)
  {
    int p = p1[m];
    // Rprintf("\nDataset: %d\n", m);
    for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < p; j++)
      {
        // if (j < 10)
        //  Rprintf("\n\n\n");

        // Rprintf("\n %d,",k);
        for (int l = 0; l < 3; l++)
        {
          // if (j < 10)
          //  Rprintf(" %.2lf ", probkappajk[m][j][k][l] / Nsample);
          probkappajk1[ll] = probkappajk[m][j][k][l] / Nsample;
          ll++;
        }
      }
    }
  }

  ll = 0;
  // Rprintf("\nTransMatMeans\n");
  for (int m = 0; m < M; m++)
  {
    // Rprintf("\nDataset: %d\n", m);
    for (int l = 0; l < 3; l++)
    {
      // Rprintf("\n");
      for (int l1 = 0; l1 < 3; l1++)
      {
        TransMatMean1[ll] = TransMatMean[m][l][l1];
        ll++;
        // Rprintf(" %.2lf ", TransMatMean[m][l][l1]);
      }
    }
  }

  ll = 0;
  // Rprintf("\nMeans\n");
  for (int m = 0; m < M; m++)
  {
    // Rprintf("\nDataset: %d\n", m);
    for (int k = 0; k < K; k++)
    {
      // Rprintf("\n %d,", k);
      for (int l = 0; l < 2; l++)
      {
        MuMean1[ll] = MuMean[m][k][l];
        // Rprintf(" %.2lf ", MuMean[m][k][l]);
        ll++;
      }
    }
  }

  // Rprintf("\nYMeans\n");
  // for (int k = 0; k < K; k++)
  //{
  //  Rprintf("\n %d,", k);
  //  Rprintf(" %.2lf ", CoefYMean[k]);
  // }
  // cat("\n")
  // ttmax=which.max(logpostsample)
  int ttmax = maxIdx(Nsample + burnin, logpostsample);
  DIC[0] = -4 * Mean(Nsample + burnin, loglik) + 2 * loglik[ttmax];

  // Rprintf("\nDIC=%lf ",DIC[0]);

  // list(probZMean=probZMean/(Nsample),probkappajk=probkappajk/Nsample,TransMatMean=TransMatMean,MuMean=MuMean,
  //    logpostsample=logpostsample,DIC=-4*mean(loglik)+2*loglik[ttmax])

  Rprintf("\n");
}
