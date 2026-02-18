#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "GammaFunction.hpp"
#include "utils.h"
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define INT_MAX 2147483647

float LGammaF(float x);
void gsl_ran_dirichlet(int K, const float alpha[], float theta[]);
static void ran_dirichlet_small(int K, const float alpha[], float theta[]);
double gsl_ran_dirichlet_lnpdf(int K, const double alpha[], const float theta[]);
double ran_beta(const double a, const double b);
unsigned int ran_binomial(double p, unsigned int n);
void ran_multinomial(const int K, const unsigned int N, const float p[], int n[]);
float GenerateNormal(const float sigma);
static double gamma_large(const double a);
static double gamma_frac(const double a);
double ran_gamma_int(const unsigned int a);
double ran_gamma_knuth(const double a, const double b);
double rexponential(const double mu);
double ers_a_inf(double a);
double nrs_a_inf(double a);
double r_lefttruncnorm(double a, double mean, double sd);
double r_righttruncnorm(double b, double mean, double sd);
