
/*
    Copyright (C) <2011>  <Hyungwon Choi, Damian Fermin>
    For troubleshooting, contact hyung_won_choi@nuhs.edu.sg.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You can obtain a copy of the GNU General Public License from
    <http://www.gnu.org/licenses/>.
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>


#define _MAX_COMP_ 10

#define _SKIP_ 10
#define _PRINT_FREQ_ 100
#define _tiny_ 1e-100


gsl_rng *r;

char **X;
int np;
int nq;

double bandwidth;
int ntick;
double *tickMarks;
double *f;
double *f1;
double *f0;

double *df;
double *df1;
double *df0;

int PRINTDENSITY;
int THEORETICAL;

double thres;

/* DP parameters */

double pi_true;
double nullmean;
double nullvar;

double *zstat;

double *fdr;
double *FDRup;
double *FDRdown;


/*************/
/* functions */
/*************/
int nrow(FILE *fp);
int ncol(FILE *fp);
int newlinechar(char *buf, int k);

/****************************************************************************************/
int read_matrix_data(FILE *fpm);
void init_param(void );  
void check_mean_shift(void );
void getBandwidth(void );

void computeProportion(void );

void fit_FDR(void );
void computeFDR(void );

/****************************************************************************************/
double vec_sum(const double *vec, int len);
double vec_max(const double *vec, int len);
double vec_max(const double *vec, int len);
double vec_min(const double *vec, int len);
double vec_mean(const double *vec, int len);
double vec_var(const double *vec, int len);
double vec_med(const double *vec, int len);
double vec_mad(const double *vec, int len);
int ranMultinom(double *p, int K);
void vec_var2(const double *vec, int len, double *mean, double *var);

/****************************************************************************************/

double logit(double p);
double invlogit(double x);


double gaussian_pdf(double x, double mu, double sigmasq);
double log_gaussian_pdf(double x, double mu, double sigmasq);
void print_result(FILE *fp);
void print_density(FILE *fp);

void evaluateDensities(void );

double fmin(double a, double b);


