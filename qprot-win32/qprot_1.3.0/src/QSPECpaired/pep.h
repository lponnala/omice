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

#define _MAX_PROT_ 10000
#define Max_buf 10000
#define Max_ID 1000
#define _MAX_SAMPLE_ 100

#define _MAX_COMP_ 50
#define _MAX_COMP2_ 50

#define _SKIP_ 20
#define _PRINT_FREQ_ 1000
#define _tiny_ 1e-100

typedef struct tagDATA {
  double scale;
  char **experiment;
  char **name;
  int *label;
  double *d;  
  double *dat;  
  double *len;
  double *sum;
  int nrow; 
  int ncol; 
  int npair;
} DATA;

typedef struct tagPARAM{
  int p;
  int q;
  int nA;
  int nB;

  double *ratio;
  double **diff;

  double **mu0sum;
  double *mu1sum;

  double *zstat;
  double *FDR;
  double pi[101];
  double mnosig[101];
  double vnosig[101];
} PARAM;


/* GLOBAL VARIABLES */
double **mu0;
double *mu1;

double *rowsum;
double *allzero0;
double *allzero1;

double *lambda;

double m_mu0;
double v_mu0;

double m_mu1;
double v_mu1;

gsl_rng *r;


int NORMALIZE;
double fdr_center;
double tmean;
double thres;

double elapsed;



/*************/
/* functions */
/*************/
int nrow(FILE *fp);
int ncol(FILE *fp);
int newlinechar(char *buf, int k);

/****************************************************************************************/
void init_data(DATA *data, int *p, int *q);
void free_data(DATA *data);
void free_param(PARAM *param);
int read_matrix_data(FILE *fpm, DATA *data, int *p, int *q);
void countFunc(int *nA, int *nB, int *label, int *q);

void normalizeData(DATA *data);

void update_mu0(PARAM *param, DATA *data);
void update_mu1(PARAM *param, DATA *data);
void update_mu_hyper(PARAM *param, DATA *data);


/****************************************************************************************/
void initialize(PARAM *param, DATA *data, int *p, int *q, int *nA, int *nB, int *total);

/****************************************************************************************/
int find_index(double mu, double sigmasq, double a);
double log_poisson(double x, double mu);
double log_gaussian(double d, double mu, double sigmasq);

void calcLambdaAll(DATA *data, PARAM *param);
void calcLambdaRow(DATA *data, PARAM *param, int r);
void calcLambdaCol(DATA *data, PARAM *param, int c);
void calcLambdaElement(DATA *data, PARAM *param, int r, int c);

void runBurnIn(DATA *data, PARAM *param, int *burn); 
/* type=-1 for full model, type=-2 for RR model, type=k, for R model w/o factor k */
void runGibbs(DATA *data, PARAM *param, int *iter);
void output(FILE *fp, DATA *data, PARAM *param, int *p, int *q, int *nA, int *nB, int *total);

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

/****************************************************************************************/


void fdr_master(DATA *data, PARAM *param) ;
void EM_FDR(PARAM *param, double *zstat, double xcap_neg, double xcap_pos, int ct);
void compute_FDR(DATA *data, PARAM *param, double xcap_neg, double xcap_pos);
