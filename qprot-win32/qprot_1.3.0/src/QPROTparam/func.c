#include "pep.h"

void countFunc(int *nA, int *nB, int *label, int *q) {
  int j;
  *nA = 0;
  *nB = 0;
  for(j=0;j<*q;j++) {
    if(label[j]==0) *nA = (*nA) + 1;
    else *nB = (*nB) + 1; 
  }
}

void initialize(PARAM *param, DATA *data, int *p, int *q, int *nA, int *nB, int *total) {
  int i,j;
  param->p = *p;
  param->q = *q;
  param->nA = *nA;
  param->nB = *nB;

  /* init_data(ratio, p, total); */
  /* assert(param->ratio = (double *) calloc(param->p, sizeof(double))); */
  
  assert(param->diff = (double **) calloc(param->p, sizeof(double *)));
  for(i=0;i<param->p;i++) assert(param->diff[i] = (double *) calloc(*total, sizeof(double)));

  assert(lambda = (double *) calloc(param->p * param->q, sizeof(double)));
  
  assert(mu0 = (double *) calloc(param->p, sizeof(double)));
  assert(mu1 = (double *) calloc(param->p, sizeof(double)));
  assert(sigmasq0 = (double *) calloc(param->p, sizeof(double)));

  assert(param->zstat = (double *) calloc(param->p, sizeof(double)));
  assert(param->FDR = (double *) calloc(param->p, sizeof(double)));

  assert(param->mu0sum = (double *) calloc(param->p, sizeof(double)));
  assert(param->mu1sum = (double *) calloc(param->p, sizeof(double)));
  assert(param->sigmasq0sum = (double *) calloc(param->p, sizeof(double)));

  for(i=0;i<param->p;i++) {
    param->mu0sum[i] = 0.0;
    param->mu1sum[i] = 0.0;
    param->sigmasq0sum[i] = 0.0;
  }

  /* set DPs */
  m_mu0 = 0.0;
  v_mu0 = 100.0;
  for(i=0;i<param->p;i++) mu0[i] = gsl_ran_gaussian(r, 1.0) + m_mu0; 

  m_mu1 = 0.0;
  v_mu1 = 100.0;
  for(i=0;i<param->p;i++) mu1[i] = gsl_ran_gaussian(r, 1.0) + m_mu1; 

  a_sigmasq0 = 1.0;
  b_sigmasq0 = 0.1;
  for(i=0;i<param->p;i++) {
    sigmasq0[i] = 0.1;
  }

  for(i=0;i<*p;i++) {
    for(j=0;j<*q;j++) {
      lambda[i*data->ncol+j] = data->label[j] == 1 ? mu0[i] + mu1[i] : mu0[i] ;
    }
  }

}

void free_param(PARAM *param) {
  free(lambda);
  free(mu0);
  free(mu1);
  free(sigmasq0);
  free(sigmasq1);
  free(param->mu0sum);
  free(param->mu1sum);
  free(param->sigmasq0sum);
  free(param->sigmasq1sum);
}

