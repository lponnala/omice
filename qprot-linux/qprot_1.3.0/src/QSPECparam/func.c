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
  double tmp;

  param->p = *p;
  param->q = *q;
  param->nA = *nA;
  param->nB = *nB;

  assert(param->diff = (double **) calloc(param->p, sizeof(double *)));
  for(i=0;i<param->p;i++) assert(param->diff[i] = (double *) calloc(*total, sizeof(double)));

  assert(rowsum = (double *) calloc(param->p, sizeof(double)));
  assert(allzero0 = (double *) calloc(param->p, sizeof(double)));
  assert(allzero1 = (double *) calloc(param->p, sizeof(double)));

  assert(lambda = (double *) calloc(param->p * param->q, sizeof(double)));
  
  assert(mu0 = (double *) calloc(param->p, sizeof(double)));
  assert(mu1 = (double *) calloc(param->p, sizeof(double)));

  assert(param->mu0sum = (double *) calloc(param->p, sizeof(double)));
  assert(param->mu1sum = (double *) calloc(param->p, sizeof(double)));

  assert(param->zstat = (double *) calloc(param->p, sizeof(double)));
  assert(param->FDR = (double *) calloc(param->p, sizeof(double)));

  for(i=0;i<param->p;i++) {
    param->mu0sum[i] = 0.0;
    param->mu1sum[i] = 0.0;
  }

  /* set DPs */
  m_mu0 = 0.0;
  v_mu0 = 100.0;
  for(i=0;i<param->p;i++) mu0[i] = gsl_ran_gaussian(r, 1.0) + m_mu0;

  m_mu1 = 0.0;
  v_mu1 = 100.0;
  for(i=0;i<param->p;i++) mu1[i] = gsl_ran_gaussian(r, 1.0) + m_mu1;

  for(i=0;i<*p;i++) {
    rowsum[i] = 0.0;
    for(j=0;j<*q;j++) {
      rowsum[i] += data->d[i];
    }
  }

  for(i=0;i<*p;i++) {
    for(j=0;j<*q;j++) {
      lambda[i*data->ncol+j] = data->label[j] == 1 ? mu0[i] + mu1[i] : mu0[i];
      if(NORMALIZE) lambda[i*data->ncol+j] += log(data->len[i]) + log(data->sum[j]);
      lambda[i*data->ncol+j] = exp(lambda[i*data->ncol+j]);
    }
  }

  tmp = 0.0;
  for(j=0;j<*q;j++) {
    if(data->label[j] == 0) tmp += 1.0;
  }
  for(i=0;i<*p;i++) {
    allzero0[i] = 0.0;
    for(j=0;j<*q;j++) {
      if(data->label[j] == 0 && data->d[i*data->ncol+j] > 0) allzero0[i] += 1.0;
    }
    allzero0[i] /= tmp;
  }

  tmp = 0.0;
  for(j=0;j<*q;j++) {
    if(data->label[j] == 1) tmp += 1.0;
  }
  for(i=0;i<*p;i++) {
    allzero1[i] = 0.0;
    for(j=0;j<*q;j++) {
      if(data->label[j] == 1 && data->d[i*data->ncol+j] > 0) allzero1[i] += 1.0;
    }
    allzero1[i] /= tmp;
  }


}

void free_param(PARAM *param) {
  free(lambda);
  free(mu0);
  free(mu1);
  free(param->mu0sum);
  free(param->mu1sum);
}




