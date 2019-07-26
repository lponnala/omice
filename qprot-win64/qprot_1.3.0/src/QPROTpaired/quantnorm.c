#include "pep.h"

void normalizeData(DATA *data) {
  int i,j,k,ct;
  double dist_low, dist_high, dist_sum;
  double **percentile;
  double *mod_percentile;
  double *stack; 
  double *tmp;
  
  assert(percentile = (double **) calloc(data->ncol, sizeof(double *)));
  for(j=0;j<data->ncol;j++) assert(percentile[j] = (double *) calloc(101, sizeof(double)));

  assert(mod_percentile = (double *) calloc(101, sizeof(double)));

  assert(stack = (double *) calloc(data->nrow, sizeof(double)));
  assert(tmp = (double *) calloc(data->ncol, sizeof(double)));

  // get percentiles in each data: copy onto stack, sort, and get percentiles
  for(j=0;j<data->ncol;j++) {
    ct = 0;
    for(i=0;i<data->nrow;i++) {
      if(data->dat[i*data->ncol+j] > 0) {
        stack[ct] = data->d[i*data->ncol+j];
        ct++;
      }
    }
    gsl_sort(stack, 1, ct);
    for(k=0;k<=100;k++) percentile[j][k] = stack[k * (ct-1) / 100];  
  }

  // equalize percentiles by median
  for(k=0;k<=100;k++) {
    for(j=0;j<data->ncol;j++) tmp[j] = percentile[j][k];
    mod_percentile[k] = vec_med(tmp, data->ncol);
  }

  // for each data point, find the nearest two percentiles and interpolate from mod_percentile
  for(j=0;j<data->ncol;j++) {
    for(i=0;i<data->nrow;i++) {
      if(data->dat[i*data->ncol+j] > 0) {
        k = 0;
        while(data->d[i*data->ncol+j] >= percentile[j][k] && k<= 100) k++;
        if(k > 100) k = 100;
        if(k < 1) k = 1;
        dist_low = data->d[i*data->ncol+j] - percentile[j][k-1];
        dist_high = percentile[j][k] - data->d[i*data->ncol+j];
        if(dist_low == 0.0 && dist_high == 0.0) data->d[i*data->ncol+j] = mod_percentile[k];
        else {   
          dist_sum = dist_low + dist_high;
          data->d[i*data->ncol+j] = (dist_high / dist_sum) * mod_percentile[k-1] + (dist_low / dist_sum) * mod_percentile[k];
        }
      }
    }
  }

}


