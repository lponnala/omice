#include "pep.h"

int find_index(double mu, double sigmasq, double a) {
  int res;
  double tmp = (a-mu) / sqrt(sigmasq);
  if(tmp > 5.0) res = 10000;
  else if(tmp < -5.0) res = 0;
  else res = 5000 + ((int) (tmp * 1000.0));
  return res;
}

double log_poisson(double x, double mu) {
  double res;
  double xx;
  xx = x >= 1.0 ? x : 1.0;
  res = -mu + xx * log(mu) - gsl_sf_lngamma(xx + 1.0);
  return res;
}

double log_gaussian(double d, double mu, double sigmasq) {
  double out;
  out = -.5 * pow(d-mu,2.0) / sigmasq - .5 * log(2.0 * M_PI * sigmasq);
  return out;
}

void calcLambdaAll(DATA *data, PARAM *param) {
  int i, j;
  for(i=0;i<param->p;i++) {
    for(j=0;j<param->q;j++) {
      lambda[i*data->ncol+j] = data->label[j] == 1 ? mu0[i] + mu1[i] : mu0[i];
      if(NORMALIZE) lambda[i*data->ncol+j] += log(data->len[i]) + log(data->sum[j]) ;
      lambda[i*data->ncol+j] = exp(lambda[i*data->ncol+j]);
    }
  }
}

void calcLambdaRow(DATA *data, PARAM *param, int r) {
  int j;
  for(j=0;j<param->q;j++) {
    lambda[r*data->ncol+j] = data->label[j] == 1 ? mu0[r] + mu1[r] : mu0[r];
    if(NORMALIZE) lambda[r*data->ncol+j] += log(data->len[r]) + log(data->sum[j]) ;
    lambda[r*data->ncol+j] = exp(lambda[r*data->ncol+j]);
  }
}

void calcLambdaCol(DATA *data, PARAM *param, int c) {
  int i;
  for(i=0;i<param->p;i++) {
    lambda[i*data->ncol+c] = data->label[c] == 1 ? mu0[i] + mu1[i] : mu0[i];
    if(NORMALIZE) lambda[i*data->ncol+c] += log(data->len[i]) + log(data->sum[c]) ;
    lambda[i*data->ncol+c] = exp(lambda[i*data->ncol+c]);
  }
}

void calcLambdaElement(DATA *data, PARAM *param, int r, int c) {
  lambda[r*data->ncol+c] = data->label[c] == 1 ? mu0[r] + mu1[r] : mu0[r];
  if(NORMALIZE) lambda[r*data->ncol+c] += log(data->len[r]) + log(data->sum[c]);
  lambda[r*data->ncol+c] = exp(lambda[r*data->ncol+c]);
}

void runBurnIn(DATA *data, PARAM *param, int *burn) {
  int i;
  fprintf(stderr, "Burn-in: ");
  for(i=0;i<*burn;i++) {
    if((i!=0) & (((i+1) % (_PRINT_FREQ_)) == 0)) {
      fprintf(stderr, "%d\t", i+1);
      // printDP(param, data);
    }
    update_mu0(param, data);
    update_mu1(param, data);
    //update_mu_hyper(param, data);

  }
  fprintf(stderr, "....done.\n");
}

void runGibbs(DATA *data, PARAM *param, int *iter) {
  int i, k, total, count; 
  fprintf(stderr, "Iteration: ");
  total = *iter / _SKIP_;
  count = 0;
  for(i=0;i<*iter;i++) {
    if((i!=0) & (((i+1) % (_PRINT_FREQ_)) == 0)) {
      fprintf(stderr, "%d\t", i+1);
    }
    
    update_mu0(param, data);
    update_mu1(param, data);
    // update_mu_hyper(param, data);

    if((i+1) % _SKIP_ == 0  &&  count < total) {
      /********* Likelihood Evaluation for All Models  *********/        
      for(k=0;k<param->p;k++) {
        param->diff[k][count] = (mu1[k]); // - mu0[k]);
        param->mu0sum[k] += mu0[k];
        param->mu1sum[k] += mu1[k];
      }
      count++;
    }
  }
  for(k=0;k<param->p;k++) {
    param->mu0sum[k] /= ((double) count);
    param->mu1sum[k] /= ((double) count);
  }
  fprintf(stderr, "....done.\n"); 
}

/********************************** REWRITE THIS FUNCTION ***********************************/

void output(FILE *fp, DATA *data, PARAM *param, int *p, int *q, int *nA, int *nB, int *total) {
  /* Declaration */
  int i,j;
  double diff[*p];
  double sd[*p];
  double tmp;
   
  /* Print Out */
  fprintf(fp, "Protein\tLen\t");
  for(j=0;j<*q;j++) fprintf(fp, "S%d\t", data->label[j]);
  fprintf(fp, "LogFoldChange\tZstatistic\n");;

  for(i=0;i<*p;i++) {
    diff[i] = vec_mean(param->diff[i], *total);
    sd[i] = sqrt(vec_var(param->diff[i], *total));
    param->zstat[i] = diff[i] / sd[i];
  }

  if(0) {
    tmp = vec_med(param->zstat, *p);
    for(i=0;i<*p;i++) {
      param->zstat[i] -= tmp;
    }
  }


  for(i=0;i<*p;i++) {
    fprintf(fp, "%s\t%d\t", data->name[i], (int) (data->len[i]));
    for(j=0;j<*q;j++) fprintf(fp, "%d\t", (int) data->d[i*data->ncol+j]); 
    fprintf(fp, "%.3f\t%.4f\n", diff[i], param->zstat[i]);
  }
}

