#include "pep.h"

int find_index(double mu, double sigmasq, double a) {
  int res;
  double tmp = (a-mu) / sqrt(sigmasq);
  if(tmp > 5.0) res = 10000;
  else if(tmp < -5.0) res = 0;
  else res = 5000 + ((int) (tmp * 1000.0));
  return res;
}

double log_trunc_gaussian(double x, double mu, double sigma2, double min) {
  double res;
  if(x >= min) {
    res = - .5 * (x-mu) * (x-mu) / sigma2 - .5 * log(sigma2 * 2.0 * M_PI);
  }
  else {
    res = gsl_cdf_gaussian_P(min-mu, sqrt(sigma2));
    if(res < _tiny_) res = _tiny_;
    res = log(res);
  }
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
      lambda[i*data->ncol+j] = data->label[j] == 1 ? mu0[i][j/2] + mu1[i] : mu0[i][j/2];
    }
  }
}

void calcLambdaRow(DATA *data, PARAM *param, int r) {
  int j;
  for(j=0;j<param->q;j++) {
    lambda[r*data->ncol+j] = data->label[j] == 1 ? mu0[r][j/2] + mu1[r] : mu0[r][j/2];
  }
}

void calcLambdaCol(DATA *data, PARAM *param, int c) {
  int i;
  for(i=0;i<param->p;i++) {
    lambda[i*data->ncol+c] = data->label[c] == 1 ? mu0[i][c/2] + mu1[i] : mu0[i][c/2];
  }
}

void calcLambdaElement(DATA *data, PARAM *param, int r, int c) {
  lambda[r*data->ncol+c] = data->label[c] == 1 ? mu0[r][c/2] + mu1[r] : mu0[r][c/2];
}

void runBurnIn(DATA *data, PARAM *param, int *burn) {
  int i;
  fprintf(stderr, "Burn-in: ");
  for(i=0;i<*burn;i++) {
    if((i!=0) & (((i+1) % (_PRINT_FREQ_)) == 0)) {
      fprintf(stderr, "%d\t", i+1);
    }
    update_mu0(param, data);
    update_mu1(param, data);
    update_sigmasq0(param, data);
  }

  fprintf(stderr, "....done.\n");
}

void runGibbs(DATA *data, PARAM *param, int *iter) {
  int i, k, l, total, count; 
  fprintf(stderr, "Iteration: ");
  total = *iter / _SKIP_;
  count = 0;
  for(i=0;i<*iter;i++) {
    if((i!=0) & (((i+1) % (_PRINT_FREQ_)) == 0)) {
      fprintf(stderr, "%d\t", i+1);
    }
    
    update_mu0(param, data);
    update_mu1(param, data);
    update_sigmasq0(param, data);

    if((i+1) % _SKIP_ == 0  &&  count < total) {
      /********* Likelihood Evaluation for All Models  *********/        
      for(k=0;k<param->p;k++) {
        param->diff[k][count] = mu1[k];
        for(l=0;l<data->npair;l++) {
            param->mu0sum[k][l] += mu0[k][l];
        }
        param->mu1sum[k] += mu1[k];
        param->sigmasq0sum[k] += sigmasq0[k];
      }
      count++;
    }
  }
  for(k=0;k<param->p;k++) {
    for(l=0;l<data->npair;l++) param->mu0sum[k][l] /= ((double) count);
    param->mu1sum[k] /= ((double) count);
    param->sigmasq0sum[k] /= ((double) count);
  }
  fprintf(stderr, "....done.\n"); 
}

/********************************** REWRITE THIS FUNCTION ***********************************/

void output(FILE *fp, DATA *data, PARAM *param, int *p, int *q, int *nA, int *nB, int *total) {
  /* Declaration */
  int i,j;
  double diff[*p];
  double sd[*p];
   
  /* Print Out */
  fprintf(fp, "Protein\t");
  for(j=0;j<*q;j++) fprintf(fp, "S%d\t", data->label[j]);
  fprintf(fp, "LogFoldChange\tZstatistic\n"); ;

  for(i=0;i<*p;i++) {
    diff[i] = vec_mean(param->diff[i], *total);
    sd[i] = sqrt(vec_var(param->diff[i], *total));
    param->zstat[i] = diff[i] / sd[i];
  }

  for(i=0;i<*p;i++) {
    fprintf(fp, "%s\t", data->name[i]);
    for(j=0;j<*q;j++) fprintf(fp, "%.4f\t", data->dat[i*data->ncol+j]); 
    fprintf(fp, "%.3f\t%.4f\n", diff[i], param->zstat[i]);
  }
}



