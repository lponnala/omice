#include "pep.h"

/********* mu0 *********/

void update_mu0(PARAM *param, DATA *data) {
  int j, k, pid, accept;
  double cur_mu0, tmp_lambda, tmp, Delta, mhratio;

  for(pid=0;pid<param->p;pid++) {

     for(k=0;k<data->npair;k++) {
    
        cur_mu0 = mu0[pid][k];
        Delta = gsl_ran_gaussian(r, 0.5); 
        mhratio = 0.0;
        for(j=(2*k);j<(2*k+2);j++) {
              tmp_lambda = lambda[pid*data->ncol+j];
              tmp = data->d[pid * data->ncol + j];
              mhratio += log_poisson(tmp, tmp_lambda*exp(Delta));
              mhratio -= log_poisson(tmp, tmp_lambda);
        }
        mhratio += log_gaussian(mu0[pid][k] + (Delta), m_mu0, v_mu0) 
                 - log_gaussian(mu0[pid][k], m_mu0, v_mu0);  
        accept = gsl_ran_flat(r, 0.0, 1.0) <= GSL_MIN(1.0, exp(mhratio)) ? 1 : 0 ;

        if(accept) {
           mu0[pid][k] += Delta;
           for(j=(2*k);j<(2*k+2);j++) {
              lambda[pid*data->ncol+j] *= exp(Delta);
           }
        } 

     }

  }

}


void update_mu1(PARAM *param, DATA *data) {
  int j,pid, accept;
  double cur_mu1, tmp_lambda, tmp, Delta, mhratio;

  for(pid=0;pid<param->p;pid++) {
    cur_mu1 = mu1[pid];
    Delta = gsl_ran_gaussian(r, 0.5); 
    mhratio = 0.0;
    for(j=0;j<param->q;j++) {
      if(data->label[j] == 1) {
        tmp_lambda = lambda[pid*data->ncol+j];
        tmp = data->d[pid * data->ncol + j];
        mhratio += log_poisson(tmp, tmp_lambda*exp(Delta));
        mhratio -= log_poisson(tmp, tmp_lambda);
      }
    }
    mhratio += log_gaussian(mu1[pid] + (Delta), m_mu1, v_mu1) 
             - log_gaussian(mu1[pid], m_mu1, v_mu1);  
    accept = gsl_ran_flat(r, 0.0, 1.0) <= GSL_MIN(1.0, exp(mhratio)) ? 1 : 0 ;

    if(accept) {
      mu1[pid] += Delta;
      for(j=0;j<param->q;j++) {
        if(data->label[j] == 1) lambda[pid*data->ncol+j] *= exp(Delta);
      }
    }
  }

}




void update_mu_hyper(PARAM *param, DATA *data) {
  int k;
  double tmp_mean;
  double tmp_var;
  
  tmp_mean = 0.0;
  tmp_var = 0.0;
 
  for(k=0;k<data->npair;k++) {
     tmp_mean += vec_mean(mu0[k], param->p);
     tmp_var += vec_var(mu0[k], param->p);
  }

  m_mu0 = tmp_mean / ((double) data->npair);
  v_mu0 = tmp_var / ((double) data->npair) * 4.0;

  m_mu1 = vec_mean(mu1, param->p);
  v_mu1 = vec_var(mu1, param->p) * 4.0;

}





