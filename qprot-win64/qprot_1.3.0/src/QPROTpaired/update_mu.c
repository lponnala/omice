#include "pep.h"

/********* mu0 *********/

void update_mu0(PARAM *param, DATA *data) {
  int j, k, pid, accept;
  double cur_mu0, tmp_lambda, tmp_sigmasq0, tmp, Delta, mhratio;

  for(pid=0;pid<param->p;pid++) {

     for(k=0;k<data->npair;k++) {

        cur_mu0 = mu0[pid][k];
        Delta = gsl_ran_gaussian(r, 0.5); 
        mhratio = 0.0;
        for(j=(2*k);j<(2*k+2);j++) {
              tmp_lambda = lambda[pid*data->ncol+j];
              tmp = data->d[pid * data->ncol + j];
              tmp_sigmasq0 = sigmasq0[pid];
              mhratio += log_trunc_gaussian(tmp, tmp_lambda+Delta, tmp_sigmasq0, data->dmin[2*pid+data->label[j]]);
              mhratio -= log_trunc_gaussian(tmp, tmp_lambda, tmp_sigmasq0, data->dmin[2*pid+data->label[j]]);
        }
        mhratio += log_gaussian(mu0[pid][k] + (Delta), m_mu0, v_mu0) 
                 - log_gaussian(mu0[pid][k], m_mu0, v_mu0);  
        accept = gsl_ran_flat(r, 0.0, 1.0) <= GSL_MIN(1.0, exp(mhratio)) ? 1 : 0 ;

        if(accept) {
           mu0[pid][k] += Delta;
           for(j=(2*k);j<(2*k+2);j++) {
              lambda[pid*data->ncol+j] += Delta;
           }
        }

     }
  }

}



void update_mu1(PARAM *param, DATA *data) {
  int j,pid, accept;
  double cur_mu1, tmp_lambda, tmp_sigmasq0, tmp, Delta, mhratio;

  for(pid=0;pid<param->p;pid++) {
    cur_mu1 = mu1[pid];
    Delta = gsl_ran_gaussian(r, 0.5); 
    mhratio = 0.0;
    for(j=0;j<param->q;j++) {
      if(data->label[j] == 1) {
        tmp_lambda = lambda[pid*data->ncol+j];
        tmp = data->d[pid * data->ncol + j];
        tmp_sigmasq0 = sigmasq0[pid];
        mhratio += log_trunc_gaussian(tmp, tmp_lambda+Delta, tmp_sigmasq0, data->dmin[2*pid+data->label[j]]);
        mhratio -= log_trunc_gaussian(tmp, tmp_lambda, tmp_sigmasq0, data->dmin[2*pid+data->label[j]]);
      }
    }
    mhratio += log_gaussian(mu1[pid] + (Delta), m_mu1, v_mu1) 
             - log_gaussian(mu1[pid], m_mu1, v_mu1);  
    accept = gsl_ran_flat(r, 0.0, 1.0) <= GSL_MIN(1.0, exp(mhratio)) ? 1 : 0 ;

    if(accept) {
      mu1[pid] += Delta;
      for(j=0;j<param->q;j++) {
        if(data->label[j] == 1) lambda[pid*data->ncol+j] += Delta;
      }
    }
  }

}







