#include "pep.h"

/********* ALPHA_prey *********/


void update_sigmasq0(PARAM *param, DATA *data) {
  int j,pid, accept;
  double cur_sigmasq0, tmp_lambda, tmp_sigmasq0, tmp, Delta, mhratio;

  for(pid=0;pid<param->p;pid++) {

    cur_sigmasq0 = sigmasq0[pid];
    Delta = gsl_ran_gaussian(r, 0.1);
    if((cur_sigmasq0 + Delta) <= 0.0) {
      accept = 0;
    }
    else {
      mhratio = 0.0;
      for(j=0;j<param->q;j++) {
        tmp_lambda = lambda[pid*data->ncol+j];
        tmp = data->d[pid * data->ncol + j];
        tmp_sigmasq0 = sigmasq0[pid];
        mhratio += log_trunc_gaussian(tmp, tmp_lambda, tmp_sigmasq0+Delta, data->dmin[2*pid+data->label[j]]);
        mhratio -= log_trunc_gaussian(tmp, tmp_lambda, tmp_sigmasq0, data->dmin[2*pid+data->label[j]]);
      }

      mhratio += - b_sigmasq0 / (cur_sigmasq0 + Delta) - (a_sigmasq0 + 1.0) * log(cur_sigmasq0 + Delta);
      mhratio -= - b_sigmasq0 / (cur_sigmasq0) - (a_sigmasq0 + 1.0) * log(cur_sigmasq0);

      accept = gsl_ran_flat(r, 0.0, 1.0) <= GSL_MIN(1.0, exp(mhratio)) ? 1 : 0 ;

      if(accept) {
        sigmasq0[pid] += Delta;
      }
    }

  }

}

void update_sigmasq_hyper(PARAM *param, DATA *data) {

  int i, accept;
  double mhratio;
  double tmp_a, tmp_b, cur_sigmasq0;

  tmp_a = a_sigmasq0 + gsl_ran_gaussian(r, 0.1);
  tmp_b = b_sigmasq0 + gsl_ran_gaussian(r, 0.1);

  if(tmp_a > 1.0 && tmp_b > 0.0) {

    mhratio = 0.0;
    for(i=0;i<param->p;i++) {
      cur_sigmasq0 = sigmasq0[i];
      mhratio += - tmp_b / (cur_sigmasq0) - (tmp_a + 1.0) * log(cur_sigmasq0) + tmp_a * log(tmp_b) - gsl_sf_lngamma(tmp_a);
      mhratio -= - b_sigmasq0 / (cur_sigmasq0) - (a_sigmasq0 + 1.0) * log(cur_sigmasq0) + a_sigmasq0 * log(b_sigmasq0) - gsl_sf_lngamma(a_sigmasq0);
    }
    accept = gsl_ran_flat(r, 0.0, 1.0) <= GSL_MIN(1.0, exp(mhratio)) ? 1 : 0 ;

    if(accept) {
      a_sigmasq0 = tmp_a;
      b_sigmasq0 = tmp_b;
    }
  }

}






