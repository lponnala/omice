#include "pep.h"

void init_data(DATA *data, int *p, int *q) {
  int i,j;
  data->nrow = *p;
  data->ncol = *q;
  assert(data->label = (int *) calloc(*q, sizeof(int)));
  assert(data->d = (double *) calloc((*p)*(*q), sizeof(double)));
  assert(data->dat = (double *) calloc((*p)*(*q), sizeof(double)));
  assert(data->dmin = (double *) calloc(2*(*p), sizeof(double)));
  assert(data->experiment = (char **) calloc(*q, sizeof(char *)));   
  for(j=0;j<*q;j++) assert(data->experiment[j] = (char *) calloc(Max_ID, sizeof(char)));
  assert(data->name = (char **) calloc(*p, sizeof(char *))); 
  for(i=0; i<*p; i++) {
    assert(data->name[i] = (char *) calloc(Max_ID, sizeof(char))); 
  }
}

void free_data(DATA *data) {
  int i;
  free(data->label);
  for(i=0; i<data->nrow; i++) {
    free(data->name[i]); 
  }  
  for(i=0;i<data->ncol;i++) free(data->experiment[i]); 
  free(data->experiment);
  free(data->d);
  free(data->dat);
  free(data->name); 
}

void read_vector_data(FILE *fp, double vec[], int *q) {
  int j;
  char buf[Max_buf];
  for(j=0;j<*q;j++) {
    fscanf(fp,"%s",buf);
    vec[j] = atof(buf);
  }
}

int read_matrix_data(FILE *fpm, DATA *data, int *p, int *q) {
  int i,j; 
  int m;
  int n0, n1;
  int ct0, ct1, ct;
  double gmin, tmp, tmpsum0, tmpsum1;
  char buf[Max_buf];
  double *nonzero;
  double nonzero_mean;
  double nonzero_var;
  assert(nonzero = (double *) calloc((*p)*(*q), sizeof(double)));

  init_data(data, p, q); 
  
  fscanf(fpm,"%s",buf);
  for(j=0;j<data->ncol;j++) {
    fscanf(fpm,"%s",buf);
    data->label[j] = atoi(buf);
    if(data->label[j] != 0 && data->label[j] != 1) {
      fprintf(stderr, "Label should be either zero or one.\n");
      return 1;
    }
  }
  
  for(i=0;i<data->nrow;i++) {
    fscanf(fpm,"%s",buf);
    strcpy(data->name[i], buf);
    for(j=0;j<data->ncol;j++) {
      m = i*data->ncol+j;
      fscanf(fpm,"%s",buf);
      data->dat[m] = atof(buf);
      data->d[m] = atof(buf);
      if(data->dat[m] == 0.0) data->d[m] = GSL_NEGINF;
      else data->d[m] = log(data->d[m]);
    }
  }

  /* centering data */
  m = 0;
  for(i=0;i<(data->nrow * data->ncol);i++) {
    if(data->dat[i] > 0.0) {
      nonzero[m] = data->d[i];
      m++;
    }
  }
  nonzero_mean = vec_mean(nonzero, m);
  nonzero_var = vec_var(nonzero, m);
  // fprintf(stderr, "nonzero_med=%.5f\n", nonzero_mean);
  for(i=0;i<(data->nrow * data->ncol);i++) {
    if(data->dat[i] > 0.0) {
      data->d[i] = (data->d[i] - nonzero_mean); //  / sqrt(nonzero_var);
    }
  }

  if(NORMALIZE) normalizeData(data);

  m = 0;
  for(i=0;i<(data->nrow * data->ncol);i++) {
    if(data->dat[i] > 0.0) {
      nonzero[m] = data->d[i];
      m++;
    }
  }
  gmin = vec_perc(nonzero, 0.025, m);


  ct0 = 0;
  ct1 = 0;
  for(j=0;j<*q;j++) {
     if(data->label[j] == 0) ct0++;
     else ct1++;
  }

  for(i=0;i<*p;i++) {
    n0 = 0; 
    n1 = 0;
    for(j=0;j<*q;j++) {
      if(data->label[j] == 0 && data->dat[i*data->ncol+j] > 0.0) n0++;
      if(data->label[j] == 1 && data->dat[i*data->ncol+j] > 0.0) n1++;
    }

    
    data->dmin[2*i] = GSL_POSINF;
    for(j=0;j<*q;j++) {
      if(data->dat[i*data->ncol+j] > 0.0 && data->d[i*data->ncol+j] < data->dmin[2*i] && data->label[j] == 0) {
        data->dmin[2*i] = data->d[i*data->ncol+j];
      }
    }    

    data->dmin[2*i+1] = GSL_POSINF;
    for(j=0;j<*q;j++) {
      if(data->dat[i*data->ncol+j] > 0.0 && data->d[i*data->ncol+j] < data->dmin[2*i+1] && data->label[j] == 1) {
        data->dmin[2*i+1] = data->d[i*data->ncol+j];
      }
    }    
 
    if(data->dmin[2*i] == GSL_POSINF) {
      if(n1 > 2) data->dmin[2*i] = gmin;  
      else data->dmin[2*i] = data->dmin[2*i+1];  
    }
    else if(data->dmin[2*i+1] == GSL_POSINF) {
      if(n0 > 2) data->dmin[2*i+1] = gmin; 
      else data->dmin[2*i+1] = data->dmin[2*i];  
    }
    else {
      tmp = data->dmin[2*i] > data->dmin[2*i+1] ? data->dmin[2*i+1] : data->dmin[2*i];
      data->dmin[2*i] = tmp;
      data->dmin[2*i+1] = tmp;
    }
    
  }

  tmpsum0 = 0.0;
  tmpsum1 = 0.0;
  for(j=0;j<*q;j++) {
    if(data->label[j] == 0) tmpsum0 += 1.0;
    if(data->label[j] == 1) tmpsum1 += 1.0;
  }

  for(i=0;i<*p;i++) {
    tmp = 0.0;
    ct = 0;
    for(j=0;j<*q;j++) {
      if(data->label[j] == 0 && data->dat[i*data->ncol+j] > 0.0) {
        tmp += data->d[i*data->ncol+j];
        ct++;
      }
    }
    tmp /= tmpsum0;
    /*
    if(ct > 0) data->dmin[2*i] = tmp;
    else data->dmin[2*i] = gmin;
    */

    tmp = 0.0;
    ct = 0;
    for(j=0;j<*q;j++) {
      if(data->label[j] == 1 && data->dat[i*data->ncol+j] > 0.0) {
        tmp += data->d[i*data->ncol+j];
      }
    }
    tmp /= tmpsum1;
    /*
    if(ct > 0) data->dmin[2*i+1] = tmp;
    else data->dmin[2*i+1] = gmin;
    */
  }

  free(nonzero);
  return 0;
}



