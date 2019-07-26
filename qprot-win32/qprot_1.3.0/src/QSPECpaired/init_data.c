#include "pep.h"

void init_data(DATA *data, int *p, int *q) {
  int i,j;
  data->nrow = *p;
  data->ncol = *q;
  assert(data->label = (int *) calloc(*q, sizeof(int)));
  assert(data->len = (double *) calloc(*p, sizeof(double)));
  assert(data->sum = (double *) calloc(*q, sizeof(double)));
  assert(data->d = (double *) calloc((*p)*(*q), sizeof(double)));
  assert(data->dat = (double *) calloc((*p)*(*q), sizeof(double)));
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
  int i,j,m; 
  char buf[Max_buf];
  double tmp_sum;
  double tmp0, tmp1;

  init_data(data, p, q); 
  
  fscanf(fpm,"%s",buf);
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

    fscanf(fpm,"%s",buf);
    data->len[i] = atof(buf);

    for(j=0;j<data->ncol;j++) {
      fscanf(fpm,"%s",buf);
      m = i*data->ncol+j;
      data->dat[m] = atof(buf);
      data->d[m] = atof(buf);
    }
  }

  tmp0 = 0.0;
  for(j=0;j<data->ncol;j++) {
    if(data->label[j] == 0) tmp0 += 1.0;
  }
  tmp1 = 0.0;
  for(j=0;j<data->ncol;j++) {
    if(data->label[j] == 1) tmp1 += 1.0;
  }

  for(j=0;j<data->ncol;j++) {
    tmp_sum = 0.0;
    for(i=0;i<data->nrow;i++) {
      tmp_sum += data->d[i*data->ncol+j] / data->len[i];
    }
    data->sum[j] = tmp_sum;
  }
  tmp_sum = vec_mean(data->sum, data->ncol);
  // for(j=0;j<data->ncol;j++) data->sum[j] /= tmp_sum;

  return 0;
}





