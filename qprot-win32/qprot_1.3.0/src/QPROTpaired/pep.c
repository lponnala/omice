#include "pep.h"
int nrow(FILE *fp) {
  char buf[10000];
  int n = 0;
  while(fgets(buf, sizeof(buf), fp) != NULL) n++;
  return n;
}

int newlinechar(char *buf, int k) {
  int i;
  int found = 0;
  for(i=0;i<k;i++) {
    if(buf[i] == '\n') {
      found = 1;
      break;
    }
  }
  return found;
}

int ncol(FILE *fp) {
  char buf[10000];
  int i,cont = 0;
  fgets(buf, sizeof(buf), fp);
  for(i=0;i<10000;i++) {
    if(buf[i] == '\t') cont++;
    if(buf[i] == '\0') break;
  }
  return cont;
}

int main(int argc, char **argv) {

  int p, q, n0, n1, burn, iter, total;

  /* Data */
  FILE *fpm;
  DATA data; 
  /* DATA ratio; */
  PARAM param;
  char outfile[100];
  
  const gsl_rng_type *T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  if (argc != 5) {
    fprintf(stderr, "usage: qprot-paired <matrixData> <nburnin> <niter> <normalize?(0/1)>\n");
    return 1;
  }  

  FILE *fpi = fopen(argv[1], "r");
  p = nrow(fpi)-1;
  rewind(fpi);
  q = ncol(fpi);
  fclose(fpi);

  fpm = fopen(argv[1], "r");
  burn = atoi(argv[2]);
  iter = atoi(argv[3]);
  NORMALIZE = atoi(argv[4]) > 0 ? 1 : 0;

  if(fpm == NULL) { fprintf(stderr, "Matrix data %s does not exist\n", argv[1]); return 1; }
  if(&burn == NULL) { fprintf(stderr, "Number of burnin iterations was not provided\n"); return 1; }
  if(&iter == NULL) { fprintf(stderr, "Number of main interations was not provided\n"); return 1; }
  if(&NORMALIZE == NULL) { fprintf(stderr, "Number of main interations was not provided\n"); return 1; }
 
  /* Read Data */
  read_matrix_data(fpm, &data, &p, &q); 
  fprintf(stderr, "%d Proteins and %d Experiments\n", p, q);
  countFunc(&n0, &n1, data.label, &q);

  /* Estimation */
  total = iter / _SKIP_;
  initialize(&param, &data, &p, &q, &n0, &n1, &total);

  /* Full Model */
  runBurnIn(&data, &param, &burn);
  runGibbs(&data, &param, &iter);

  /* Output */
  strcpy(outfile, argv[1]);
  strcat(outfile, "_qprot");
  FILE *fp_out = fopen(outfile, "w");

  output(fp_out, &data, &param, &p, &q, &n0, &n1, &total);  

  /* print_data(&data); */
  // free_param(&param);
  // free_data(&data);
  /* free_data(&ratio); */
  fclose(fpm);
  fclose(fp_out);
  return 0;
}


