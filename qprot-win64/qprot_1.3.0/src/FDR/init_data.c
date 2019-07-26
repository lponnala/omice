#include "pep.h"


int read_matrix_data(FILE *fpm) {
	int i,j,N,ctr; 
	char buf[10000];
	// double tmp;

	fprintf(stderr, "Data has %d rows and %d columns\n", np, nq);
	N = (np + 1) * (nq);
	assert(X = (char **) calloc(N, sizeof(char *)));
	for(i=0;i<N;i++) assert(X[i] = (char *) calloc(100, sizeof(char)));
 
	assert(zstat = (double *) calloc(np, sizeof(double)));
	assert(fdr = (double *) calloc(np, sizeof(double)));
	assert(FDRup = (double *) calloc(np, sizeof(double)));
	assert(FDRdown = (double *) calloc(np, sizeof(double)));
  
	ctr = 0;
	for(j=0;j<nq;j++) {
		fscanf(fpm,"%s", buf);
		strcpy(X[ctr], buf);
		ctr++;
	}

	for(i=0;i<np;i++) {
		for(j=0;j<nq;j++) {
			fscanf(fpm,"%s", buf);
			strcpy(X[ctr], buf);
			ctr++;
		}
	}

	for(i=0;i<np;i++) {
		strcpy(buf, X[(i+1)*(nq)+(nq-1)]);   /* last column as statistic */
		zstat[i] = atof(buf);
	}

	return 0;
}


/* Initialize parameters */

void getBandwidth(void ) {
	double tmp_var = vec_var(zstat, np);
	bandwidth = 1.06 * sqrt(tmp_var) / pow(((double) np), 0.2) ;
}


void init_param(void ) {
	int i; 

   	thres = 0.001;

	double xmin, xmax;

	pi_true = 0.1;
	nullmean = 0.0;
	nullvar = 1.0;

	ntick = 10000;
	assert(tickMarks = (double *) calloc(ntick+1, sizeof(double)));
  
	assert(f = (double *) calloc(ntick+1, sizeof(double)));
	assert(f0 = (double *) calloc(ntick+1, sizeof(double)));
	assert(f1 = (double *) calloc(ntick+1, sizeof(double)));

	assert(df = (double *) calloc(np, sizeof(double)));
	assert(df0 = (double *) calloc(np, sizeof(double)));
	assert(df1 = (double *) calloc(np, sizeof(double)));

	xmin = vec_min(zstat, np) - 0.1;
	xmax = vec_max(zstat, np) + 0.1; 

	for(i=0;i<=ntick;i++) tickMarks[i] = xmin + ((double) i) * (xmax-xmin) / ((double) ntick);

	// for(i=0;i<=ntick;i++) fprintf(stderr, "%.4f\n", tickMarks[i]);

	for(i=0;i<=ntick;i++) {
		f0[i] = gaussian_pdf(tickMarks[i], nullmean, nullvar);
		f1[i] = 0.5 * gaussian_pdf(tickMarks[i], -4.0, 1.0) + 0.5 * gaussian_pdf(tickMarks[i], 4.0, 1.0);
		f[i] = pi_true * f1[i] + (1.0-pi_true) * f0[i];
	}

	for(i=0;i<np;i++) {
		df[i] = 0.0;
		df0[i] = 0.0;
		df1[i] = 0.0;
	}

	getBandwidth();
	
}





