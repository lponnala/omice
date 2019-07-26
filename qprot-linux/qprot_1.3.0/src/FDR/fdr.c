#include "pep.h"


double log_gaussian_pdf(double x, double mu, double sigmasq) {
    double out;
    out = -.5 * pow(x-mu,2.0) / sigmasq - .5 * log(2.0 * M_PI * sigmasq);
    return out;
}

double gaussian_pdf(double x, double mu, double sigmasq) {
    return exp(log_gaussian_pdf(x, mu, sigmasq));
}


void check_mean_shift(void ) {
	int i, j;
	int ncount;
	int tmp_vec[ntick];
	double tmp;

	for(i=0;i<ntick;i++) tmp_vec[i] = 0;
	for(j=0;j<np;j++) {
		i = 0;
		while(tickMarks[i] < zstat[j]) i++;
		i--;
		(tmp_vec[i])++;		
	}

	// find max
	/* mid = 0;
	cur_max = tmp_vec[mid];
	for(i=1;i<ntick;i++) {
		if(tmp_vec[i] > cur_max) {
			mid = i;
			cur_max = tmp_vec[mid];
		}
	}
	nullmean = (tickMarks[mid] + tickMarks[mid+1]) * 0.5;
	*/

	ncount = 0;
	tmp = 0.0;
	for(i=0;i<np;i++) {
		if(fabs(zstat[i] - nullmean) < 3.0) {
			tmp += zstat[i];	
			ncount++;
		}
	}
	nullmean = tmp / ((double) (ncount));

	ncount = 0;
	tmp = 0.0;
	for(i=0;i<np;i++) {
		if(fabs(zstat[i] - nullmean) < 3.0) {
			tmp += pow(zstat[i] - nullmean, 2.0);			
			ncount++;
		}
	}
	nullvar = tmp / ((double) (ncount-1));

	for(i=0;i<=ntick;i++) {
		f0[i] = gaussian_pdf(tickMarks[i], nullmean, nullvar);
	}
	for(i=0;i<np;i++) {
		df0[i] = gaussian_pdf(zstat[i], nullmean, nullvar);
	}

}


void computeProportion(void ) {
	// get conservative estimates of pi_true: just use the data within a range: (-5,5) for now
	int i;
	double tmpsum0, tmpsum;

	tmpsum0 = 0.0;
	tmpsum = 0.0;
	for(i=0;i<ntick;i++) {
		if(tickMarks[i] > -3.0 && tickMarks[i+1] < 3.0) {
			tmpsum0 += 0.5 * (f0[i] + f0[i+1]) * (tickMarks[i+1] - tickMarks[i]);
			tmpsum += 0.5 * (f[i] + f[i+1]) * (tickMarks[i+1] - tickMarks[i]);
		}
	}

	pi_true = 1 - tmpsum / tmpsum0;
	if(pi_true < 0.01) pi_true = 0.01;
	fprintf(stderr, "conservative estimate of pi_true is %.6f\n", pi_true);

}


void fit_FDR(void ) {

	int i, j;
	double numer, denom;

	fprintf(stderr, "Model Fitting\n");

	check_mean_shift();
        
        // pi_true = vec_mean(Z, np);
    
        denom = 0.0;
        for(j=0;j<np;j++) denom += 1;

        for(i=0;i<=ntick;i++) {
		numer = 0.0;
		for(j=0;j<np;j++) {
			numer += gaussian_pdf( (zstat[j] - tickMarks[i]) / bandwidth, 0.0, 1.0);
		}
		f[i] = numer / ( denom * bandwidth );
		// if(f[i] < f0[i]) f[i] = f0[i];
        }

	evaluateDensities(); // df0 and df only, not df1
	computeProportion();

	for(i=0;i<=ntick;i++) {
		f1[i] = ( f[i] - (1.0-pi_true) * f0[i] ) / pi_true;
		if(f1[i] < 0.0) f1[i] = 0.0;
	}

} 


void computeFDR(void ) {
	int i,j;
	double tmpsum, tmpsum0;

	// Local
	for(j=0;j<np;j++) {
		fdr[j] = (1.0 - pi_true) * df0[j] / df[j];
		fdr[j] = fdr[j] > 1.0 ? 1.0 : fdr[j];
	}

	// Global - Down
	for(j=0;j<np;j++) {
		i = 0;
		tmpsum = 0.0;
		tmpsum0 = 0.0;
		while(tickMarks[i+1] < zstat[j]) {
			tmpsum += 0.5 * (f[i+1] + f[i]) * (tickMarks[i+1] - tickMarks[i]);
			tmpsum0 += 0.5 * (f0[i+1] + f0[i]) * (tickMarks[i+1] - tickMarks[i]);
			i++;
		}
		
		tmpsum += 0.5 * (df[j] - f[i]) * (zstat[j] - tickMarks[i]);
		tmpsum0 += 0.5 * (df0[j] - f0[i]) * (zstat[j] - tickMarks[i]);
	
		FDRdown[j] = (1.0 - pi_true) * tmpsum0 / tmpsum;
		FDRdown[j] = FDRdown[j] > 1.0 ? 1.0 : FDRdown[j];
	}

	// Global - Up
	for(j=0;j<np;j++) {
		i = ntick;
		tmpsum = 0.0;
		tmpsum0 = 0.0;
		while(tickMarks[i-1] > zstat[j]) {
			tmpsum += 0.5 * (f[i] + f[i-1]) * (tickMarks[i] - tickMarks[i-1]);
			tmpsum0 += 0.5 * (f0[i] + f0[i-1]) * (tickMarks[i] - tickMarks[i-1]);
			i--;
		}
		
		tmpsum += 0.5 * (f[i] - df[j]) * (tickMarks[i] - zstat[j]);
		tmpsum0 += 0.5 * (f0[i] - df0[j]) * (tickMarks[i] - zstat[j]);

		FDRup[j] = (1.0 - pi_true) * tmpsum0 / tmpsum;
		FDRup[j] = FDRup[j] > 1.0 ? 1.0 : FDRup[j];
	}


}



/*********** Density evaluation ************/
void evaluateDensities(void ) {
	int i, j;
	for(j=0;j<np;j++) {
		// f0
		df0[j] = gaussian_pdf(zstat[j], nullmean, nullvar);
		// f
		i = 0;
		while(tickMarks[i] < zstat[j]  && i <= ntick) i++;
		i--;
		df[j] = f[i] + (f[i+1] - f[i]) / (tickMarks[i+1] - tickMarks[i]) * (zstat[j] - tickMarks[i]);
	}
}


void print_result(FILE *fp) {
	int i,j;

	for(j=0;j<nq;j++) fprintf(fp, "%s\t", X[j]);
	fprintf(fp, "fdr\tFDRup\tFDRdown\n");

	for(i=0;i<np;i++) {
		for(j=0;j<nq;j++) fprintf(fp, "%s\t", X[(i+1)*(nq) + j]);
		fprintf(fp, "%.6f\t%.6f\t%.6f\n", fdr[i], FDRup[i], FDRdown[i]);
	}
}


// Fix this function later

void print_density(FILE *fp) {
	int i;
	for(i=0;i<=ntick;i++) {
		fprintf(fp, "%.4f\t%.6f\t%.6f\t%.6f\n", tickMarks[i], (1.0-pi_true) * f0[i], pi_true * f1[i], f[i]);
	}
}


