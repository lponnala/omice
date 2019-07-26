#include <Rcpp.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <ctime>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calcSTNdistr(NumericMatrix A, double n_A, double n_B, double n_iter, double xbar0_replace, Function calcSTNdx) {
// NumericVector calcSTNdistr(NumericMatrix A, double n_A, double n_B, double n_iter, double xbar0_replace, Function calcs, RObject best_model) {
	int nA = (int) n_A;
	int nB = (int) n_B;
	int iter = (int) n_iter;
	
	boost::random::mt19937 gen;
	// boost::random::mt19937 gen(std::time(0));
	boost::random::uniform_int_distribution<> distA(0, nA-1);
	
	NumericMatrix xbar_star(iter*A.nrow(),2);
	for (int i=0; i<A.nrow(); i++) {
		for (int k=0;  k<iter; k++) {
			int idx = (i*iter) + k;
			xbar_star(idx,0) = 0;
			xbar_star(idx,1) = 0;
			// A
			for (int j=0; j<nA; ++j) { xbar_star(idx,0) += A(i, distA(gen)); }
			if (xbar_star(idx,0)==0) { xbar_star(idx,0) = xbar0_replace; } else { xbar_star(idx,0) /= nA; }
			// B
			for (int j=0; j<nB; ++j) { xbar_star(idx,1) += A(i, distA(gen)); }
			if (xbar_star(idx,1)==0) { xbar_star(idx,1) = xbar0_replace; } else { xbar_star(idx,1) /= nB; }
		}
	}
	NumericVector stn_dist(xbar_star.nrow());
	stn_dist = calcSTNdx(xbar_star);
	/*
	for (int i=0; i<xbar_star.nrow(); ++i) {
		double sdA = as<double>(calcs(best_model,xbar_star(i,1)));
		double sdB = as<double>(calcs(best_model,xbar_star(i,2)));
		stn_dist[i] = (xbar_star(i,2) - xbar_star(i,1))/(sdA + sdB);
	}
	*/
	return stn_dist;
}


// [[Rcpp::export]]
NumericMatrix calcXbarStar(NumericMatrix A, double n_A, double n_B, double n_iter, double xbar0_replace) {
	int nA = (int) n_A;
	int nB = (int) n_B;
	int iter = (int) n_iter;
	
	boost::random::mt19937 gen;
	// boost::random::mt19937 gen(std::time(0));
	boost::random::uniform_int_distribution<> distA(0, nA-1);
	
	NumericMatrix xbar_star(iter*A.nrow(),2);
	for (int i=0; i<A.nrow(); i++) {
		for (int k=0;  k<iter; k++) {
			int idx = (i*iter) + k;
			xbar_star(idx,0) = 0;
			xbar_star(idx,1) = 0;
			// A
			for (int j=0; j<nA; ++j) { xbar_star(idx,0) += A(i, distA(gen)); }
			if (xbar_star(idx,0)==0) { xbar_star(idx,0) = xbar0_replace; } else { xbar_star(idx,0) /= nA; }
			// B
			for (int j=0; j<nB; ++j) { xbar_star(idx,1) += A(i, distA(gen)); }
			if (xbar_star(idx,1)==0) { xbar_star(idx,1) = xbar0_replace; } else { xbar_star(idx,1) /= nB; }
		}
	}
	return xbar_star;
}


/*
// [[Rcpp::export]]
NumericVector calcXbarStarA(NumericMatrix A, int iter, NumericVector ind, double xbar0_replace) {
	NumericVector xbar_star(iter*A.nrow());
	for (int i = 0; i < A.nrow(); i++) {
		for (int k = 0; k < iter; k++) {
			int idx = (i*iter) + k;
			xbar_star[idx] = 0;
			for (int j=0; j < A.ncol(); j++) {
				xbar_star[idx] += A(i, ind[(k*A.ncol()) + j]-1);
			}
			xbar_star[idx] /= A.ncol();
			if (xbar_star[idx]==0) { xbar_star[idx] = xbar0_replace; }
		}
	}
	return xbar_star;
}
*/

/*
// [[Rcpp::export]]
NumericVector calcSTN(NumericMatrix A, double n_A, double n_B, double n_iter, double xbar0_replace, Function calcs, RObject best_model) {
	int nA = (int) n_A;
	int nB = (int) n_B;
	int iter = (int) n_iter;
	
	// boost::random::mt19937 gen;
	boost::random::mt19937 gen(std::time(0));
	
	boost::random::uniform_int_distribution<> distA(0, nA-1);
	
	NumericVector model_stn_dist(iter*A.nrow());
	
	for (int i=0; i<A.nrow(); i++) {
		double sumA = 0; double xbarA = 0; double sdA = 0;
		double sumB = 0; double xbarB = 0; double sdB = 0;
		for (int k=0;  k<iter; k++) {
			// A
			for (int j=0; j<nA; j++) {
				sumA += A(i, distA(gen));
			}
			xbarA = sumA / nA;
			if (xbarA==0) { xbarA = xbar0_replace; }
			sdA = as<double>(calcs(best_model, xbarA));
			// B
			for (int j=0; j<nB; j++) {
				sumB += A(i, distA(gen));
			}
			xbarB = sumB / nB;
			if (xbarB==0) { xbarB = xbar0_replace; }
			sdB = as<double>(calcs(best_model, xbarB));
			// STN
			model_stn_dist[(i*iter) + k] = (xbarB - xbarA) / (sdA + sdB);
		}
	}
	return model_stn_dist;	
}
*/

/*
calculate p-values of numbers in x, using distribution in f
f need not be sorted
*/
// [[Rcpp::export]]
NumericVector calcPvals(NumericVector x, NumericVector f) {
	NumericVector pVal(x.size());
	for (int i = 0; i < x.size(); i++) {
		double p = (double) std::count_if(f.begin(), f.end(), [thresh = x[i]](double elem) { return (elem > thresh); }) / f.size();
		pVal[i] = 2*std::min(p,1-p);
	}
	return pVal;
}


/*
calculate p-value of x, using distribution in f
f need not be sorted

// [[Rcpp::export]]
double calcPval(double x, NumericVector f) {
	double p = (double) std::count_if(f.begin(), f.end(), [x](double elem) { return (elem > x); }) / f.size();
	return 2*std::min(p,1-p);
}
*/

/*
calculate p-values of numbers in x, using distribution in f
f need not be sorted, it will be sorted inside here

// [[Rcpp::export]]
NumericVector calcPvals(NumericVector x, NumericVector f) {
	std::sort(f.begin(),f.end());
	NumericVector pVal(x.size());
	for (int i = 0; i < x.size(); i++) {
		double p = 0;
		for (int j = 0; j < f.size(); j++) {
			if (f[j] > x[i]) {
				p = (double) (f.size() - j)/f.size();
				break;
			}
		}
		pVal[i] = 2 * std::min(p, 1-p);
	}
	return pVal;
}
*/
