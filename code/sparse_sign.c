#include <math.h> 
#include "mex.h"

int log_base_d(long a, int d) {
    int output = 0;
    while (a > 0) {
        a /= d;
        output += 1;
    }
    return output-1;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

	mwSize d = (mwSize) mxGetScalar(prhs[0]);
	mwSize m = (mwSize) mxGetScalar(prhs[1]);
	mwSize zeta = (mwSize) mxGetScalar(prhs[2]);
	mwSize nnz = m*zeta;

    int idx_per_rand = log_base_d((long) RAND_MAX + 1, d);
    int bit_per_rand = log_base_d((long) RAND_MAX + 1, 2);

	if (zeta > d) zeta = d;

	double lowval = -1/sqrt((double) zeta);
    double increment = -2*lowval;

    plhs[0] = mxCreateSparse(d,m,nnz,false);
    double *vals  = mxGetPr(plhs[0]);
    mwIndex *rows = mxGetIr(plhs[0]);
    mwIndex *colstarts = mxGetJc(plhs[0]);

	// Set values
    int myrand = rand();
	for (int i = 0; i+bit_per_rand < nnz; i += bit_per_rand) {
        for (int j = i; j < i+bit_per_rand; ++j) {
	        vals[j] = (myrand % 2) * increment + lowval;
            myrand = myrand >> 1;
        }
        myrand = rand();
	}
    for (int i = bit_per_rand*(nnz/bit_per_rand); i < nnz; ++i) {
        vals[i] = (myrand % 2) * increment + lowval;
        myrand = myrand >> 1;
    }

	// Set column starts
	for (int i = 1; i < m+1; ++i){
	    colstarts[i] = i*zeta;
	}

	// Set row indices
    myrand = rand();
    int ir = 0;
    for (int i = 0; i < m*zeta; i += zeta) {
	    int idx = 0;
	    while (idx < zeta) {
		    rows[i+idx] = myrand % d;
            ir++;
            if (ir == idx_per_rand) {
                ir = 0;
                myrand = rand();
            } else {
                myrand /= d;
            }
		    int j = 0;
		    for (; j < idx; ++j) {
		        if (rows[i+idx] == rows[i+j]) break;
		    }
		    idx += (int) (j == idx);
	    }
	}

}