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

  int subcolsize = (d-1)/zeta+1;
  int subcolsize2 = d - (zeta-1)*subcolsize;
  int idx_per_rand = log_base_d((long) RAND_MAX + 1, subcolsize);
  int bit_per_rand = log_base_d((long) RAND_MAX + 1, 2);

  if (zeta > d) zeta = d;

  double lowval = -1/sqrt((double) zeta);
  double increment = -2*lowval;

  plhs[0] = mxCreateSparse(d,m,nnz,false);
  double *vals  = mxGetPr(plhs[0]);
  mwIndex *rows = mxGetIr(plhs[0]);
  mwIndex *colstarts = mxGetJc(plhs[0]);

  // Set values
  int myrand;
  for (int i = 0; i < nnz; i++) {
    if (i % bit_per_rand == 0) myrand = rand();
    vals[i] = (myrand % 2) * increment + lowval;;
    myrand >>= 1;
  }

  // Set column starts
  for (int i = 1; i < m+1; ++i){
    colstarts[i] = i*zeta;
  }

  // Set row indices
  for (int i = 0; i < m*zeta; i += 1) {
    if (i % idx_per_rand == 0) myrand = rand();
    int j = i % zeta;
    int size = (j == zeta - 1) ? subcolsize : subcolsize2;
    rows[i] = j * subcolsize + myrand % size; 
    myrand /= size;
  }

}
