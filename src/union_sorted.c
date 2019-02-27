/* Given two (ascending) sorted, unique vectors merge them into sorted vector of unique values
 *
 * Inputs: v1 = (n1 x 1) sorted, unique vector of integers
 *         v2 = (n2 x 1) sorted, unique vector of integers
 * Output: sorted, unique merge of two vector
 *
 * by SeHyoun Ahn
 *
 * LICENSE:
 *   Copyright 2017-2019 SeHyoun Ahn
 *   BSD 2-clause see <https://github.com/sehyoun/>
 *
 */
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  mwSize n1, n2;
  mwSize ind1, ind2, num_remain, ind_out;
  double *prv1, *prv2, *prout;

  n2 = mxGetM(prhs[1]);
  prv2 = mxGetPr(prhs[1]);
  n1 = mxGetM(prhs[0]);

  num_remain = n1 + n2;

  plhs[0] = mxCreateDoubleMatrix(num_remain, 1, mxREAL);
  prout = mxGetPr(plhs[0]);

  // Allows first v1 to be empty, but requires v2 to be non-empty
  if (n1 == 0) {
    for (ind2=0; ind2<n2; ++ind2) {
      prout[ind2] = prv2[ind2];
    }
    return;
  }

  prv1 = mxGetPr(prhs[0]);
  ind1 = 0;
  ind2 = 0;
  ind_out = 0;
  while (num_remain > 0) {
    if ( ind1 >= n1 || prv1[ind1] > prv2[ind2]) {
      prout[ind_out] = prv2[ind2];
      --num_remain;
      ++ind_out;
      ++ind2;
    } else if ( ind2 >= n2 || prv1[ind1] < prv2[ind2]) {
      prout[ind_out] = prv1[ind1];
      --num_remain;
      ++ind_out;
      ++ind1;
    } else {
      prout[ind_out] = prv1[ind1];
      num_remain -= 2;
      ++ind_out;
      ++ind1;
      ++ind2;
    }
  }

  prout = mxRealloc(prout, ind_out*sizeof(*prout));
  mxSetPr(plhs[0], prout);
  mxSetM(plhs[0], ind_out);
}
