/* Find Neighbors
 *
 * Args:
 *   n2bd(ind, :): boundary points of node to find neighbors of
 *   n2bd : boundary points of the afv_grid
 *   direction : direction to consider neighbors
 *   is_left : whether looking for left neighbors
 *   pos_candidate : indices of nodes under consideration
 *
 * Output: Indices of neighboring nodes
 *
 * Note:
 *   For numerical safety, this program uses small_num buffer in making comparisons, e.g., <, >, >=, and etc. This small_num has to be consistent with other parts of the programs. For example, it needs to be check with the small_num used in <edge_weights.m>
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
  mwSize n_dim, num_candidate, direction, is_left;
  mwSize ind_candidate, ind_out, ind_dim, ind_start;
  mwSize num_n2bd;
  double *pos_candidate;
  double *val, *val_candidate, *prout, *newptr;
  bool flag;
  const double small_num = 1e-13;

  n_dim = (mwSize) mxGetN(prhs[0])/2;
  val = mxGetPr(prhs[0]);

  val_candidate = mxGetPr(prhs[1]);
  num_n2bd = (mwSize) mxGetM(prhs[1]);

  direction = (mwSize) mxGetScalar(prhs[2]) - 1;

  is_left = (mwSize) mxGetScalar(prhs[3]);

  num_candidate = mxGetM(prhs[4]);
  pos_candidate = mxGetPr(prhs[4]);

  plhs[0] = mxCreateDoubleMatrix(num_candidate, 1, mxREAL);
  prout = mxGetPr(plhs[0]);

  ind_out = 0;
  for (ind_candidate = 0; ind_candidate<num_candidate; ++ind_candidate) {

    ind_start = (mwSize) pos_candidate[ind_candidate] - 1;
    flag = true;

    for (ind_dim = 0; ind_dim<n_dim; ++ind_dim) {
      if (ind_dim == direction) {
        if (is_left == 1) {
          if ((val[ind_dim] - val_candidate[ind_start + num_n2bd*(ind_dim + n_dim)] > small_num) \
              || (val[ind_dim] - val_candidate[ind_start + num_n2bd*(ind_dim + n_dim)] < -small_num)) {
            flag = false;
            break;
          }
        } else {
          if ((val[ind_dim+n_dim] - val_candidate[ind_start + num_n2bd*ind_dim] > small_num) \
              || (val[ind_dim+n_dim] - val_candidate[ind_start + num_n2bd*ind_dim] < -small_num)) {
            flag = false;
            break;
          }
        }
      } else {
        if ((val[ind_dim] >= val_candidate[ind_start + num_n2bd*(ind_dim + n_dim)] - small_num) \
            || (val[ind_dim + n_dim] <= val_candidate[ind_start + num_n2bd*ind_dim] + small_num)) {
          flag = false;
          break;
        }
      }
    }

    if (flag) {
      prout[ind_out] = pos_candidate[ind_candidate];
      ++ind_out;
    }
  }

  newptr = mxRealloc(prout, ind_out*sizeof(*prout));
  mxSetPr(plhs[0], newptr);
  mxSetM(plhs[0], ind_out);
}
