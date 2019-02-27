function [e_midpoints] = edge_midpoints(obj, ind)
% Compute mid-points of the edges
%
% Args:
%   ind (vector of indices): indices of the edges to compute the mid-points
%
% Returns:
%   e_midpoints (vector of doubles): mid-points of the edges
%
% See also:
%   - :func:`compute_edge_midpoints <afv_grid.compute_edge_midpoints>`

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  if nargin < 2
    e_midpoints = (obj.e2bd(:, 1:obj.n_dim) +  obj.e2bd(:, (obj.n_dim+1):2*obj.n_dim))/2;
  else
    e_midpoints = (obj.e2bd(ind, 1:obj.n_dim) +  obj.e2bd(ind, (obj.n_dim+1):2*obj.n_dim))/2;
  end

end
