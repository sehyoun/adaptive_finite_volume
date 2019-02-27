function [n_midpoints] = node_midpoints(obj, ind)
% Compute mid-points of the nodes
%
% Args:
%   ind (vector of indices): (default: ':') indices of the nodes to compute the mid-points
%
% Returns:
%   : n_midpoints (vector of doubles): mid-points of the nodes
%

% LICENSE:
%    Copyright 2017-2019 SeHyoun Ahn
%    BSD 2-clause see <https://github.com/sehyoun/>
%

  if nargin < 2
    n_midpoints = (obj.n2bd(:, 1:obj.n_dim) + obj.n2bd(:, (obj.n_dim+1):2*obj.n_dim))/2;
  else
    n_midpoints = (obj.n2bd(ind, 1:obj.n_dim) + obj.n2bd(ind, (obj.n_dim+1):2*obj.n_dim))/2;
  end
end
