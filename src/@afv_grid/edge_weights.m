function [e_weight] = edge_weights(obj,ind)
% Compute weight of the given edge
%
% Args:
%   ind (vector of indices): indices of the edges to compute the surface area
%
% Returns:
%   e_weight (vector of doubles): surface area of the edges
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  small_num = 20*eps;
  if nargin < 2
    e_weight = (obj.e2bd(:, (obj.n_dim+1):2*obj.n_dim) ...
               - obj.e2bd(:, 1:obj.n_dim));
  else
    e_weight = (obj.e2bd(ind, (obj.n_dim+1):2*obj.n_dim) ...
               - obj.e2bd(ind, 1:obj.n_dim));
  end

  e_weight( abs(e_weight) < small_num ) = 1;
  e_weight = prod(e_weight, 2);
end
