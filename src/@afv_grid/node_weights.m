function [n_weight] = node_weights(obj,ind)
% Compute the weight of the given node
%
% Args:
%   ind (vector of indices): indices of the nodes to compute the volume
%
% Returns:
%   n_weight (vector of doubles): volume of the nodes
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  if nargin < 2
    ind = 1:obj.num_n;
  end

  n_weight = prod((obj.n2bd(ind, (obj.n_dim+1):2*obj.n_dim) - obj.n2bd(ind, 1:obj.n_dim)), 2);
end
