function [r_dist] = compute_diffusion_distance(obj, ind)
% INTERNAL METHOD: Calculate the diffusion distance to approximate :math:`\frac{dg}{dx}` for diffusion in finite volume method
%
% Args:
%   ind (vector of indices): indices of the node to compute the diffusion distance
%
% Returns:
%   r_dist (vector of doubles): the distance to be used in diffusion computation.
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  n_midpoints = obj.node_midpoints();
  if nargin == 1
    lower = n_midpoints(obj.e2n(:,1) + (obj.e2dir-1).*obj.num_n);
    upper = n_midpoints(obj.e2n(:,2) + (obj.e2dir-1).*obj.num_n);
  else
    lower = n_midpoints(obj.e2n(ind,1) + (obj.e2dir(ind)-1).*obj.num_n);
    upper = n_midpoints(obj.e2n(ind,2) + (obj.e2dir(ind)-1).*obj.num_n);
  end

  r_dist = upper - lower;
end
