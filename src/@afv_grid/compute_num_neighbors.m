function compute_num_neighbors(obj)
% Compute the number of neighbors
%
% Args:
%   n2n : (implicit) :attr:`n2n <afv_grid.n2n>`
%
% Returns:
%   : [num_nb] implicit in class
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  for iter_node = 1:obj.num_n
    for iter_dim = 1:obj.n_dim
      obj.num_nb(iter_node, iter_dim, 1) = length(obj.n2n{iter_node, iter_dim, 1});
      obj.num_nb(iter_node, iter_dim, 2) = length(obj.n2n{iter_node, iter_dim, 2});
    end
  end
end
