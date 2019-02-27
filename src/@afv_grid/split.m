function split(obj, ind, x_cuts, flag_compute_left)
% Split nodes
%
% Args:
%   ind (vector of indices): indices of the nodes to split
%   x_cuts (cell of vector of double): for each (node, dimension) the cut points to split the node
%   flag_compute_left (bool): (default: false) whether to compute left neighbors as well
%
% Returns:
%   : updates all node structure and edge structure of the grid
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  if nargin < 4
    flag_compute_left = false;
  end
  n_node = length(ind);
  num_cuts = zeros(n_node, obj.n_dim);
  for i = 1:n_node
    for j = 1:obj.n_dim
      num_cuts(i, j) = length(x_cuts{i, j});
    end
  end

  num_new = prod(num_cuts+1, 2);
  n2bd_new = zeros(sum(num_new), 2*obj.n_dim);
  n2n_new = cell(sum(num_new), obj.n_dim, 2);

  x_knots = cell(obj.n_dim, 1);
  curr = 0;
  for iter_node = 1:n_node
    for iter_dim  = 1:obj.n_dim
      x_knots{iter_dim} = [obj.n2bd(ind(iter_node),iter_dim); x_cuts{iter_node, iter_dim}(:); obj.n2bd(ind(iter_node),obj.n_dim + iter_dim)];
    end
    n2bd_new(curr+1:curr+num_new(iter_node), :) = obj.compute_node_bds(x_knots);
    n2n_new(curr+1:curr+num_new(iter_node), :, :) = repmat(obj.n2n(ind(iter_node), :, :), num_new(iter_node), 1, 1);
    curr = curr + num_new(iter_node);
  end

  ind_new = [ind; (obj.num_n+1:obj.num_n+curr-n_node)'];
  obj.n2bd(ind_new, :) = n2bd_new;
  obj.n2n(ind_new, :, :) = n2n_new;
  obj.num_n = obj.num_n + curr - n_node;

  obj.find_neighbor_structure(ind_new, flag_compute_left);
  obj.extract_edges();
end
