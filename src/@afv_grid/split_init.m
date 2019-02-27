function split_init(obj, node_ind, x_cuts)
% Split a given node which is a starting point
%
% Args:
%   node_ind (integer): index of the cell to split
%   x_cuts (cell of vectors): points to introduce new edges
%
% Returns:
%   : (implicit in class) internally update states to introduce new cells for the given cell
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  if nargin < 3
    x_cuts = num2cell(obj.node_midpoints(node_ind));
  end

  % Get knots
  x_knots = cell(obj.n_dim, 1);
  for iter_dim  = 1:obj.n_dim
    x_knots{iter_dim} = [obj.n2bd(node_ind, iter_dim); x_cuts{iter_dim}(:); obj.n2bd(node_ind, obj.n_dim + iter_dim)];
  end

  % Size of split
  n_nodes = cellfun(@(x) length(x)+1, x_cuts);
  n_nodes = n_nodes(:);
  n_nodes_tot = prod(n_nodes);

  % Node index
  loc_hypercube = reshape(1:n_nodes_tot, n_nodes');

  obj.n2bd(loc_hypercube(:),:) = obj.compute_node_bds(x_knots);
  obj.n2n = -1*ones(n_nodes_tot, obj.n_dim, 2);

  % Use hypercube shifting to fill in n2n values
  for iter_dim = 1:obj.n_dim
    s = struct();
    s.type = '()';
    for i = obj.n_dim:-1:1
      s.subs{i} = ':';
    end

    for iter_right = 0:1
      s.subs{iter_dim} = 1*(1-iter_right) + iter_right*n_nodes(iter_dim);
      aux_n2n = loc_hypercube;
      aux_n2n = circshift(subsasgn(aux_n2n, s, -1), 2*iter_right -1, iter_dim);
      obj.n2n(:, iter_dim, 2 - iter_right) = aux_n2n(:);
    end
  end

  obj.num_n = obj.num_n + numel(loc_hypercube) - 1;
  obj.n2n = num2cell(obj.n2n);
  for iter_dim = 1:obj.n_dim
    for iter_n = 1:obj.num_n
      if obj.n2n{iter_n,iter_dim,1} == -1
        obj.n2n{iter_n,iter_dim,1} = [];
      end
      if obj.n2n{iter_n,iter_dim,2} == -1
        obj.n2n{iter_n,iter_dim,2} = [];
      end
    end
  end
end
