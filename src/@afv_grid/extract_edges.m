function extract_edges(obj)
% Extract edges from nodal connections. This function only finds internal edges
%
% Args:
%   n2n : (implicit) :attr:`n2n <afv_grid.n2n>`
%   n2bd : (implicit) :attr:`n2bd <afv_grid.n2bd>`
%
% Returns:
%   : [e2bd, e2n, e2dir, e_weights, num_e] (implicit in class)
%
%     - **e2bd** : (in class) boundaries of the edges
%     - **e2n** : (in class) connectivity of the edges to the nodes
%     - **e2dir** : (in class) normal direction of the edge
%     - **e_weights** : (in class) surface area of the edge
%     - **num_e** : (in class) total number of edges
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  n_node = obj.num_n;
  num_n2bd = size(obj.n2bd, 1);
  n_dim = obj.n_dim;

  % Find edges
  ind_left = cell(n_dim, 1);
  ind_right = cell(n_dim, 1);
  ind_dir = cell(n_dim, 1);
  for iter_dim = 1:n_dim
    aux = cell(n_node,1);
    for iter_node = n_node:-1:1
      aux{iter_node} = repmat(iter_node, length(obj.n2n{iter_node, iter_dim, 2}), 1);
    end
    ind_left{iter_dim} = cell2mat(aux);
    ind_right{iter_dim} = cell2mat(obj.n2n(1:n_node, iter_dim, 2));
    ind_dir{iter_dim} = repmat(iter_dim, length(ind_left{iter_dim}), 1);
  end
  ind_left = cell2mat(ind_left);
  ind_right = cell2mat(ind_right);
  ind_dir = cell2mat(ind_dir);

  % Parse edge Informations
  obj.num_e = length(ind_left);

  obj.e2bd = zeros(obj.num_e, 2*n_dim);
  obj.e2bd(:, 1:obj.n_dim) = max(obj.n2bd(ind_left, 1:obj.n_dim), obj.n2bd(ind_right, 1:obj.n_dim));
  obj.e2bd(:, obj.n_dim+1:2*obj.n_dim) = min(obj.n2bd(ind_left, obj.n_dim+1:2*obj.n_dim), obj.n2bd(ind_right, obj.n_dim+1:2*obj.n_dim));
  obj.e2bd( (1:obj.num_e)' + obj.num_e*(ind_dir - 1) ) = obj.n2bd(ind_right + num_n2bd*(ind_dir-1));
  obj.e2bd( (1:obj.num_e)' + obj.num_e*(ind_dir - 1 + n_dim) ) = obj.n2bd(ind_right + num_n2bd*(ind_dir-1)) + 1;    % Adjust to compute weight, will be unwound later.

  weight = obj.e2bd(:, n_dim+1 : 2*n_dim) - obj.e2bd(:, 1:n_dim);
  obj.e_weights = prod(weight, 2);

  obj.e2n = zeros(obj.num_e, 2);
  obj.e2n(:, 1) = ind_left;
  obj.e2n(:, 2) = ind_right;

  obj.e2dir = ind_dir;

  obj.e2bd( (1:obj.num_e)' + obj.num_e*(ind_dir - 1 + n_dim) ) = obj.n2bd(ind_right + num_n2bd*(ind_dir-1));
end
