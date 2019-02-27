function add_new_nodes(obj, n2bd_new, ind_to_consider, flag_compute_left)
% Add new nodes to the grid
%
% Args:
%   n2bd_new (* x 2n_dim matrix): boundaries of new node points
%   ind_to_consider (vector of indices): (optional) indices of nodes that have interface with new nodes
%   flag_compute_left (bool): (default: false) whether to compute left neighbors as well.
%
% Returns:
%   : [n2n] (implicit in class)
%
%     - **n2n** : (in class) neighbor structure
%
% Note:
%   The computing neighbor structure from scratch is an expensive operation. Therefore, one should use new_ind whenever possible. This behavior is guaranteed if one computes using :func:`split_init <afv_grid.split_init>` or :func:`split <afv_grid.split>`, so one should use that function whenever it is feasible.
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  num_new = size(n2bd_new, 1);
  obj.n2bd(obj.num_n+1:obj.num_n+num_new, :) = n2bd_new;
  new_ind = (obj.num_n+1:obj.num_n+num_new)';
  if nargin < 4
    flag_compute_left = false;
  end

  if nargin > 2    % Update with subset of nodes
    % Update existing points
    num_to_consider = length(ind_to_consider);
    ind_to_consider = sort(ind_to_consider);
    for iter_consider = 1:num_to_consider
      iter_node = ind_to_consider(iter_consider);
      for iter_dim = 1:obj.n_dim
        indices = union_sorted(obj.n2n{iter_node, iter_dim, 2}(:), new_ind);
        obj.n2n{iter_node, iter_dim, 2} = obj.find_neighbor(iter_node, iter_dim, false, indices);

        if flag_compute_left
          indices = union_sorted(obj.n2n{iter_node, iter_dim, 1}(:), new_ind);
          obj.n2n{iter_node, iter_dim, 1} = obj.find_neighbor(iter_node, iter_dim, true, indices);
        end
      end
    end

    % Update new points
    indices = union_sorted(ind_to_consider(:), new_ind);
    for iter_node = obj.num_n+1 : obj.num_n+num_new
      for iter_dim = 1:obj.n_dim
        obj.n2n{iter_node, iter_dim, 2} = obj.find_neighbor(iter_node, iter_dim, false, indices);

        if flag_compute_left
          obj.n2n{iter_node, iter_dim, 1} = obj.find_neighbor(iter_node, iter_dim, true, indices);
        end
      end
    end
  else
    % Update existing points
    for iter_node = 1:obj.num_n
      for iter_dim = 1:obj.n_dim
        indices = union_sorted(obj.n2n{iter_node, iter_dim, 2}(:), new_ind);
        obj.n2n{iter_node, iter_dim, 2} = obj.find_neighbor(iter_node, iter_dim, false, indices);

        if flag_compute_left
          indices = union_sorted(obj.n2n{iter_node, iter_dim, 1}(:), new_ind);
          obj.n2n{iter_node, iter_dim, 1} = obj.find_neighbor(iter_node, iter_dim, true, indices);
        end
      end
    end

    % Update new points
    for iter_node = obj.num_n+1 : obj.num_n+num_new
      for iter_dim = 1:obj.n_dim
        obj.n2n{iter_node, iter_dim, 2} = obj.find_neighbor(iter_node, iter_dim, false);

        if flag_compute_left
          obj.n2n{iter_node, iter_dim, 1} = obj.find_neighbor(iter_node, iter_dim, true);
        end
      end
    end
  end
  obj.num_n = obj.num_n + num_new;
end
