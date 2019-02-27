function find_neighbor_structure(obj, new_ind, flag_compute_left)
% Construct neighbor structure
%
% Args:
%   new_ind (**sorted** vector of indices): (optional) indices of new nodes to append to the neighbor structure
%   flag_compute_left (bool): (default: false) whether to compute left neighbors as well
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

  if nargin < 3
    flag_compute_left = false;
  end

  if nargin > 1    % only update with new nodes
    new_ind = sort(new_ind(:));
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
  else    % update from scratch
    obj.n2n = cell(obj.num_n, obj.n_dim, 2);
    start_time = tic;
    for iter_node = 1:obj.num_n
      for iter_dim = 1:obj.n_dim
        obj.n2n{iter_node, iter_dim, 2} = obj.find_neighbor(iter_node, iter_dim, false);

        if flag_compute_left
          obj.n2n{iter_node, iter_dim, 1} = obj.find_neighbor(iter_node, iter_dim, true);
        end
      end
      if ~mod(iter_node - 1, floor(obj.num_n/100))
        fprintf('%3.2f %% done (estimated time to finish %10.1f seconds) \n', (iter_node/obj.num_n)*100, (obj.num_n-iter_node)./iter_node*toc(start_time));
      end
    end
  end
end
