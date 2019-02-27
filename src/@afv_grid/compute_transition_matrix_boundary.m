function [A_eFP] = compute_transition_matrix_boundary(obj, ind, dir, flow, is_left_edge)
% Compute the transition matrix corresponding to the boundary of the given node
%
% Args:
%   ind (vector of indices): indices of the node to compute the edge dynamics
%   dir (vector of integers): coordinate directions of the edges/flows to consider
%   flow (vector of doubles): flow/drift rate across the edges
%   is_left_edge (vector of bool): whether the edge is a left edge or not
%
% Returns:
%   A_eFP (:attr:`num_n <afv_grid.num_n>`:math:`\times`:attr:`num_n <afv_grid.num_n>` sparse matrix of doubles): the transition matrix of the FPK-equation for the given flows across the edges.
%
% Note:
%   There is :func:`compute_transition_matrix_modified <afv_grid.compute_transition_matrix_modified>` function that handles the construction of the transition matrix for internal edges facing other cells. This function is supposed to be used for creating the flows for the boundary conditions.
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  num_ind = length(ind);
  weights = obj.n2bd(ind, obj.n_dim+1:2*obj.n_dim) - obj.n2bd(ind, 1:obj.n_dim);
  if ~isempty(dir)
    weights(sub2ind(size(weights), (1:num_ind)', dir)) = 1;
  end

  weights = prod(weights,2).*flow.*(2*is_left_edge-1);

  weights = flow.*(2*is_left_edge-1);

  A_eFP = sparse(ind(:), ind(:), weights(:), obj.num_n, obj.num_n);
end
