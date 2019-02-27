function [midpts] = compute_edge_midpoints(obj, ind, dir, is_left_edge)
% Compute the midpoint of the given edge of a node
%
% Args:
%   ind (n-vector of indices): indices of the **node** to compute the mid-point of the edge
%   dir (n-vector of integers): coordinate direction of the edge to consider
%   is_left_edge (n-vector of bools): whether computing the mid-points of the left edge or not
%
% Returns:
%   midpts (n :math:`\times` :attr:`n_dim <afv_grid.n_dim>` of doubles): midpoints of the edge
%
% Note:
%   This function is to calculate external edges of the node that does not face a different node. For internal nodes, consider using function :func:`edge_midpoints <afv_grid.edge_midpoints>`
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  midpts = (obj.n2bd(ind, obj.n_dim+1:2*obj.n_dim) + obj.n2bd(ind, 1:obj.n_dim))/2;
  midpts(sub2ind(size(midpts), (1:length(ind))', dir)) = obj.n2bd(ind, dir + (1-is_left_edge)*obj.n_dim);
end
