function [nbs] = find_neighbor(obj, ind, dir, is_left, varargin)
% Find the neighboring node for the given (single) node
%
% Args:
%   ind (integer): index of the node to find the neighbors of
%   dir (integer): direction to find the neighbors in
%   is_left (bool): whether to look for left neighbors or not
%   varargin{1} (vector of indices): indices of nodes to consider for neighbors (assumed unique)
%   n2bd : (implicit in class) :attr:`n2bd <afv_grid.n2bd>`
%   n_dim : (implicit in class) :attr:`n_dim <afv_grid.n_dim>`
%
% Returns:
%   : nbs (vector of indices): indices of the neighbors
%
% Note:
%   :meth:`find_neighbor_structure` runs this function for all nodes.
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  if nargin == 5
    nbs = find_neighbor_mex(obj.n2bd(ind, :), obj.n2bd, dir, is_left.*1, varargin{1});
  else
    nbs = find_neighbor_mex(obj.n2bd(ind, :), obj.n2bd, dir, is_left.*1, (1:obj.num_n)');
  end

  %% Without C-acceleration
  % if nargin == 5
  %   nbs = varargin{1}(obj.n2bd(ind, obj.n_dim.*(1-is_left) + dir) == obj.n2bd(varargin{1}, obj.n_dim.*is_left + dir));
  % else
  %   nbs = find(obj.n2bd(ind, obj.n_dim.*(1-is_left) + dir) == obj.n2bd(1:obj.num_n, obj.n_dim.*is_left + dir));
  % end

  % for iter_dim = 1:obj.n_dim
  %   if iter_dim ~= dir
  %     nbs = nbs(max(obj.n2bd(ind, iter_dim), obj.n2bd(nbs, iter_dim)) < min(obj.n2bd(ind, iter_dim + obj.n_dim), obj.n2bd(nbs, iter_dim + obj.n_dim)));
  %   end
  % end
  % nbs = nbs(:);
end
