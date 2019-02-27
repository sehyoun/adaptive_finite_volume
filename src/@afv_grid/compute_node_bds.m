function output = compute_node_bds(obj, x_knots)
% INTERNAL FUNCTION: computes boundaries of new nodes
%
% Args:
%   x_knots (cell of vector of doubles): the boundary knots that define the tensor cut split of a node
%
% Returns:
%   output (matrix of doubles): boundaries of the nodes
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  eff_n = cellfun(@(x) length(x) - 1, x_knots);

  prev_eff_n = cumprod(eff_n);
  prev_eff_n = circshift(prev_eff_n,1);
  prev_eff_n(1) = 1;

  post_eff_n = cumprod(eff_n,'reverse');
  post_eff_n = circshift(post_eff_n,-1);
  post_eff_n(end) = 1;

  output = zeros(prod(eff_n), 2*obj.n_dim);
  for i = 1:obj.n_dim
    repeat_knots = repmat(x_knots{i}(1:end-1)',prev_eff_n(i),1);
    repeat_more = repmat(repeat_knots(:), post_eff_n(i),1);
    output(:,i) = repeat_more;

    repeat_knots = repmat(x_knots{i}(2:end)',prev_eff_n(i),1);
    repeat_more = repmat(repeat_knots(:), post_eff_n(i),1);
    output(:,obj.n_dim + i) = repeat_more;
  end
end
