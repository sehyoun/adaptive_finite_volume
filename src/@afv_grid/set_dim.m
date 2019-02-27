function set_dim(obj, n)
% Sets dimensionality of the problem. This function works as the initializer for the class
%
% Args:
%   n (int): dimensionality of the problem
%
% Returns:
%   : (implicit in class)
%
% Note:
%   Due to object-oriented design, the set of codes assumes a certain internal consistency of states. Hence, one should always set dimentionality of the problem through :func:`set_dim <afv_grid.set_dim>`.
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  obj.n_dim = n;
  obj.int_bdries = cell(n,2);

  % set connectivity
  obj.n2bd = kron([0,1], ones(1,obj.n_dim));
  obj.e2bd = kron([0,1], ones(1,obj.n_dim));

  obj.num_e = 1;
  obj.num_n = 1;

end
