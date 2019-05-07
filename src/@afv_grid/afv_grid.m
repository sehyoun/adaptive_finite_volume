classdef afv_grid < handle & matlab.mixin.Copyable
% adaptively refined grids for fokker-planck equations
%
% Documentation is available at <https://sehyoun.com>
%
% by SeHyoun Ahn, March 2018
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  properties
    n_dim
    n_points

    int_bdries

    n2e
    n2bd
    n2n
    num_nb

    e2n
    e2bd
    e2dir

    num_e
    num_n

    e_weights
    n_weights

    drift
    diffusion

    x_i
    boundaries
  end

  methods
    function obj = afv_grid(n_dim)
      % Empty Constructor
      %
      % Args:
      %   n_dim (integer): Dimensionality of the grid
      %
      % See Also:
      %   :func:`set_dim <afv_grid.set_dim>`
      %

      if nargin < 1
        obj.n_dim = 2;
      else
        obj.n_dim = n_dim;
      end
      obj.set_dim(obj.n_dim);
    end
  end
end
