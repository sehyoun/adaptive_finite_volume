% Tutorial 1: Grid structure
%
% One can see a marked up version of the tutorial at <https://sehyoun.com/adaptive_finite_volume>
%
% License:
%    Copyright 2017-2019 SeHyoun Ahn
%    BSD 2-clause see <https://github.com/sehyoun/>
%

addpath('../src');

% Initialize Grid
grid = afv_grid;

% Set dimension to 2
grid.set_dim(2);

% At this stage, the grid is one cell with boundaries [0,1] x [0,1]
% we can see this by seeing

disp(grid.num_n);
disp(grid.n2bd);

x_knots{1} = 0.5;
x_knots{2} = [1/3; 2/3];

grid.split_init(1, x_knots);

disp(grid.num_n);
disp(grid.n2bd);

disp(grid.n2n);

x_knots{1,1} = 0.25;
x_knots{1,2} = 1/6;
grid.split(1, x_knots);

grid.extract_edges();

grid.edge_midpoints()
