%% Tutorial 2: Ornstein-Uhlenbeck Process using adaptive grid data structure, but with adaptivity
%
% One can see a marked up version of the tutorial at
%    <https://sehyoun.com/adaptive_finite_volume>
%

% License:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 0-clause <https://opensource.org/licenses/0BSD>
%

clear;
format long;

% Add location of the toolbox
addpath('../src');


%% Set parameters
n_dim = 2;
n_points = 20;
int_sig = 0.1;    % Standard Deviation of the Brownina Motion
make_plots = true;

%% Set Parameters of the Ornstein-Ihlenbeck Process
mu = 0.495.*ones(1, n_dim);
theta = 1.*ones(n_dim,1);
sigma = int_sig.^2.*ones(n_dim, 1);

%% Define the Exact Steady-State Distribution
g_true = @(x) 1/sqrt(pi^n_dim.*int_sig.^(2*n_dim)).*exp(-sum((x-mu).*(x-mu)./int_sig.^2, 2));

%% Define Adaptivity metric
drop_metric = @(g,x,drift) g < sqrt(eps).*max(g);

fprintf('Processing with %d grid points (per dimension).\n', n_points);

%% Initialize Grid
grid = afv_grid(n_dim);

%% Split Grid
cut_points = cell(n_dim,1);
cut_points_1d = linspace(0, 1, n_points+1);
cut_points_1d = cut_points_1d(2:end-1)';
for iter_dim = 1:n_dim
  cut_points{iter_dim} = cut_points_1d;
end
grid.split_init(1, cut_points);


figure(1);
clf;
x_i = grid.node_midpoints();
scatter(x_i(:,1), x_i(:,2), '.');
xlabel('x');
ylabel('y');
title('Grid');
saveas(gcf,'afv_tutorial2_fig_1.png');
t_start = tic;

% Extract interaction structure
grid.extract_edges();

x_i = grid.edge_midpoints();
hold on;
scatter(x_i(:,1), x_i(:,2), 'r.');
legend('Node Midpoints', 'Edge Midpoints');
saveas(gcf,'afv_tutorial2_fig_2.png');

%% Set Drift and Diffusion
grid.drift = zeros(grid.num_e,1);
grid.diffusion = zeros(grid.num_e,1);
for iter_dim = 1:n_dim
  cur_ind = find(grid.e2dir(1:grid.num_e) == iter_dim);
  grid.drift(cur_ind) = -theta(iter_dim).*(x_i(cur_ind, iter_dim) - mu(iter_dim));
  grid.diffusion(cur_ind) = sigma(iter_dim);
end

%% Prepare the Transition Matrix
A_FP = grid.compute_transition_matrix_modified();

[g, ~] = eigs(A_FP, 1, 'sm');
g = g./sum(g);

%% Plot Approximations
x_i = grid.node_midpoints();
if make_plots
  figure(1);
  clf
  scatter3(x_i(:,1),x_i(:,2),g./grid.node_weights);
  hold on;
  scatter3(x_i(:,1),x_i(:,2),g_true(x_i),'.');
  legend('FV scheme','Exact Soln')
  title('Comparison of Calculation');
  saveas(gcf, 'afv_tutorial2_fig_3.png');
end

%% Set Drift and Diffusion
int_sig = 0.2;
mu = 0.35.*ones(1, n_dim);
theta = 1.*ones(n_dim,1);
sigma = int_sig.^2.*ones(n_dim, 1);

x_i = grid.edge_midpoints();
for iter_dim = 1:n_dim
  cur_ind = find(grid.e2dir(1:grid.num_e) == iter_dim);
  grid.drift(cur_ind) = -theta(iter_dim).*(x_i(cur_ind, iter_dim) - mu(iter_dim));
  grid.diffusion(cur_ind) = sigma(iter_dim);
end

%% Prepare the Transition Matrix
A_FP = grid.compute_transition_matrix_modified();

[g, ~] = eigs(A_FP, 1, 'sm');
g = g./sum(g);

%% Plot Approximations
x_i = grid.node_midpoints();
g_true = @(x) 1/sqrt(pi^n_dim.*int_sig.^(2*n_dim)).*exp(-sum((x-mu).*(x-mu)./int_sig.^2, 2));
if make_plots
  figure(1);
  clf;
  scatter3(x_i(:,1),x_i(:,2),g./grid.node_weights);
  hold on;
  scatter3(x_i(:,1),x_i(:,2),g_true(x_i),'.');
  legend('FV scheme','Exact Soln')
  title('Comparison of Calculation');
  saveas(gcf, 'afv_tutorial2_fig_4.png');
end
