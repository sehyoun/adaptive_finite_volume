%% Tutorial 3: Ornstein-Uhlenbeck Process with local adaptivity
%
% One can see a marked up version of the tutorial at <https://sehyoun.com/adaptive_finite_volume>
%
% License:
%    Copyright 2017-2019 SeHyoun Ahn
%    BSD 2-clause see <https://github.com/sehyoun/>
%

clear;
format long;

% Add location of the toolbox
addpath('../src');

%% Set parameters
flow_rate = 1;
death_rates = linspace(0.2, 2, 5);
n_dim = 2;
n_points = 40;
int_sig = 0.1;    % Standard Deviation of the Brownina Motion
make_plots = true;

%% Set Parameters of the Ornstein-Ihlenbeck Process
mu = 0.495.*ones(1, n_dim);
theta = 1.*ones(n_dim,1);
sigma = int_sig.^2.*ones(n_dim, 1);

% %% Define the Exact Steady-State Distribution
g_true = @(x) 1/sqrt(pi^n_dim.*int_sig.^(2*n_dim)).*exp(-sum((x-mu).*(x-mu)./int_sig.^2, 2));

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


t_start = tic;

% Extract interaction structure
grid.extract_edges();


for iter_rates = 1:length(death_rates)
  x_i = grid.edge_midpoints();
  %% Set Drift and Diffusion
  death_rate = death_rates(iter_rates);
  grid.drift = zeros(grid.num_e,1);
  grid.diffusion = zeros(grid.num_e,1);
  for iter_dim = 1:n_dim
    cur_ind = find(grid.e2dir(1:grid.num_e) == iter_dim);
    grid.drift(cur_ind) = -theta(iter_dim).*(x_i(cur_ind, iter_dim) - mu(iter_dim));
    grid.diffusion(cur_ind) = sigma(iter_dim);
  end

  %% Prepare the Transition Matrix
  A_FP = grid.compute_transition_matrix_modified();

  left_boundary = grid.n2bd(:,1) == 0 & (grid.n2bd(:,2) < 0.5);
  ind_left = find(left_boundary);
  n_ind = length(ind_left);
  flow = flow_rate*ones(n_ind, 1);
  direction = ones(n_ind, 1);
  is_left = true(n_ind, 1);

  A_FP = A_FP + grid.compute_transition_matrix_boundary(1:grid.num_n, ones(grid.num_n, 1), -death_rate*ones(grid.num_n, 1), ones(grid.num_n, 1));

  g = A_FP\(grid.compute_transition_matrix_boundary(ind_left, direction, flow, is_left)*ones(grid.num_n, 1));

  % [g, ~] = eigs(A_FP, 1, 'sm');
  g = g./sum(g);

  %% Plot Approximations
  x_i = grid.node_midpoints();
  if make_plots
    figure(1);
    clf
    scatter3(x_i(:,1), x_i(:,2), g./grid.node_weights, '.');
    hold on;
    title(sprintf('Ornstein-Uhlenbeck with Boundary Flow\n Death rate: %1.2f', death_rate));
    xlabel('x_1');
    ylabel('x_2');
    saveas(gcf, sprintf('afv_tutorial4_fig_%d.png',iter_rates));
  end
end
