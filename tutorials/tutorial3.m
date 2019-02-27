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
n_adapt = 7;
adapt_fraction = 0.05;
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
for iter_adapt = 1:n_adapt
  grid.extract_edges();
  x_i = grid.edge_midpoints();


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
    scatter3(x_i(:,1), x_i(:,2), g./grid.node_weights);
    hold on;
    scatter3(x_i(:,1), x_i(:,2), g_true(x_i),'.');
    legend('FV scheme','Exact Soln')
    title('Comparison of Calculation');
    saveas(gcf, sprintf('afv_tutorial3_fig_%d.png', iter_adapt-1));
  end

  drift_i = -theta'.*(x_i - mu);
  drift = max(abs(drift_i), [], 2);
  metric = g.*drift;

  [~, ind_adapt] = sort(metric, 'descend');
  ind_adapt = ind_adapt(1:floor(adapt_fraction*length(ind_adapt)));

  grid.split(ind_adapt, mat2cell(x_i(ind_adapt,:), ones(length(ind_adapt),1), ones(n_dim, 1)));
end
