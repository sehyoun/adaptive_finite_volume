function [A_FP] = compute_transition_matrix_modified(obj, weighter)
% Compute the transition matrix corresponding to the Fokker-Planck equation
%
% .. math::
%
%   \frac{dg}{dt} = -\frac{d}{dx}\left(s(x)\cdot g(x)\right) + \nu \frac{d^2 g}{dx^2}
%
% using the modified upwind scheme (equation 2.3) of :cite:`axelsson1979modified`.
%
% Args:
%   e2n : implicit in class
%   drift : implicit in class
%   diffusion : implicit in class
%   node_weights : implicit in class
%   e_weights : implicit in class
%
% Returns:
%   A_FP (:attr:`num_n <afv_grid.num_n>`:math:`\times`:attr:`num_n <afv_grid.num_n>` sparse matrix of doubles): the transition matrix from the FPK-equation
%
% Note:
%
%   - Many separate parts are needed before this function can be called
%     properly. Instead of setting them manually, try to use the given
%     functions that guarantee the internal structure.
%   - This function is for internal edges, for external edges, use
%     :func:`compute_transition_matrix_boundary <afv_grid.compute_transition_matrix_boundary>`.
%
% References:
%   - :cite:`axelsson1979modified`
%

% Reference:
%   Axelsson, O., & Gustafsson, I. (1979). A modified upwind scheme for
%   convective transport equations and the use of a conjugate gradient method
%   for the solution of non-symmetric systems of equations. IMA Journal of
%   Applied Mathematics (Institute of Mathematics and Its Applications), 23(3),
%   321–337. https://doi.org/10.1093/imamat/23.3.321
%

% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  obj.n_weights = obj.node_weights();
  obj.e_weights = obj.edge_weights();

  n_weight_l = obj.n_weights(obj.e2n(:,1));
  n_weight_r = obj.n_weights(obj.e2n(:,2));

  weight_l = obj.e_weights./n_weight_l;
  weight_r = obj.e_weights./n_weight_r;

  sigma = obj.diffusion./obj.compute_diffusion_distance() + 0*obj.drift;
  drift = obj.drift;
  sigma(sigma>0) = sigma(sigma>0)./(1 + abs(drift(sigma>0))./sigma(sigma>0));

  row = zeros(obj.num_e, 4);
  col = zeros(obj.num_e, 4);
  val = zeros(obj.num_e, 4) + 0*drift(1);


  row(:,1) = obj.e2n(:,1);
  col(:,1) = obj.e2n(:,1);
  val(:,1) = -sigma.*weight_l/2 -drift.*weight_l.*(drift>0);

  row(:,2) = obj.e2n(:,1);
  col(:,2) = obj.e2n(:,2);
  val(:,2) = sigma.*weight_r/2 -drift.*weight_r.*(drift<0);

  row(:,3) = obj.e2n(:,2);
  col(:,3) = obj.e2n(:,1);
  val(:,3) = sigma.*weight_l/2 + drift.*weight_l.*(drift>0);

  row(:,4) = obj.e2n(:,2);
  col(:,4) = obj.e2n(:,2);
  val(:,4) = -sigma.*weight_r/2 + drift.*weight_r.*(drift<0);

  A_FP = sparse(row(:), col(:), val(:), obj.num_n, obj.num_n);
end
