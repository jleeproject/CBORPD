function [yMu, yK, yStd] = gpPosteriorSimple_SeparateLengthscale(Xte, Xtr, L, alpha, thetas, variance, meanFunc, covFunc)

  meanXte = meanFunc(Xte);
%   Ktetr = sqExpKernel(bw, scale, Xte, Xtr);
%   Ktete = sqExpKernel(bw, scale, Xte, Xte);
  Ktetr = gpKernelSqExp_SeparateLengthscale(variance, thetas, Xte, Xtr);
  Ktete = gpKernelSqExp_SeparateLengthscale(variance, thetas, Xte, Xte);

  % Predictive Mean
  yMu = meanXte + Ktetr * alpha;
  % Predictive Variance
  V = L \ (Ktetr)';
  yK = Ktete - V'*V;
  yStd = sqrt(real(diag(yK)));

end