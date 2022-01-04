function [yMu, yK, yStd] = gpPosteriorSimple(Xte, Xtr, L, alpha, bw, scale, meanFunc, covFunc)

  meanXte = meanFunc(Xte);
%   Ktetr = sqExpKernel(bw, scale, Xte, Xtr);
%   Ktete = sqExpKernel(bw, scale, Xte, Xte);
  Ktetr = covFunc(bw, scale, Xte, Xtr);
  Ktete = covFunc(bw, scale, Xte, Xte);

  % Predictive Mean
  yMu = meanXte + Ktetr * alpha;
  % Predictive Variance
  V = L \ (Ktetr)';
  yK = Ktete - V'*V;
  yStd = sqrt(real(diag(yK)));

end