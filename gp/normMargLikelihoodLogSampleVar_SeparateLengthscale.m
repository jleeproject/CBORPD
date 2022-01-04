function nlml = normMargLikelihoodLogSampleVar_SeparateLengthscale(variance, thetas,  X, Y, meanFunc, noiseVar, covFunc)
% Computes the normalised Marginal Likelihood for a single GP
  numData = size(X, 1);
%   covFunc = gpKernelSqExp_SeparateLengthscale;
  Ky = gpKernelSqExp_SeparateLengthscale(variance, thetas, X, X) + diag(noiseVar);
  Y_ = Y - meanFunc(X);
  L = stableCholeskyDec(Ky);
  alpha = L' \ (L \ Y_);
  nlml = -0.5 * Y_' * alpha - sum( log(diag(L)) ) -  log(2*pi)*numData/2;
end

