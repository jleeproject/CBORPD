function nlml = normMargLikelihoodStochasticGivenNoise_SeparateLengthscale(variance, thetas,  X, Y, meanFunc, noiseVar, covFunc)
% Computes the normalised Marginal Likelihood for a single GP
  numData = size(X, 1);
%   covFunc = gpKernelSqExp_SeparateLengthscale;
  if(numel(noiseVar)==1)
%       noiseVar = noiseVar*ones(size(Xtr,1),1);
    Ky = gpKernelSqExp_SeparateLengthscale(variance, thetas, X, X) + noiseVar*eye(size(X,1));
  elseif(numel(noiseVar) ==size(X,1))
    Ky = gpKernelSqExp_SeparateLengthscale(variance, thetas, X, X) + diag(noiseVar);
  else
      error('Dimension dismatch. need to define. (gpPosteriorStochasticGivenNoise)');
  end


  Y_ = Y - meanFunc(X);
  L = stableCholeskyDec(Ky);
  alpha = L' \ (L \ Y_);
  nlml = -0.5 * Y_' * alpha - sum( log(diag(L)) ) -  log(2*pi)*numData/2;
end

