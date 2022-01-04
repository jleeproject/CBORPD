function nlml = normMargLikelihoodStochasticGivenNoise(bw, scale, X, Y, meanFunc, noiseVar, covFunc)
  if(numel(noiseVar)==1)
    Ky = covFunc(bw, scale, X, X) + noiseVar * eye(size(X,1));
  elseif(numel(noiseVar) ==size(X,1))
    Ky = covFunc(bw, scale, X, X) + diag(noiseVar);
  else
      error('Dimension dismatch. need to define. (normMargLikelihoodStochasticGivenNoise)');
  end
  
  numData = size(X, 1);
  Y_ = Y - meanFunc(X);
  L = stableCholeskyDec(Ky);
  alpha = L' \ (L \ Y_);
  nlml = -0.5 * Y_' * alpha - sum( log(diag(L)) ) -  log(2*pi)*numData/2;
end

