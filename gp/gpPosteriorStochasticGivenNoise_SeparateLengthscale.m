function [funcH, funcCovModifiaiblePost] = ...
  gpPosteriorStochasticGivenNoise_SeparateLengthscale(Xte, Xtr, Ytr, thetas, variance, meanFunc, noiseVar, covFunc)
% Outputs the posterior mean (nTe x 1), standard deviations (nTe x 1) and covariance
% matrix (nTe x nTe) of the test data Xte. In addition, returns a function Handle
% for the GPs.

  numTrData = size(Xtr, 1);

  if numTrData == 0
    Ytr = zeros(0, 1);
  end

%   Ktrtr = covFunc(thetas, variance, Xtr, Xtr) + diag(noiseVar);
%   Ktrtr = gpKernelSqExp_SeparateLengthscale(variance, thetas, Xtr, Xtr) + diag(noiseVar);
  if(numel(noiseVar)==1)
%       noiseVar = noiseVar*ones(size(Xtr,1),1);
        Ktrtr = gpKernelSqExp_SeparateLengthscale(variance, thetas, Xtr, Xtr) + noiseVar*eye(size(Xtr,1));
  elseif(numel(noiseVar) ==size(Xtr,1))
        Ktrtr = gpKernelSqExp_SeparateLengthscale(variance, thetas, Xtr, Xtr) + diag(noiseVar);
  else
      error('Dimension dismatch. need to define. (gpPosteriorStochasticGivenNoise)');
  end
  
  Y_ = Ytr - meanFunc(Xtr);
  L = stableCholeskyDec(Ktrtr);
  alpha = L' \ (L \ Y_);
%             gpHPs.covFunc = @gpKernelSqExp;

  % obtain the function handle
  funcH = @(X) gpPosteriorSimple_SeparateLengthscale(X, Xtr, L, alpha, thetas, variance, meanFunc, covFunc);
  funcCovModifiaiblePost = @(nugget) gpPosteriorCovModifierWrapper(Ktrtr, nugget, Xtr, Y_, @gpPosteriorSimple_SeparateLengthscale, thetas, variance, meanFunc, covFunc);
  % Compute outputs for the test data
%   if ~isempty(Xte)
%     [teMean, teK, teStd] = funcH(Xte);
%   else
%     teMean = []; teK = []; teStd = [];
%   end

end



