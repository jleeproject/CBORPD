function [funcH, funcCovModifiaiblePost] = ...
  gpPosteriorStochasticGivenNoise(Xte, Xtr, Ytr, bw, scale, meanFunc, noiseVar, covFunc)
% Outputs the posterior mean (nTe x 1), standard deviations (nTe x 1) and covariance
% matrix (nTe x nTe) of the test data Xte. In addition, returns a function Handle
% for the GPs.

  numTrData = size(Xtr, 1);

  if numTrData == 0
    Ytr = zeros(0, 1);
  end

%   if(numel(noiseVar)==1)
%       noiseVar = noiseVar*ones(size(Xtr,1),1);
%   elseif(numel(noiseVar) ~=size(Xtr,1))
%       error('Dimension dismatch. need to define. (gpPosteriorStochasticGivenNoise)');
%   end
  if(numel(noiseVar)==1)
%       noiseVar = noiseVar*ones(size(Xtr,1),1);
        Ktrtr = covFunc(bw, scale, Xtr, Xtr) + noiseVar*eye(size(Xtr,1));
  elseif(numel(noiseVar) ==size(Xtr,1))
        Ktrtr = covFunc(bw, scale, Xtr, Xtr) + diag(noiseVar);
  else
      error('Dimension dismatch. need to define. (gpPosteriorStochasticGivenNoise)');
  end
      
  Y_ = Ytr - meanFunc(Xtr);
  L = stableCholeskyDec(Ktrtr);
  alpha = L' \ (L \ Y_);
%             gpHPs.covFunc = @gpKernelSqExp;

  % obtain the function handle
  funcH = @(X) gpPosteriorSimple(X, Xtr, L, alpha, bw, scale, meanFunc, covFunc);
  funcCovModifiaiblePost = @(nugget) gpPosteriorCovModifierWrapper(Ktrtr, nugget, Xtr, Y_, @gpPosteriorSimple, bw, scale, meanFunc, covFunc);

  
end



