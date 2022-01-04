function [funcH, funcCovModifiaiblePost] = ...
  gpPosteriorLogSampleVar(Xte, Xtr, Ytr, bw, scale, meanFunc, noiseVar, covFunc)
% Outputs the posterior mean (nTe x 1), standard deviations (nTe x 1) and covariance
% matrix (nTe x nTe) of the test data Xte. In addition, returns a function Handle
% for the GPs.

  numTrData = size(Xtr, 1);

  if numTrData == 0
    Ytr = zeros(0, 1);
  end

  Ktrtr = covFunc(bw, scale, Xtr, Xtr) + diag(noiseVar);
  Y_ = Ytr - meanFunc(Xtr);
  L = stableCholeskyDec(Ktrtr);
  alpha = L' \ (L \ Y_);
%             gpHPs.covFunc = @gpKernelSqExp;

  % obtain the function handle
  funcH = @(X) gpPosteriorSimple(X, Xtr, L, alpha, bw, scale, meanFunc, covFunc);
%   funcCovModifiaiblePost = @(X, covFunc) gpPosteriorSimple(X, Xtr, L, alpha, bw, scale, meanFunc, covFunc);
% gpPosteriorCovModifierWrapper(Ktrtr, Xtr, Y_, gpPosterior, varargin)
  
%% function that adds nugget 
  funcCovModifiaiblePost = @(nugget) gpPosteriorCovModifierWrapper(Ktrtr, nugget, Xtr, Y_, @gpPosteriorSimple, bw, scale, meanFunc, covFunc);
  
  % Compute outputs for the test data
%   if ~isempty(Xte)
%     [teMean, teK, teStd] = funcH(Xte);
%   else
%     teMean = []; teK = []; teStd = [];
%   end

end



