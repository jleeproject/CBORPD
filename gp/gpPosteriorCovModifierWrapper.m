%% Adding Nuggets
function funcH = gpPosteriorCovModifierWrapper(Ktrtr, proportion, Xtr, Y_, gpPosterior, varargin)
    nugget = mean(diag(Ktrtr)).*proportion;
  Ktrtr = Ktrtr + nugget.*eye(size(Ktrtr,1));
%   Y_ = Ytr - meanFunc(Xtr);
  L = stableCholeskyDec(Ktrtr);
  alpha = L' \ (L \ Y_);
%             gpHPs.covFunc = @gpKernelSqExp;

  % obtain the function handle
  funcH = @(X) gpPosterior(X, Xtr, L, alpha, varargin{:});
% funcH = @(X) gpPosterior(X, Xtr, L, alpha, bw, scale, meanFunc, covFunc);
end