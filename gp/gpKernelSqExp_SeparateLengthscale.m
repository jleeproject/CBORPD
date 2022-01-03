function K = gpKernelSqExp_SeparateLengthscale(variance, thetas, X, Y)
% Returns the Kernel Matrix for a Gaussian Kernel of bandwidth h.
% X is an nxd matrix. K is an nxn matrix.
% If Y is nonempty then returns the gaussian kernel for XxY

% thetas = 1/2 * legnthscale matrix
  if ~exist('Y', 'var') 
    Y = X;
  end
  try
%     d = X-Y;
%     K = variance .* exp( - d' * diag(thetas) * d);
    D = distSquaredGPModule(X.*sqrt(thetas), Y.*sqrt(thetas));
  catch e
      fprintf('[ERROR] [gpKernelSqExp] %s\n',e);
      showErrors(e);
  end
%   K = variance * exp( -D / (2*lengthscales^2) );
  K = variance * exp( -D  );
end

