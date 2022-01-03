function K = gpKernelMatern52(bw, scale, X, Y)
% Returns the Kernel Matrix for a Gaussian Kernel of bandwidth h.
% X is an nxd matrix. K is an nxn matrix.
% If Y is nonempty then returns the gaussian kernel for XxY

  if ~exist('Y', 'var') 
    Y = X;
  end
  try
    D = distSquaredGPModule(X, Y);
  catch e
      disp(e);
  end
%   K=sigmaf2 .* ( 1 + sqrt(5).*D./sigmal + 5 .*D.^2./3./sigmal.^2).*exp(-sqrt(5).*D./sigmal)
  K = scale .* ( 1 + sqrt(5).*D./bw + 5 .*D.^2./3./bw.^2).*exp(-sqrt(5).*D./bw);
%   K = scale * exp( -D / (2*bw^2) );
end

