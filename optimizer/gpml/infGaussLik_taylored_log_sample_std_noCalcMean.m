function [post nlZ dnlZ] = infGaussLik_taylored_log_sample_std(hyp, mean, cov, lik, x, y, opt)

% Exact inference for a GP with Gaussian likelihood.
%
% Compute a parametrization of the posterior, the negative log marginal
% likelihood and its derivatives w.r.t. the hyperparameters. The function takes
% a specified covariance function (see covFunctions.m) and likelihood function
% (see likFunctions.m), and is designed to be used with gp.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2018-08-22.
%                                      File automatically generated using noweb.
%
% See also INFMETHODS.M, APX.M.

if nargin<7, opt = []; end                          % make sure parameter exists
if iscell(lik), likstr = lik{1}; else likstr = lik; end
if ~ischar(likstr), likstr = func2str(likstr); end
if ~strcmp(likstr,'likGauss_taylored_log_sample_std')               % NOTE: no explicit call to likGauss
  error('Exact inference only possible with Gaussian likelihood');
end



[n, D] = size(x);
[m,dm] = feval(mean{:}, hyp.mean, x);           % evaluate mean vector and deriv
% sn2 = exp(2*hyp.lik); 
% -------------------------------------------------------------
% ------------------- JAESUNG LEE -----------------------------
% -------------------------------------------------------------
% When y = log(sample standard deviation)
% Used Sample sizes for historical samples
samplesize = hyp.samplesize;
% sn2 = (samplesize-1)/8 * gamma( (samplesize-1)/2)^4 / gamma(samplesize/2)^4;
% Vector of variances of each sample
sn2 = (samplesize-1)./8 .* exp(4.*(gammaln( (samplesize-1)./2) - gammaln(samplesize./2)));
% -------------------------------------------------------------
W = ones(n,1)./sn2;            % noise variance of likGauss
K = apx(hyp,cov,x,opt);                        % set up covariance approximation
[ldB2,solveKiW,dW,dhyp,post.L] = K.fun(W); % obtain functionality depending on W
% ldB2 = log(det(B))/2, where B = I + K*diag(W)*I
% solveKiW(r) = inv(K+inv(W))*r
% dhyp > struct that has .cov as element : derivative of covariance function
% K.P(r) = r

alpha = solveKiW(y-m);                          % inv(K+inv(W))*(y-m)
post.alpha = K.P(alpha);                       % return the posterior parameters
post.sW = sqrt(W);                              % sqrt of noise precision vector
if nargout>1                               % do we want the marginal likelihood?
%   nlZ = (y-m)'*alpha/2 + ldB2 + n*log(2*pi*sn2)/2;    % -log marginal likelihood = (y-m)'*inv(K+inv(W))*(y-m)/2 + ...
  nlZ = (y-m)'*alpha/2 + ldB2 + sum(log(2*pi*sn2))/2;    % -log marginal likelihood = (y-m)'*inv(K+inv(W))*(y-m)/2 + ...
  if nargout>2                                         % do we want derivatives?
    dnlZ = dhyp(alpha); 
%     dnlZ.mean = -dm(alpha);
%     dnlZ.lik = -sn2*(alpha'*alpha) - 2*sum(dW)/sn2 + n; %% JAESUNG
%     dnlZ.lik = -sn2*(alpha'*alpha) - 2*sum(dW./sn2) + n; %% JAESUNG
  end
end
