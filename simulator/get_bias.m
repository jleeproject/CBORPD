function bias = get_bias(samplesize)
%     logdiff = @(samplesize) (gammaln(samplesize/2) - gammaln( (samplesize-1)/2));
%     bias =  (  log(sqrt(2./(samplesize-1))) + logdiff(samplesize) + 0.5.*( 1 - (samplesize-1)./(2 .* exp(2.*logdiff(samplesize))) )  );
    bias = -0.5./(samplesize-1);
end