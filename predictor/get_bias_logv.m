function bias = get_bias_logv(samplesize)
    bias = log(2) + psi( (samplesize-1)/2 ) - log(samplesize-1);
end