function opt_samplesize= getAdjustedSamplesize(orig_std,max_samplesize);
    orig_var = orig_std.^2;
    % orig_var = ones(dim_ss,1)*rng_var';
    if(nargin<2)
        max_samplesize = 200;
    else
    end
    
    samplesizes = [2:max_samplesize]';
    samp_var = get_var_logs(samplesizes);
    fin_var = 1./(1./orig_var + 1./samp_var);

    [~,idxmax] = max( (sqrt(orig_var) - sqrt(fin_var))./samplesizes );

    opt_samplesize = samplesizes(idxmax);
end