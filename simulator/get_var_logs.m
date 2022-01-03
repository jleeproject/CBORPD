function var_logs = get_var_logs(samplesize)
%     var_logs = (samplesize-1)./8 .* exp(4.*(gammaln( (samplesize-1)./2) - gammaln(samplesize./2)));
    var_logs = 0.5./(samplesize-1);
end