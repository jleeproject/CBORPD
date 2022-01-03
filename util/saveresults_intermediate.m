% function saveresults(str)
try
    str = A_RUN_PURPOSE;
    setVersion;
    str_filename = sprintf('results_%s_v%s_%s_intermediate.mat',datestr(now,'yymmdd'), strrep(str_version,'.','_'), str);
    eval(sprintf('save %s;',str_filename));
    fprintf('saved %s\n',str_filename);
catch err
    showErrors(err)
end