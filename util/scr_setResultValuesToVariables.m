fnames =fieldnames(result);
nfields = numel(fnames);
for idx_field=1:nfields
    fname = fnames{idx_field};
    str_eval = sprintf('%s = result.%s;', fname, fname);
    eval(str_eval);
end