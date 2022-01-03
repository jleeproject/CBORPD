function mkdirIfNotExist(newdir)
    sel_dir = dir(newdir);
    if(numel(sel_dir)==0)
        mkdir(newdir);        
    else
        if(~sel_dir.isdir())
            mkdir(newdir);
        end
    end;
end