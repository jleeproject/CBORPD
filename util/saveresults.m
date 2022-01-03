% function saveresults(str)
    str = A_RUN_PURPOSE;
    setVersion;
    str_filename = sprintf('results_%s_v%s_%s.mat',datestr(time_exp_start,'yymmdd'), strrep(str_version,'.','_'), str);
%     fprintf('You want to save file into %s ??\n', str_filename);
%     pause();
    if askYesNo(sprintf('You want to save file into ''%s'' ??',str_filename))
        eval(sprintf('save %s;',str_filename));
        fprintf('saved %s\n',str_filename);
    elseif askYesNo('You want to save into different name??')
        str = input('What name?','s');
        str_filename = sprintf('results_%s_v%s_%s.mat',datestr(now,'yymmdd'), strrep(str_version,'.','_'), str);
        if askYesNo(sprintf('You want to save file into ''%s'' ??',str_filename))
            eval(sprintf('save %s;',str_filename));
            fprintf('saved %s\n',str_filename);
        end
    else
    end
%     save(str_filename,'*','-v7.3');
% end