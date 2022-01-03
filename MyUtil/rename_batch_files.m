path_target = 'G:\My Drive\01.PhD_Program\01.Research\07.Bayesian Optimization\SDOE'

% dir(path_target);

% path_upper_dir =  'G:\My Drive\01.PhD_Program\01.Research\07.Bayesian Optimization\Noisy BO\Cropped_all\Transformed'

str_prefix = '[Article;SDOE]';
str_regex = '(.*)_cro.*';

% str_prefix = '';
% str_regex = '(.*) ( PDFDrive.*';
ext_file = 'pdf';
% str_regex = '(.*)_k2op.*';
% dirs = dir(sprintf('%s',path_upper_dir));

% for j= 1:size(dirs,1)
%     if(dirs(j).isdir && ~strcmp(dirs(j).name(1),'.')) 
%         path_target = sprintf('%s\\%s',dirs(j).folder,dirs(j).name);
        filenames = dir(sprintf('%s/*.%s',path_target,ext_file));
        for i=1:size(filenames,1)
            fname = filenames(i).name;
            [stt fin toksExt match toks names splits]=regexp(fname,str_regex);

            if(numel(match)==0)
                fprintf('Pass: Already done: %s\n',fname);
                continue;
            end
            old_name = match{1};

            sel_name = toks{1}{1};
            if(length(sel_name)>length(str_prefix) && strcmp(sel_name(1:length(str_prefix)),str_prefix))
                fprintf('Pass: Already done: %s\n',sel_name);
                continue;
            end
        %     fprintf('%d:%s\n',i,fname);
            if(numel(str_prefix)>0)
                new_name = sprintf('%s %s.%s',str_prefix, sel_name,ext_file);
            else
                new_name = sprintf('%s.%s',sel_name,ext_file);
            end

            old_path = sprintf('%s\\%s',path_target,old_name);
            new_path = sprintf('%s\\%s',path_target,new_name);

            [SUCCESS,MESSAGE,MESSAGEID] = movefile(old_path,new_path,'f');
            if(SUCCESS==1)
                fprintf('%d. Renamed : ''%s'' from ''%s''\n',i, new_name, old_name);
            else
                fprintf('filename: %s.\n',old_path);
                fprintf('%d. %d, %s, %s\n',i,SUCCESS,MESSAGE, MESSAGEID);
            end
        end

%     end
% end
