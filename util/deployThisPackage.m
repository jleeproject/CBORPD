setVersion
target_dir = 'backup/new';
version = sprintf('Version_%s',strrep(str_version,'.','_'));
% version = 'Version_x_x_x';

target_filename = sprintf('%s/%s.zip', target_dir, version);


omit_dir_list = {'save','temp','backup','test', 'representativeResults'...
    ,'benchmarks', 'benchmarks_orig_visual', 'benchmarks_original', 'figs', 'save_recent', 'save_recent2', 'save_recent3','save2'};
omit_files = {'_class_template.m','_enum_template.m','result_itermediate.mat'};
omit_ext = {'asv', 'mat','png', 'PNG'};

files = dir();

filestozip = cell(numel(files),1);
filestozip{1} = 'D:\Project_Matlab\MyUtil';
filecnt = 1;
for i=1:numel(files);
    file = files(i);
    if(startsWith(file.name,'.'))
        continue;
    else
        if(file.isdir)
            if(~sum(contains(omit_dir_list,file.name)))
                filecnt = filecnt+1;
                filestozip{filecnt} = getFilenameFromFileStruct(file);
            end
        else
            mark = false;
            for j = 1:numel(omit_ext)
                if(endsWith(file.name, omit_ext{j}))
                    mark = true;
%                     continue;
                end
            end
            if(mark);continue;end
            if(~sum(contains(omit_files,file.name)))
                filecnt = filecnt+1;
                filestozip{filecnt} = getFilenameFromFileStruct(file);
            end
        end
    end
end


filestozip = filestozip(1:filecnt);

zip(target_filename, filestozip);
fprintf('Finished: Packaging files to ''%s''\n', target_filename);