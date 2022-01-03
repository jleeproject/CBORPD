if isempty(strfind(pwd, '/'));
    filesep = '\';
else
    filesep = '/';
end

if(numel(dir('../MyUtil'))>0)
    addpath('../MyUtil');
end

omit_dir_list = {'save','temp','backup', 'representativeResults'...
    ,'benchmarks_orig_visual','benchmarks_original'};
files = dir();

for i=1:numel(files);
    file = files(i);
    if(startsWith(file.name,'.'))
        continue;
    else
        if(file.isdir)
            if(~sum(contains(omit_dir_list,file.name)))
                addpath(genpath(sprintf('%s%s%s',file.folder,filesep, file.name)));
                fprintf('added path:%s%s%s\n',file.folder, filesep, file.name)
            end
        end
    end
end
