% files = dir();
if isempty(strfind(pwd, '/')) ;    
    filesep = '\';
else
    filesep = '/';
end
files = dir('..');
for i=1:numel(files);
    file = files(i);
    if(startsWith(file.name,'.'))
        continue;
    else
        if(file.isdir)
            rmpath(genpath(sprintf('%s%s%s',file.folder, filesep, file.name)));
            fprintf('removed path:%s%s%s\n',file.folder, filesep, file.name)
        end
    end
end
