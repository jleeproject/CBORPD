function [answer, parantFolder, isDir ]= isDirExisting(filepath)
% If filepath is a folder, then, return true.
% If filepath is a file, then, check the parant folder.

    if(exist(filepath,'dir') == 7)
        answer = true;
        isDir = true;
    else
        isDir = false;
        [parent, currentName,currentExt]=fileparts(filepath);
        parantFolder = parent;
        if(exist(parent,'dir') == 7)
            answer = true;
        else
            answer = false;
        end
    end
end
