function out = getFilenameFromFileStruct(fileStruct)
    out = sprintf('%s\\%s',fileStruct.folder,fileStruct.name);
end