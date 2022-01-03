reader = Reader('./info/version.txt'); 
try
    cell_str_version = reader.read();
    str_version = getElementFromArrayOrCell( cell_str_version{1},1);
    versionBO = sprintf('BO Version %s',...
        str_version    );
catch e
    e
end
reader.close();