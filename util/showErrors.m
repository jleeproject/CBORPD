function showErrors(err)
    estack = err.stack;
    fprintf(2, '[ERROR] %s \n', err.identifier)
    fprintf(2, '%s \n', err.message)
    
    if(numel(err.cause)>0)
        fprintf(2, '\n Caused by : \n')
        for i=1:numel(err.cause)
    %         fprintf(2, 'cause: %s \n', err.cause{i})
            showErrors(err.cause{i});
        end
        fprintf(2, '\n in : \n')
    end
    for i=1:numel(estack)
    %     fprintf('file %s, ', estack(i).file);
        fprintf(2, 'Error in <a href="matlab: opentoline(which(''%s''),%d)">%s</a> (<a href="matlab: opentoline(which(''%s''),%d)">line %d</a>)\n\n', estack(i).file, estack(i).line, estack(i).name,  estack(i).file, estack(i).line, estack(i).line);
%         fprintf(2, 'Error in <a href="matlab: opentoline(which(''%s''),%d)"><strong>%s</strong></a> (<a href="matlab: opentoline(which(''%s''),%d)">line %d</a>)\n\n', estack(i).file, estack(i).line, estack(i).name,  estack(i).file, estack(i).line, estack(i).line);
    %     fprintf('line %s\n', estack(i).line);
    end

end
