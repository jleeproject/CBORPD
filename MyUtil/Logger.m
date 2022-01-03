classdef Logger < handle
    %LOGGER manage logs;
    % log( string )
    % info( string )
    % warn( string )
    % error( string )
    
    properties
        fileID
        writeOnFile = true;
        tics = zeros(10,1,'uint64');
        nTics = 10;
        namesTics =cell(10,1);
        isDebugging = false;
        filename = '';
    end
    
    methods
        function this = Logger(param_filepath)

            if(exist('param_filepath','var'))
                [dirExist, parantFolder, isDir ] = isDirExisting(param_filepath);
                if(isDir)
                   this.error('Log file name is specified as a Directory name.') 
                end
                if(~dirExist )
                    mkdir(parantFolder);                
                end
                this.filename = param_filepath;
                this.fileID = fopen(param_filepath,'w');
                this.writeOnFile = true;
                this.log(sprintf('Save Log on File:\t %s\n',param_filepath));
            else
                this.fileID = nan;
                this.writeOnFile = false;
            end
        end
        
        function log(this, str)
%             disp(str);
            if(isnumeric(str))
                if(size(str,1)>0 || size(str,2)>0)
                    n = size(str,1);
                    c=fix(clock);
                    if(this.writeOnFile)
                        for i = 1 : n
                            fprintf(this.fileID, '[Log %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),num2str(str(i,:)));
                            fprintf('[Log %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),num2str(str(i,:)));
                        end
                    else
                        for i = 1 : n
                            fprintf('[Log %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),num2str(str(i,:)));
                        end
                    end
                    return;
                end
                str = num2str(str);
            end
%             disp(this.fileID);
%             fileID = fopen(this.filepath,'w');
            c=fix(clock);
            if(this.writeOnFile)
                fprintf('[Log %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
                fprintf(this.fileID, '[Log %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
            else
                fprintf('[Log %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
            end
%             fclose(fileID);
        end

        function info(this, str)
            if(isnumeric(str))
                str = num2str(str);
            end
            
%             disp(['WARNING:' str] );
            c=fix(clock);
%             fprintf(this.fileID, 'WARNING:%2d:%2d:%2d:%s\n', c(1),c(2),c(3), str);
            if(this.writeOnFile)
                fprintf( '[Info %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
                fprintf(this.fileID, '[Info %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
            else
                fprintf( '[Info %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
            end

        end
        
        function warn(this, str)
            if(isnumeric(str))
                str = num2str(str);
            end
            
%             disp(['WARNING:' str] );
            c=fix(clock);
%             fprintf(this.fileID, 'WARNING:%2d:%2d:%2d:%s\n', c(1),c(2),c(3), str);
            if(this.writeOnFile)
                fprintf( '[WARNING %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
                fprintf(this.fileID, '[WARNING %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
            else
                fprintf( '[WARNING %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
            end

        end
        function error(this, str)
            if(isnumeric(str))
                str = num2str(str);
            end
%             disp(['ERROR:' str] );
            c=fix(clock);
%             fprintf(this.fileID, '[ERROR %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
%             fprintf(this.fileID, 'ERROR:%2d:%2d:%2d:%s\n', c(1),c(2),c(3), str);
            if(this.writeOnFile)
                fprintf('[ERROR %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
                fprintf(this.fileID, '[ERROR %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
            else
                fprintf('[ERROR %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
            end

        end
        
        function close(this)
            fclose(this.fileID);
        end
        
        function writeError(this, errObj)
            error(this, errObj.message);
            nErrorLines = size(errObj.stack,1);
            for i = 1:nErrorLines
                error(this, sprintf(' > Stack:%2d',i) );
                error(this, errObj.stack(1).file);
                error(this, errObj.stack(i).name);
                error(this, errObj.stack(i).line);
            end
%             errObj.stack.file
%             errObj.stack.name
%             errObj.stack.line
        end

        function seperate(this, str)
            if(this.writeOnFile)
                if(exist('str','var'))
                    fprintf(this.fileID, '=============== %s =====================================\n',str);
                    fprintf('=============== %s =====================================\n',str);
                else
                    fprintf(this.fileID, '====================================================\n');
                    fprintf('====================================================\n');
                end
            else
                if(exist('str','var'))
                    fprintf('=============== %s =====================================\n',str);
                else
                    fprintf('====================================================\n');
                end
            end
        end

        function newline(this)
            if(this.writeOnFile)
                fprintf(this.fileID, '\n');
            else
                fprintf('\n');
            end
        end

        function print(this, str)
            if(this.writeOnFile)
                fprintf(this.fileID, '%s\n',str);
                fprintf('%s\n',str);
            else
                fprintf('%s\n',str);
            end
        end
        
        function show(this)
            fprintf('logger.info(  sprintf(''%%'',)  )\n');
        end
        
        function startTrack(this, num, str)
            if(~exist('num','var')); num = 10;end;
            if(ischar(num)); str = num; num = 10; end;
            if(exist('str','var')); 
                this.log(sprintf('[ START ] %s',str)); 
                this.namesTics{num,1} = str;
            else
                this.namesTics{num,1} = '';
            end;
            if(num<=this.nTics)
                this.tics(num,1) = tic();
            else
                this.error('Out of tic index.')
            end
        end
        
        function elapsedTime = endTrack(this, num, str)
            if(~exist('num','var')); num = 10;end;
            if(ischar(num)); str = num; num = 10; end;
            if(~exist('str','var')); 
%                 str = '';
                str = sprintf(': %s ',this.namesTics{num,1});
            else;
                str= sprintf(': %s ',str);
            end;
            if(num<=this.nTics)
                elapsedTime = toc(this.tics(num,1));
%                 elpasedTime = elapsedT;
%                 if(elapsedT<=60)
%                     this.log(sprintf('[ Elapsed Time %s] %2.1f sec in total',str,elapsedT ))
%                 elseif(elapsedT <=3600)
%                     totElaspedT = elapsedT;
%                     min = floor(elapsedT / 60);
%                     sec = (elapsedT - 60*min);
%                     this.log(sprintf('[ Elapsed Time %s] %2d min, %.1f sec (%.1fSec in total)',str,min, sec, totElaspedT ))
%                 elseif(elapsedT <=3600*24)
%                     totElaspedT = elapsedT;
%                     hr = floor(elapsedT / 3600);
%                     elapsedT = elapsedT - hr*3600;
%                     min = floor(elapsedT / 60);
%                     sec = (elapsedT - 60*min);
%                     this.log(sprintf('[ Elapsed Time %s] %2d hr, %2d min, %.1f sec (%.1fSec in total)',str,hr, min, sec, totElaspedT ))
%                 else
%                     totElaspedT = elapsedT;
%                     day = floor(elapsedT / 3600/24);
%                     elapsedT = elapsedT - day*3600*24;
%                     hr = floor(elapsedT / 3600);
%                     elapsedT = elapsedT - hr*3600;
%                     min = floor(elapsedT / 60);
%                     sec = (elapsedT - 60*min);
%                     this.log(sprintf('[ Elapsed Time %s] %2d day, %2d hr, %2d min, %.1f sec ( %.1fSec in total)',str,day, hr, min, sec, totElaspedT )   )
%                 end
                this.log( sprintf('[ Elapsed Time %s] %s',str, showPrettyElapsedTime(elapsedTime)) );
            else
                this.error('Out of tic index.')
            end
        end
        
        function debug(this, str)
            if(this.isDebugging)
                if(isnumeric(str))
                    str = num2str(str);
                end

    %             disp(['WARNING:' str] );
                c=fix(clock);
    %             fprintf(this.fileID, 'WARNING:%2d:%2d:%2d:%s\n', c(1),c(2),c(3), str);
                if(this.writeOnFile)
                    fprintf(this.fileID, '[Debug %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
                else
                    fprintf( '[Debug %2d/%2d/%2d %2d:%2d:%2d] %s\n', c(1),c(2),c(3),c(4),c(5),c(6),str);
                end
            end

        end
        
        function enableDebugMode(this)
            this.isDebugging = true;
        end
        
        function disableDebugMode(this)
            this.isDebugging = false;
        end
        
        function printVariable(this, variable, str)
            varName = inputname(2);
            if(nargin==2)
                hasStr = false;
                str = '';
            else
                hasStr = true;
            end
            
            
            if( isa(variable,'logical') )
                if(variable)
                    if(hasStr)
                        this.print(sprintf(' ? O :%s', str));
                    else
                        this.print(sprintf(' ? O :%s', varName));
                    end
                else
                    if(hasStr)
                        this.print(sprintf(' ? X :%s', str));
                    else
                        this.print(sprintf(' ? X :%s ', varName));
                    end
                end
                    
            elseif(  isa(variable,'char')  )
                if(hasStr)
                    this.print(sprintf(' %s : %s', str, variable));
                else
                    this.print(sprintf(' @ %s : %s', varName, variable));
                end
            elseif(  isa(variable,'double')  )
                strValues = '';
                if(max(size(variable))>1)
                    [d1,d2] = size(variable);
                    for i = 1:d1
                        for j=1:d2
                            if(j==1)
                                strValues = strcat(strValues,sprintf('%g', variable(i,j)));
                            else
                                strValues = strcat(strValues,sprintf(', %g', variable(i,j)));
                            end
                        end
                        strValues = strcat(strValues,sprintf(' ; '));
                    end
                else
                    strValues = sprintf('%g', variable);
                end
            
                if(hasStr)
                    this.print(sprintf(' - %s : %s', str, strValues));
                else
                    this.print(sprintf(' - %s : %s', varName, strValues));
                end
            elseif(  isa(variable,'function_handle')  )
                if(hasStr)
                    this.print(sprintf('@ [ Function ] %s : %s', str, func2str(variable)));
                else
                    this.print(sprintf('@ [ Function ] %s : %s', varName, func2str(variable)));
                end
            else
            end
            

        end
        
        function printSingleValue(this, variable, hasStr, str)

        end       
    end
    
end

