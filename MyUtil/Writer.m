classdef Writer
    %LOGGER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fileID
    end
    
    methods
        function this = Writer(filepath)
            this.fileID = fopen(filepath,'w');
        end
        
        function write(this, str)
%             disp(str);
            if(isnumeric(str))
                if(size(str,1)>1 || size(str,2)>1)
                    n = size(str,1);
                    m = size(str,2);
%                     if(isa(str,'double'))
                    for i = 1 : n
                        for j = 1 : m
                            if(mod(str,1)>0)
                                fprintf(this.fileID, '%.12f\t', str(i,j));
                            else
                                fprintf(this.fileID, '%d\t', str(i,j));
                            end
                        end
                        fprintf(this.fileID, '\n');
                    end
%                         fprintf(this.fileID, '%s\n', num2str(str(i,:)));
%                     elseif(isa(str,'double'))
%                     elseif(isa(str,'double'))
%                     end

                    return;
                end
                str = num2str(str);
            end
            fprintf(this.fileID, str);
        end
        
        function close(this)
            fclose(this.fileID);
        end
        
    end
    
end

