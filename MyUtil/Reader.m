classdef Reader
    
    properties
        fileID
    end
    
    methods
        function this = Reader(filepath)
            this.fileID = fopen(filepath,'r');
        end
        
%         function str = readLine(this)
%             str=fscanf(this.fileID, '%s');
%         end
        
        function cell_str = read(this)
            cell_str=textscan(this.fileID, '%s');
        end
        
        function close(this)
            fclose(this.fileID);
        end
        
    end
    
end

