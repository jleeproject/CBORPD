classdef CellList < handle
    %CELLLIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        list = cell(0,1);
        length = 0;
    end
    
    methods
        function this = CellList(size)
            this.list = cell(size,1);
            this.length = 0;
        end
        
        function add(this, obj)
            this.length = this.length +1;
            this.list{this.length} = obj;
        end
        
        function res = getList(this)
            res = this.list{1:this.length};
        end
        
        function res = size(this)
            res = this.length;
        end
        
        function res = get(this, index)
            if(index<=this.length)
                res= this.list{index};
                return;
            else
                res = nan;
                return;
            end
        end
        
        function res = set(this, index, obj)
            if(index<=this.length)
                this.list{index} = obj;
                return;
            else
                res = nan;
                return;
            end
        end
        
    end
    
end


