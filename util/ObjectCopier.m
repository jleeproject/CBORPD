classdef ObjectCopier  
    methods (Static)  
        function that = copy(this)
            if(numel(this)==0)
                that = this;
            elseif(isa(this,'matlab.mixin.Copyable'))
                that = this.copy;
            elseif(isa(this,'struct'))
                that = struct();
                fnames = fieldnames(this);
                for i=1:numel(fnames)
                    fname = fnames{i};
                    val = getfield(this,fname);
                    that = setfield(that,fname,val);
                end
            elseif(isa(this,'cell'))
                that = cell(size(this));
                for i = 1:size(this,1)
                    for j=1:size(this,2)
                        that{i,j} = ObjectCopier.copy(this{i,j});
                    end
                end
            else
%                 if(isa(this,'class'))
                    str = class(this);
%                 else
%                     str = '';
%                 end
                str_err = sprintf('Undefined type: %s',str)
                error(str_err);
            end
        end
    end
end

