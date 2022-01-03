function out = getElementFromArrayOrCell(array, index, varargin)
%     if(nargin>2)
%         if(isa(array,'cell'))
%             out = array{index,varargin{:}};
%         elseif(isnumeric(array))
%             out = array(index);
%         else
%             out = array(index);
%         end
% 
%     else
%         if(numel(array)<index)
%             throwError(sprintf('Index exceeds the number of array elements (%d).',index))
%         end
        if(isa(array,'cell'))
            out = array{index,varargin{:}};
        elseif(isnumeric(array))
            out = array(index,varargin{:});
        else
            out = array(index,varargin{:});
        end
%     end
end