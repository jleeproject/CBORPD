classdef class_template
    %_CLASS_TEMPLATE Summary of this class goes here
    %   Detailed explanation goes here
    
%     properties (Constant)
%         Property1
%     end
%     
%     properties (Abstract)
%         Property1
%     end
    
    properties
        Property1
    end
    
%     methods (Abstract)
%         method2(this)
%     end
    methods
        function obj = class_template(inputArg1,inputArg2)
            % Construct an instance of this class
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        
        function outputArg = method1(obj,inputArg)
            outputArg = obj.Property1 + inputArg;
        end
    end
%     methods (Static)
%         function obj = getObject(type, varargin)
%             switch type
%                 case TypeInfillOptimizer.Direct
%                     obj = InfillOptimizerDirect(varargin{:});
%                 case TypeInfillOptimizer.Grid
%                     obj = InfillOptimizerGrid(varargin{:});
%                 otherwise
%                     throwError('Undefined Class')
%             end
%         end
%     end
%     

end

