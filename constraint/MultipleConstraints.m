classdef MultipleConstraints < handle
    % AbsObjectiveFunction
    %   Detailed explanation goes here
    properties (Constant)
    end
    
    properties 
    end
    
    properties 
        constraints % array of Constraint object
        nConstraints % number of constraints
    end
    
%     methods (Abstract)
% %         init(this)
%     end
    
    methods
        function this = MultipleConstraints(constraints)
            this.constraints = constraints;
            this.nConstraints = numel(constraints);
            init(this)
        end
        
        function init(this)
%             for i=1:this.nConstraints
%                 TypeConstraint.UnknownMeanConstraint, TypeConstraintFuncion.Shc2, [], 10 
%             end
        end
        

    end
end

