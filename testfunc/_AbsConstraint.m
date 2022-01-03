classdef AbsConstraint < handle
    % AbsObjectiveFunction
    %   Detailed explanation goes here
    properties (Constant)
    end
    
    properties 
        isMaximize
    end
    
    properties (Abstract)
        fnEval
        xDomain
        dim
        name
        type
    end
    
    methods (Abstract)
        init(this, power)
    end
    
    methods
        function this = AbsObjectiveFunction(isMaximize, power)
            this.isMaximize = isMaximize;
            if(nargin>1)
                this.power = power;
                init(this, power);
            else
                init(this, 1);
            end
%             disp('Objective Function Super Class');
        end
        
        function out = evalWithVecX(this, xx)
            xargs = mat2cell(xx,1,[ones(size(xx))]); % based on min est. mu   
            out = this.fnEval(xargs{:}); 
        end
%         function outputArg = method1(this,inputArg)
%             outputArg = this.Property1 + inputArg;
%         end

            %--------------------- [fnEval] ---------------------
            function output = getFnEval(this)
                output = this.fnEval;
            end


            %--------------------- [xDomain] ---------------------
            function output = getXDomain(this)
                output = this.xDomain;
            end


            %--------------------- [dim] ---------------------
            function output = getDim(this)
                output = this.dim;
            end


            %--------------------- [name] ---------------------
            function output = getName(this)
                output = this.name;
            end

            %--------------------- [optSol] ---------------------
            function output = getOptSol(this)
                output = this.optSol;
            end


            %--------------------- [optValFn] ---------------------
            function output = getOptVal(this)
                output = this.optVal;
            end


    end
end

