classdef AbsFunction < handle & matlab.mixin.Copyable
    % AbsObjectiveFunction
    %   Detailed explanation goes here
    properties (Constant)
    end
    
    properties 
        isMaximize
        modFunc = [];
    end
    
    properties (Abstract)
        fnEval
        xDomain
        dim
        name
        optSol
        optVal
        power;
        type
    end
    
    methods (Abstract)
        init(this, power)
    end
    
    methods
        function this = AbsFunction(isMaximize, power, modFunc)
            this.isMaximize = isMaximize;
            if(nargin<2 || numel(power)==0)
                power = 1;
            end
            if(nargin>1)
                this.power = power;
                init(this, power);
            else
                init(this, 1);
            end
            if(nargin>2)
                if(numel(modFunc)>0)
                    postprocess(this, modFunc)
                    this.modFunc = modFunc;
                    if(numel(modFunc.newDomain)>0)
                        if(sum(size(this.xDomain) == size(modFunc.newDomain))==2)
                            this.xDomain = modFunc.newDomain;
                        else
                            error('New x domain does not match with the existing function dimension ');
                        end
                    end
                    if(numel(modFunc.newOptX)>0)
                        if( numel(modFunc.newOptX) == this.dim )
                            this.optSol = modFunc.newOptX;
                        else
                            error('New optimal x does not match with the existing dimension. ');
                        end
                    end
                end
                
                
            end
%             disp('Objective Function Super Class');
        end
        function postprocess(this, modFunc)
            modFunc.setBaseFnOutput(this.fnEval);
            this.fnEval = @modFunc.adjustedValue;
        end
        
        function out = evalWithVecX(this, xx)
            xargs = mat2arg(xx);
%             xargs = mat2cell(xx,1,[ones(size(xx))]); % based on min est. mu   
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

