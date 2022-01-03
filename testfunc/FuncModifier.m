classdef FuncModifier   < handle 
    properties
        fnOutputModifier % e.g. @(x) x.^2
        fnInputModifier
        hasOutputMod = false;
        hasInputMod = false;
        nModInputDimension = 0;
        baseFnOutput
        %
        newDomain = [];
        newOptX = [];
    end
    
    methods
        function this = FuncModifier(fnOutputModifier, fnInputModifier, newDomain, newOptX)
            this.fnOutputModifier = fnOutputModifier;
            if(numel(fnOutputModifier)>0)
                this.hasOutputMod=true;
            end
            
            if(nargin>1)
                this.fnInputModifier = fnInputModifier;
                if(numel(fnInputModifier)>0)
                    this.hasInputMod=true;
                    this.nModInputDimension = numel(fnInputModifier);
                end
            end
            if(nargin>2)
                this.newDomain = newDomain;
            end
            if(nargin>3)
                this.newOptX = newOptX;
            end
        end
        
        function setBaseFnOutput(this,baseFnOutput)
            this.baseFnOutput = baseFnOutput;
        end
        
        function out = adjustedValue(this, varargin)
            args = varargin;
            if(this.hasInputMod)
                if(this.nModInputDimension> numel(varargin))
                    error('Larger dimension of modification than the original dimension.');
                end
                for i=1:this.nModInputDimension
                    fn = this.fnInputModifier{i};
                    args{i} = fn(args{i});
                end
            end
%             if(numel(this.baseFnOutput)==0)
%                 throwError('this.baseFnOutput is not specified.')
%             end
            if(this.hasOutputMod)
                out = this.fnOutputModifier(this.baseFnOutput(args{:}));
            else
                out = this.baseFnOutput(args{:});
            end
        end
    end
end

