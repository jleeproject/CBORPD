classdef Constraint < handle
    % This class is for formulation
    %
    properties (Constant)
    end
    
    properties 
        lb
        ub
        hasLb
        hasUb
        xDomain
        %
        funcEval
        typeConstraint
    end
    
    methods
        function this = Constraint(typeConstraint, lb, ub)
            this.typeConstraint = typeConstraint;
            
            if(numel(lb)>0)
                this.hasLb = true;
                this.lb = lb;
            else
                this.hasLb = false;
            end
            %
            if(numel(ub)>0)
                this.hasUb = true;
                this.ub = ub;
            else
                this.hasUb = false;
            end
            %
%             this.funcEval = FunctionFactory.getFunction(typeFunc, [], []);
%             if(nargin>4)
%                 this.xDomain = xDomain;
%             else
%                 this.xDomain = this.funcEval.getXDomain;
%             end
            
%             this.type = type;
%             this.nConstraints = numel(constraints);
        end
        
        function setFunc(this, funcEval, xDomain)
            this.funcEval = funcEval;
            if(nargin>2)
                this.xDomain = xDomain;
            else
                this.xDomain = this.funcEval.getXDomain;
            end
        end
        
        function out = getXDomain(this)
            out = this.xDomain;
        end
        
        function out = getFuncEval(this)
            out = this.funcEval;
        end
%         function out = probConstraintsAtX(this, xx)
%             probs = zeros(nConstraints,1);
%             for i=1:nConstraints;
%                 probs(i) = constraints.probConstraintsAtX(this, xx);
%             end
%             out = prod(probs);
% %             xargs = mat2cell(xx,1,[ones(size(xx))]); % based on min est. mu   
% %             out = this.fnEval(xargs{:}); 
%         end
    end
end

