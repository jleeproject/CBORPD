classdef PoFCalculator < handle
    
    properties
        name = 'PoFCalculator';
        constraints
        nConstraints 
    end
    
    methods
        function this = PoFCalculator(param_BO)
            if(param_BO.hasConstraints)
                this.constraints = param_BO.getConstraints.constraints;
                this.nConstraints = param_BO.nConstraints;
            else
                throwError('[PofCalculator] No constraints are set in parameter BO')
            end
        end
                
        function [pof] = getPof(this, multiConstrPredictors, x_pred, delta)
            if(nargin<4)
                delta = 0;
            end
            try
                [mus, sigmas] = multiConstrPredictors.predict(x_pred);
            catch err
                showErrors(err);
            end
            sigmas =real(sigmas);
            ps = zeros(size(x_pred,1),this.nConstraints);
            for i=1:this.nConstraints
                if(this.constraints{i}.hasLb&&this.constraints{i}.hasUb)
                    lb = this.constraints{i}.lb;
                    ub = this.constraints{i}.ub;
                    ps(:,i) = normcdf( (ub-mus(:,i) + delta)./sigmas(:,i)) - normcdf( (lb-mus(:,i) - delta)./sigmas(:,i));
                elseif(this.constraints{i}.hasLb)
                    lb = this.constraints{i}.lb;
                    ps(:,i) = normcdf( 1 - normcdf( (lb-mus(:,i) - delta)./sigmas(:,i)) );
                elseif(this.constraints{i}.hasUb)
                    ub = this.constraints{i}.ub;
                    ps(:,i) = normcdf( (ub-mus(:,i) + delta )./sigmas(:,i)  ) ;
                else
                    throwError('Both Lb and Ub are not specified.');
                end
            end
            
            pof = prod(ps,2)+1e-10;
%             aug_pof = [pof zeros(size(pof,1),2)];
        end
    end
end
