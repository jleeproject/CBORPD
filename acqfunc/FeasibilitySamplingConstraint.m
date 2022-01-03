classdef FeasibilitySamplingConstraint < handle
    
    properties
        name = 'Sampling Constraint';
        constraints
        nConstraints 
        pofCalculator
    end
    
    methods
        function this = FeasibilitySamplingConstraint(param_BO)
            this.constraints = param_BO.getConstraints.constraints;
            this.nConstraints = param_BO.nConstraints;
            this.pofCalculator = PoFCalculator(param_BO);
        end
        
        function validate(this)
            if(this.isMaximizeObj)
                throwError('[ERROR] LCB:maximization, LCB is to minimize');
                return;
            end
        end
        
%         function init(this)
% %             numDims = this.objFunc.getDim();
%             if( isa(this.objFunc,'AbsFunction'))
% %                 this.xDomain = optProb.objFunc.getXDomain();
%                 numDims = this.objFunc.getDim();
%             elseif( isa(this.objFunc,'GrapheneModelSolver'))
%                 numDims = size(optProb.objFunc.decisionVariables.domains.all,1);
% %                 this.xDomain = optProb.objFunc.decisionVariables.domains.all;
%             else
%                 error('[AbsInfillOptimizer] Undefined Type.');
%             end
%             this.beta_t = @(iter) 0.2 * numDims * log(2*numDims*iter);
%         end
        
        function [out_feas, pp] = isInFeasible(this, multiConstrPredictors, x_pred, urnd)
            pofInCol = this.pofCalculator.getPof(multiConstrPredictors, x_pred);
            pp = prod(pofInCol,2);
            feas = urnd < pp;
            out_feas = 1 -2.*(feas);
        end
        
        function [pof] = pof(this, multiConstrPredictors, x_pred, varargin)
            pofInCol = this.pofCalculator.getPof(multiConstrPredictors, x_pred);
            pof = prod(pofInCol,2);
        end

        function [out_feas] = isThisValueFeasible(this, obs)
            % obs: [ # samples to eval , # constraints]
            %  - rows: over samples
            %  - cols: over constraints
            
            
            out_feas = ones(size(obs,1),1);
            for i=1:this.nConstraints
                if(this.constraints{i}.hasLb()&&this.constraints{i}.hasUb())
                    lb = this.constraints{i}.lb;
                    ub = this.constraints{i}.ub;
                    feas = lb<obs(:,i) || obs(:,i)<ub;
                elseif(this.constraints{i}.hasLb())
                    lb = this.constraints{i}.lb;
                    feas = lb<obs(:,i);
%                     ps(:,i) = normcdf( 1 - normcdf( (lb-mus(:,i))./sigmas(:,i)) );
                elseif(this.constraints{i}.hasUb())
                    ub = this.constraints{i}.ub;
                    feas = obs(:,i)<ub;
%                     ps(:,i) = normcdf( (ub-mus(:,i))./sigmas(:,i)  ) ;
                else
                    throwError('Both Lb and Ub are not specified.');
                end
                out_feas = out_feas.*feas; % if one of constraints is not satisfied, infeasible.
            end
        end
        
        function [out_feas_gap] = feasGap(this, obs)
            % obs: [ # samples to eval , # constraints]
            %  - rows: over samples
            %  - cols: over constraints
            
            
            out_feas_gap = zeros(size(obs,1),1);
            for i=1:this.nConstraints
                if(this.constraints{i}.hasLb()&&this.constraints{i}.hasUb())
                    lb = this.constraints{i}.lb;
                    ub = this.constraints{i}.ub;
                    feas_gap = max( [ lb-obs(:,i),  obs(:,i)-ub, zeros(size(obs,1),1) ], [], 2 );
                elseif(this.constraints{i}.hasLb())
                    lb = this.constraints{i}.lb;
                    feas_gap = max( [lb-obs(:,i), zeros(size(obs,1),1)], [], 2 );
%                     ps(:,i) = normcdf( 1 - normcdf( (lb-mus(:,i))./sigmas(:,i)) );
                elseif(this.constraints{i}.hasUb())
                    ub = this.constraints{i}.ub;
                    feas_gap = max( [obs(:,i)-ub, zeros(size(obs,1),1)], [], 2 );
%                     ps(:,i) = normcdf( (ub-mus(:,i))./sigmas(:,i)  ) ;
                else
                    throwError('Both Lb and Ub are not specified.');
                end
                out_feas_gap = max( [out_feas_gap  feas_gap], [], 2); % if one of constraints is not satisfied, infeasible.
            end
        end
        
    end
    
%         function [out_feas, pp] = isInFeasible(this, multiConstrPredictors, x_pred, urnd)
%             % g2<UB
%             % g1>LB
%             % OUTPUT : [row vector of UB-g2 or g1-LB]
%             % INPUT: Each Constraint's : function, lb, ub (in param_BO), predictor, 
% 
% %             [mus, sigmas] = multiConstrPredictors(x_pred);
% %             [mu, sigma] = predictor(x_pred);
% %             numDims = size(x_pred, 2);
%             [mus, sigmas] = multiConstrPredictors.predict(x_pred);
%             
%             ps = zeros(size(x_pred,1),this.nConstraints);
%             for i=1:this.nConstraints
%                 if(this.constraints{i}.hasLb()&&this.constraints{i}.hasUb())
%                     lb = this.constraints{i}.lb;
%                     ub = this.constraints{i}.ub;
%                     ps(:,i) = normcdf( (ub-mus(:,i))./sigmas(:,i)) - normcdf( (lb-mus(:,i))./sigmas(:,i));
%                 elseif(this.constraints{i}.hasLb())
%                     lb = this.constraints{i}.lb;
%                     ps(:,i) = normcdf( 1 - normcdf( (lb-mus(:,i))./sigmas(:,i)) );
%                 elseif(this.constraints{i}.hasUb())
%                     ub = this.constraints{i}.ub;
%                     ps(:,i) = normcdf( (ub-mus(:,i))./sigmas(:,i)  ) ;
%                 else
%                     throwError('Both Lb and Ub are not specified.');
%                 end
%             end
%             
%             pp = prod(ps,2);
% %             feas = pp <urnd;
%             feas = urnd < pp;
%             out_feas = 1 -2.*(feas);
% %             out = -1.*(feas) + (1-feas);
% %             if(urnd < pp)
% %                 % Feasible
% %                 out = -1;
% %             else
% %                 % Infeasible
% %                 out = 1;
% % %                 out = 10 * (pp-urnd);
% %             end
%         end
%     end
end
