classdef MultipleConstraintsEstimatorGP < handle
    
    properties
        param_BO
        gpParams
        %
        arr_types
        %
        constraints 
        nConstraints = 0;
        %
        cell_estimators
        multiPredictors
        multiTruePredictors
    end
    
    methods
        function this = MultipleConstraintsEstimatorGP(param_BO, gpParams)
            this.param_BO = param_BO;
            this.gpParams = gpParams;
        end
        
        function this = addAllConstraints(this, mutltiConstraints, typeGp)
            if(nargin<3)
                force = false;
            else
                force = true;
            end
            this.constraints = mutltiConstraints;
            this.nConstraints =  mutltiConstraints.nConstraints;
%             %
%             if(numel(arr_types) == 1 && numel(constraints)>1)
%                 arr_types = repmat(arr_types,this.nConstraints,1);
%             elseif(numel(arr_types) ~= numel(constraints))
%                 throwError(sprintf(...
%                     '[MultipleConstraintsEstimatorGP] the number of dimensions does not match: %s, %s',...
%                     'arr_types', 'constraints'));
%             end
%             this.arr_types = arr_types;
            %
            if(numel(this.gpParams) == 1 && this.nConstraints>1)
                this.gpParams = repmat(this.gpParams,this.nConstraints,1);
            elseif(numel(this.gpParams) ~= this.nConstraints)
                error(sprintf(...
                    '[MultipleConstraintsEstimatorGP] the number of dimensions does not match: %s, %s',...
                    'arr_types', 'constraints'));
            end
            %
            cell_estimators = cell(this.nConstraints,1);
            constraints = mutltiConstraints.constraints;
                for i=1:this.nConstraints
                    
                        sel_constraint = getElementFromArrayOrCell(constraints,i);
                    if(force)
                        sel_type = typeGp;
                    else

                        if (contains (this.gpParams{i}.typeGpFit.char,'_SeparateLengthscale'))

                            if(this.param_BO.typeProblem == TypeProblem.ConventionalBoStochastic)
                                if(this.param_BO.isGpFitDirect)
                                    switch sel_constraint.typeConstraint
                                        case TypeConstraint.UnknownMeanConstraint
                                            sel_type = TypeEstimatorGp.DirectStochastic_SeparateLengthscale;
                                        case TypeConstraint.UnknownObsConstraint
                                            sel_type = TypeEstimatorGp.DirectStochastic_SeparateLengthscale;
        %                                 case TypeConstraint.UnknownStdDevConstraint
        %                                     sel_type = TypeEstimatorGp.DirectLogSampleVarRandomLS;
                                        otherwise
                                            throwError('[MultipleConstraintsEstimatorGP] Unsupported types');
                                    end
    %                             elseif(this.param_BO.isGpFitGpml)
    %                                 switch sel_constraint.typeConstraint
    %                                     case TypeConstraint.UnknownMeanConstraint
    %                                         sel_type = TypeEstimatorGp.GpmlStochastic;
    %                                     case TypeConstraint.UnknownObsConstraint
    %                                         sel_type = TypeEstimatorGp.GpmlStochastic;
    %     %                                 case TypeConstraint.UnknownStdDevConstraint
    %     %                                     sel_type = TypeEstimatorGp.GpmlLogSampleVarRandomLS;
    %                                     otherwise
    %                                         throwError('[MultipleConstraintsEstimatorGP] Unsupported types');
    %                                 end
                                else
                                    error('[MultipleConstraintsEstimatorGP] Undefined.');
                                end
                            elseif(this.param_BO.typeProblem == TypeProblem.ConventionalBoDeterministic)
                                if(this.param_BO.isGpFitDirect)
                                    switch sel_constraint.typeConstraint
                                        case TypeConstraint.UnknownMeanConstraint
                                            sel_type = TypeEstimatorGp.DirectDeterministic_SeparateLengthscale;
                                        case TypeConstraint.UnknownObsConstraint
                                            sel_type = TypeEstimatorGp.DirectDeterministic_SeparateLengthscale;
        %                                 case TypeConstraint.UnknownStdDevConstraint
        %                                     sel_type = TypeEstimatorGp.DirectLogSampleVarRandomLS;
                                        otherwise
                                            throwError('[MultipleConstraintsEstimatorGP] Unsupported types');
                                    end
    %                             elseif(this.param_BO.isGpFitGpml)
    %                                 switch sel_constraint.typeConstraint
    %                                     case TypeConstraint.UnknownMeanConstraint
    %                                         sel_type = TypeEstimatorGp.GpmlDeterministic;
    %                                     case TypeConstraint.UnknownObsConstraint
    %                                         sel_type = TypeEstimatorGp.GpmlDeterministic;
    %     %                                 case TypeConstraint.UnknownStdDevConstraint
    %     %                                     sel_type = TypeEstimatorGp.GpmlLogSampleVarRandomLS;
    %                                     otherwise
    %                                         throwError('[MultipleConstraintsEstimatorGP] Unsupported types');
    %                                 end
                                else
                                    error('[MultipleConstraintsEstimatorGP] Undefined.');
                                end
                            elseif(this.param_BO.typeProblem == TypeProblem.RobustDesignOptimization)
                                if(this.param_BO.isGpFitDirect)
                                    switch sel_constraint.typeConstraint
                                        case TypeConstraint.UnknownMeanConstraint
                                            sel_type = TypeEstimatorGp.DirectGivenNoiseVariance_SeparateLengthscale;
                                        case TypeConstraint.UnknownObsConstraint
                                            sel_type = TypeEstimatorGp.DirectGivenNoiseVariance_SeparateLengthscale;
                                        case TypeConstraint.UnknownStdDevConstraint
                                            sel_type = TypeEstimatorGp.DirectLogSampleVar_SeparateLengthscale;
                                        otherwise
                                            throwError('[MultipleConstraintsEstimatorGP] Unsupported types');
                                    end
    %                             elseif(this.param_BO.isGpFitGpml)
    %                                 switch sel_constraint.typeConstraint
    %                                     case TypeConstraint.UnknownMeanConstraint
    %                                         sel_type = TypeEstimatorGp.GpmlStochasticGivenNoiseVar;
    %                                     case TypeConstraint.UnknownObsConstraint
    %                                         sel_type = TypeEstimatorGp.GpmlStochasticGivenNoiseVar;
    %                                     case TypeConstraint.UnknownStdDevConstraint
    %                                         sel_type = TypeEstimatorGp.GpmlLogSampleVarRandomLS;
    %                                     otherwise
    %                                         throwError('[MultipleConstraintsEstimatorGP] Unsupported types');
    %                                 end
                                else
                                    error('[MultipleConstraintsEstimatorGP] Undefined.');
                                end
                            else
                                throwUndefinedTypeError(this.param_BO.typeProblem);
                            end
                        else

                            if(this.param_BO.typeProblem == TypeProblem.ConventionalBoStochastic)
                                if(this.param_BO.isGpFitDirect)
                                    switch sel_constraint.typeConstraint
                                        case TypeConstraint.UnknownMeanConstraint
                                            sel_type = TypeEstimatorGp.DirectStochastic;
                                        case TypeConstraint.UnknownObsConstraint
                                            sel_type = TypeEstimatorGp.DirectStochastic;
        %                                 case TypeConstraint.UnknownStdDevConstraint
        %                                     sel_type = TypeEstimatorGp.DirectLogSampleVarRandomLS;
                                        otherwise
                                            throwError('[MultipleConstraintsEstimatorGP] Unsupported types');
                                    end
                                elseif(this.param_BO.isGpFitGpml)
                                    switch sel_constraint.typeConstraint
                                        case TypeConstraint.UnknownMeanConstraint
                                            sel_type = TypeEstimatorGp.GpmlStochastic;
                                        case TypeConstraint.UnknownObsConstraint
                                            sel_type = TypeEstimatorGp.GpmlStochastic;
        %                                 case TypeConstraint.UnknownStdDevConstraint
        %                                     sel_type = TypeEstimatorGp.GpmlLogSampleVarRandomLS;
                                        otherwise
                                            throwError('[MultipleConstraintsEstimatorGP] Unsupported types');
                                    end
                                else
                                end
                            elseif(this.param_BO.typeProblem == TypeProblem.ConventionalBoDeterministic)
                                if(this.param_BO.isGpFitDirect)
                                    switch sel_constraint.typeConstraint
                                        case TypeConstraint.UnknownMeanConstraint
                                            sel_type = TypeEstimatorGp.DirectDeterministic;
                                        case TypeConstraint.UnknownObsConstraint
                                            sel_type = TypeEstimatorGp.DirectDeterministic;
        %                                 case TypeConstraint.UnknownStdDevConstraint
        %                                     sel_type = TypeEstimatorGp.DirectLogSampleVarRandomLS;
                                        otherwise
                                            throwError('[MultipleConstraintsEstimatorGP] Unsupported types');
                                    end
                                elseif(this.param_BO.isGpFitGpml)
                                    switch sel_constraint.typeConstraint
                                        case TypeConstraint.UnknownMeanConstraint
                                            sel_type = TypeEstimatorGp.GpmlDeterministic;
                                        case TypeConstraint.UnknownObsConstraint
                                            sel_type = TypeEstimatorGp.GpmlDeterministic;
        %                                 case TypeConstraint.UnknownStdDevConstraint
        %                                     sel_type = TypeEstimatorGp.GpmlLogSampleVarRandomLS;
                                        otherwise
                                            throwError('[MultipleConstraintsEstimatorGP] Unsupported types');
                                    end
                                else
                                end
                            elseif(this.param_BO.typeProblem == TypeProblem.RobustDesignOptimization)
                                if(this.param_BO.isGpFitDirect)
                                    switch sel_constraint.typeConstraint
                                        case TypeConstraint.UnknownMeanConstraint
                                            sel_type = TypeEstimatorGp.DirectStochasticGivenNoiseVar;
                                        case TypeConstraint.UnknownObsConstraint
                                            sel_type = TypeEstimatorGp.DirectStochasticGivenNoiseVar;
                                        case TypeConstraint.UnknownStdDevConstraint
                                            sel_type = TypeEstimatorGp.DirectLogSampleVarRandomLS;
                                        otherwise
                                            throwError('[MultipleConstraintsEstimatorGP] Unsupported types');
                                    end
                                elseif(this.param_BO.isGpFitGpml)
                                    switch sel_constraint.typeConstraint
                                        case TypeConstraint.UnknownMeanConstraint
                                            sel_type = TypeEstimatorGp.GpmlStochasticGivenNoiseVar;
                                        case TypeConstraint.UnknownObsConstraint
                                            sel_type = TypeEstimatorGp.GpmlStochasticGivenNoiseVar;
                                        case TypeConstraint.UnknownStdDevConstraint
                                            sel_type = TypeEstimatorGp.GpmlLogSampleVarRandomLS;
                                        otherwise
                                            throwError('[MultipleConstraintsEstimatorGP] Unsupported types');
                                    end
                                else
                                end
                            else
                                throwUndefinedTypeError(this.param_BO.typeProblem);
                            end
                        end
                    end
                    
%                     sel_type = getElementFromArrayOrCell(arr_type,i);
                    sel_gpParams = getElementFromArrayOrCell(this.gpParams,i);
                    cell_estimators{i,1} = BoFactory.getEstimatorGP(sel_type, this.param_BO, sel_gpParams, sel_constraint.getXDomain());
                end
            this.cell_estimators = cell_estimators;
        end
        
        function [multiPredictors, gpParams, multiTruePredictors] = estimateGp(this, x_train, y_train, cumul_samplesizes, varargin)
            predictors = cell(this.nConstraints, 1);
            gpParams = cell(this.nConstraints, 1);
            truePredictors = cell(this.nConstraints, 1);
            
            for i=1:this.nConstraints
                [predictors{i,1}, gpParams{i,1}, truePredictors{i,1}] = this.cell_estimators{i}.estimateGp(x_train, y_train(:,i), cumul_samplesizes, varargin{:});
            end
            multiPredictors = ConstraintsPredictorGpWrapper(predictors);
            multiTruePredictors = ConstraintsPredictorGpWrapper(truePredictors);
        end
    end
end

