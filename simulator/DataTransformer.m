classdef DataTransformer
    
    properties
        sampleTransformer
        multiConsObj
        hasConstraints = false;
        typeProb = TypeProblem.RobustDesignOptimization;
        samples
        sampleCon
    end
    
    methods
        function this = DataTransformer(multiConsObj, samples, samplesize, typeProblem, sampleCon)
            this.multiConsObj = multiConsObj;
            if(numel(multiConsObj)>0)
                this.hasConstraints = true;
            else
                this.hasConstraints = false;
            end
            this.sampleTransformer = SamplesTransformer(samples, samplesize);
            this.samples = samples;
            if(nargin>3)
                this.typeProb = typeProblem;
            end
            if(nargin>4)
                this.sampleCon = sampleCon;
            end
        end
        
        function out = getDataForObjFunc(this)
            switch this.typeProb
                case TypeProblem.RobustDesignOptimization
                    out = this.sampleTransformer.getUnbiasSampleLogVar();
                case TypeProblem.ConventionalBoDeterministic
%                     out = this.sampleTransformer.getSamples();
                    out = this.samples;
                case TypeProblem.ConventionalBoStochastic
%                     out = this.sampleTransformer.getSamples();
                    out = this.samples;
                otherwise
                    throwUndefinedTypeError(this.typeProb);
            end
        end
        
        function out = getDataForConstraints(this)
            if(this.hasConstraints)

                switch this.typeProb
                    case TypeProblem.ConventionalBoDeterministic
                        out = this.sampleCon;
                    case TypeProblem.ConventionalBoStochastic
                        out = this.sampleCon;
                    case TypeProblem.RobustDesignOptimization
                    nCons = this.multiConsObj.nConstraints;
                    cell_out = zeros(1, nCons);
                    for i=1:nCons
                        if(this.typeProb == TypeProblem.ConventionalBoDeterministic)
                            cell_out(1,i) = this.samples;
                        elseif(this.typeProb == TypeProblem.ConventionalBoStochastic)
                            cell_out(1,i) = this.samples;
                        else
                            sel_con = getElementFromArrayOrCell(this.multiConsObj.constraints,i);
                            type_con = sel_con.typeConstraint;
                            switch type_con
                                case TypeConstraint.KnownMeanConstraint
                                    val = this.sampleTransformer.getSampleMean();
                                case TypeConstraint.KnownObsConstraint
                                    val = this.sampleTransformer.getSamples();
                                case TypeConstraint.KnownStdDevConstraint
                                    val = this.sampleTransformer.getUnbiasSampleLogVar();
                                case TypeConstraint.UnknownMeanConstraint
                                    val = this.sampleTransformer.getSampleMean();
                                case TypeConstraint.UnknownObsConstraint
        %                             val = this.sampleTransformer.getSamples();
                                    val = this.samples;
                                case TypeConstraint.UnknownStdDevConstraint
                                    val = this.sampleTransformer.getUnbiasSampleLogVar();
                                otherwise
                                    error('[DataTransformer] Undefined class');
                            end
                            if(numel(val)>0)
                                cell_out(1,i) = val;
                            end
                        end
                    end
                    out = cell_out;
                otherwise
                    throwUndefinedTypeError(this.typeProb);
                end
            else
                error('[DataTransformer] no constraints are set');
            end
        end
        
        function out = getOriginalData(this)
            out = this.sampleTransformer.getSamples();
        end
        
    end
end

