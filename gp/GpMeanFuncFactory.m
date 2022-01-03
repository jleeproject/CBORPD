classdef GpMeanFuncFactory
   
    
    methods (Static)
        function out = getGpMeanForObjFunc(conservative,typeAcquisition, isMaximize)
            if(conservative)
                switch typeAcquisition
                    case TypeAcquisition.EI
                        out = GpMeanFuncFactory.getGpMean(TypeGpMeanConservatives.Mean);
                    case TypeAcquisition.UCB
                        out = GpMeanFuncFactory.getGpMean(TypeGpMeanConservatives.MaxPlus2Sigma);
                    case TypeAcquisition.LCB
                        out = GpMeanFuncFactory.getGpMean(TypeGpMeanConservatives.MinMinus2Sigma);
                    otherwise 
                        throwUndefinedTypeError(GpMeanFuncFactory, typeAcquisition)
                end
            else
                out = GpMeanFuncFactory.getGpMean(TypeGpMeanConservatives.Zero);
            end
        end
        
        function out = getGpMeanForConsFunc(conservative)
            out = GpMeanFuncFactory.getGpMean(TypeGpMeanConservatives.Mean);
        end
        
        function out = getGpMean(typeGpMean)
            switch typeGpMean
                case TypeGpMeanConservatives.Max
                    out = @(meanV, minV, maxV, rangeV) maxV;
                case TypeGpMeanConservatives.MaxPlus1Sigma
                    out = @(meanV, minV, maxV, rangeV) maxV + 1* rangeV;
                case TypeGpMeanConservatives.MaxPlus2Sigma
                    out = @(meanV, minV, maxV, rangeV) maxV + 2* rangeV;
                case TypeGpMeanConservatives.Mean
                    out = @(meanV, minV, maxV, rangeV) meanV;
                case TypeGpMeanConservatives.Min
                    out = @(meanV, minV, maxV, rangeV) minV;
                case TypeGpMeanConservatives.MinMinus1Sigma
                    out = @(meanV, minV, maxV, rangeV) minV - 1* rangeV;
                case TypeGpMeanConservatives.MinMinus2Sigma
                    out = @(meanV, minV, maxV, rangeV) minV - 2* rangeV;
                case TypeGpMeanConservatives.Zero
                    out = @(meanV, minV, maxV, rangeV) 0;
                otherwise 
                    throwUndefinedTypeError(GpMeanFuncFactory, typeGpMean)
            end
        end
    end
end

