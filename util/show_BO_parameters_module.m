if(param_BO.isInfo())
    fprintf('Infill Method = %s\n',class(infillOptimizer) );
end


fprintf('[[ START TIME : %s ]]\n', showPrettyDateTime(now(),'yyyy/mm/dd HH:MM:SS') );
fprintf('[[ REPEAT : %d/%d]]\n', repeat, nRepeats);
fprintf('[[ obj function : ''%s'' ]]\n', setting.optProb.objFunc.name);
fprintf('[[ Problem : ''%s'' ]]\n', setting.typeProblem);
% fprintf('[[ Constraints : ''%s'' ]]\n', constraints);
if(exist('meanFunc','var') && numel(meanFunc)>0)
    fprintf('[[ Mean (Constraint) function : ''%s'' ]]\n', meanFunc.name);
end
if(typeSamplesize == TypeSampleSize.Adjust)
    fprintf('[[ Sample Size : Adjusted ]]\n');
else
    fprintf('[[ Sample Size : %d ]]\n', givenSamplesize);
end


switch param_BO.getTypeAcquisition()
    case TypeAcquisition.EI
        fprintf('[[ Acquisition Type : %s]]\n', 'EI');
    case TypeAcquisition.LCB
        fprintf('[[ Acquisition Type : %s]]\n', 'LCB');
    case TypeAcquisition.UCB
        fprintf('[[ Acquisition Type : %s]]\n', 'UCB');
    otherwise
        throwError('Undefined type');
end

switch param_BO.is_conservative
    case true
        fprintf('[[ Conservative ]]\n');
    case false
        fprintf('[[ No Conservative ]]\n');
    otherwise
        throwError('Undefined type');
end

switch param_BO.getMethodSearch
    case TypeInfillOptimizer.DirectConstSampling
        fprintf('[[ Search Method : %s]]\n', 'diRect: Sampling [Constraints]');
    case TypeInfillOptimizer.DirectUnconst
        fprintf('[[ Search Method : %s]]\n', 'diRect');
    case TypeInfillOptimizer.GridConstSampling
        fprintf('[[ Search Method : %s]]\n', 'Grid (GPML)');
    case TypeInfillOptimizer.GridUnconst
        fprintf('[[ Search Method : %s]]\n', 'Grid (GPML)');
    case TypeInfillOptimizer.DirectConstWEI
        fprintf('[[ Search Method : %s]]\n', 'diRect: WEI [Constraints]');
    case TypeInfillOptimizer.DirectConstTesting
        fprintf('[[ Search Method : %s]]\n', 'diRect: Testing[Constraints]');
    case TypeInfillOptimizer.DirectConstTesting2
        fprintf('[[ Search Method : %s]]\n', 'diRect: Testing2[Constraints]');
    otherwise
        fprintf('[[ Search Method : %s]]\n', param_BO.getMethodSearch.char);        
%         error('Undefined type');
end



switch param_BO.getTypeGpFit()
    case TypeEstimatorGp.DirectLogSampleVar
        fprintf('[[ GP Est Method : %s]]\n', 'diRect');
    case TypeEstimatorGp.GpmlLogSampleVar
        fprintf('[[ GP Est Method : %s]]\n', 'GPML');
    case TypeEstimatorGp.DirectLogSampleVarRandomLS
        fprintf('[[ GP Est Method : %s]]\n', 'diRect with random Lengthscale');
    case TypeEstimatorGp.GpmlLogSampleVarRandomLS
        fprintf('[[ GP Est Method : %s]]\n', 'GPML with random Lengthscale');
    otherwise
        if(param_BO.isGpFitGpml)
            fprintf('[[ GP Est Method : %s]]\n', 'GPML');
        elseif(param_BO.isGpFitDirect)
            fprintf('[[ GP Est Method : %s]]\n', 'diRect');
        else
            throwError('Undefined type');
        end
end

if(param_BO.isInfillSearchDirect || param_BO.isGpFitDirect)
    fprintf(' - ''diRect'' [[ Max Evaluation : %d]]\n', param_BO.getMaxOptIter());
end

fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n');
fprintf('[[ # of initial samples : %d ]]\n', param_BO.getNInitSample());
fprintf(' - ''inital samples'' [[ Sample size : %d ]]\n', param_BO.getSamplesizeForInitStage());
fprintf('[[ Sample Budget : %d ]]\n', budget);
fprintf('[[ Power : %.4g ]]\n', power);
fprintf('--------------------------------------------------------------\n');
