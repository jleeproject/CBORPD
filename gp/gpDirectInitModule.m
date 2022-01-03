function gpHPs = gp_direct_init(init_samplesizes, init_x, init_y, bounds, param_BO, meanFuncConservative)
params = struct();
params.acqStrategy = 'GP-UCB';
% samplesizes = init_samplesizes;
numDims = size(bounds, 1);

% rangeY
bounds;
max_eval_diRect = param_BO.getMaxEvalDiRect();
is_conservative = param_BO.isConservative();

% stdY = std(init_y);
meanY = mean(init_y);
maxY = max(init_y);
minY = min(init_y);
rangeY = range(init_y);
% params.rangeY = rangeY;

% Removed fidelity properties.

  % Hyper parameters for the GP.
  % ============================================================================
  if ~isfield(params, 'bwRange')
    params.bwRange = sqrt(numDims) * [1e-3, 10] * mean(bounds(:,2) - bounds(:,1));
  end
  if ~isfield(params, 'scaleRange')
%     params.scaleRange = [1 10] * rangeY;
    params.scaleRange = [1e-6 100] * rangeY^2;
  end
  
  
  if ~isfield(params, 'gpMeanFuncs'),
    if strcmp(params.acqStrategy, 'MF-GP-UCB') | ...
      strcmp(params.acqStrategy, 'Illus-MF-GP-UCB') | ...
      strcmp(params.acqStrategy, 'GP-UCB') 
%       priorMeanVal = maxY + 2*rangeY; % set it to be a large value. works best for UCB.
        if(is_conservative && param_BO.isTypeAcquisitionUCB())
%             priorMeanVal = minY - 2*rangeY; % set it to be a large value. works best for UCB.
            priorMeanVal = meanFuncConservative(meanY, minY, maxY, rangeY);
        else
            priorMeanVal = meanY; % Do the usual thing. works best for EI and the rest.
        end
    else
      priorMeanVal = meanY; % Do the usual thing. works best for EI and the rest.
    end
    params.gpMeanFuncs = @(t) priorMeanVal * ones(size(t,1), 1);
	  % t : Dimension of mean matrix. 
  end
%   if ~isfield(params, 'gpNoiseVars'),
% %     params.gpNoiseVars = 1e-4 * stdY^2;
%       params.gpNoiseVars =  get_var_logs(samplesizes);
%   end
  if ~isfield(params, 'multGPLearnStrategy'),
%     params.multGPLearnStrategy = 'jointLearn'; % seems to work best for all.
    params.multGPLearnStrategy = 'sepLearn'; % seems to work best for all.
  end
  if ~isfield(params, 'gpDiRectOpts')
%     params.gpDirectOpts.maxevals = 200;
    params.gpDirectOpts.maxevals = max_eval_diRect;
  end

%   % Parameters for the initial GP
%   % -----------------------------------------------------------------------------
%   if ~isfield(params, 'initGPBw')
%     params.initGPBw = params.bwRange(1);
%   end
%   if ~isfield(params, 'initGPScale') 
%     params.initGPScale = std( cell2mat(initVals) );
%   end
%   if ~isfield(params, 'initGPFuncHs')
%     [~,~,~,currFuncH] = GPRegression([], initPts, initVals, ...
%       params.initGPBw, params.initGPScale, params.gpMeanFuncs, ...
%       params.gpNoiseVars );
%     initFuncHs = currFuncH;
%     params.initGPFuncHs = initFuncHs;
%   end

  % Finally, is the budget time based ?
%   if ~isfield(params, 'budgetType')
%     params.budgetType = 'givenCost';
%   end
  if ~isfield(params, 'diRectParams')
%     diRectParams.maxevals = ceil(7 * min(5,numDims)^2 * sqrt(min(iter, 1000)));
%     diRectParams.maxevals = 200;
    diRectParams.maxevals = max_eval_diRect;
    diRectParams.maxits = inf;
%     fprintf('t = %d, diREctEvals: %d\n', boIter, diRectParams.maxevals);
  end
  params.diRectParams = diRectParams;


gpHPs.bwRange = params.bwRange;
gpHPs.scaleRange = params.scaleRange;
gpHPs.meanFuncs = params.gpMeanFuncs;
gpHPs.diRectOpts = params.gpDirectOpts;
% gpHPs.noiseVars = params.gpNoiseVars;
gpHPs.multGPLearnStrategy = params.multGPLearnStrategy;
gpHPs.zetas  = 0;
gpHPs.diRectParams = params.diRectParams;
%   if (threshExceededCounter >= NUM_THRESHOLD_EXCEEDS)
%     fprintf('  ***  Re-learning GP hyper-params since thresholds exceeded !\n');
%     gpHPs.bwRange(2) = max(1.1*gpHPs.bwRange(1), 0.9*gpHPs.bwRange(2));
%     threshExceededCounter = 0;
%   end
% 
% boQueries = params.initPts;
% boVals = params.initVals;
%   [funcHs, learnedHPs] = multipleGPRegressionML(boQueries, boVals, gpHPs);
