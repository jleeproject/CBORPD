function [funcH, hyperParams, funcCovModifiaiblePost] = gpDirectDeterministic(xx, yy, hyperParams, param_BO)

    meanFuncs = hyperParams.meanFuncs;
%     noiseVars = hyperParams.noiseVars;
    noiseVars = 1e-2;
    bwRange = hyperParams.bwRange;
    scaleRange = hyperParams.scaleRange;
    diRectOpts = hyperParams.diRectOpts;
    covFunc = hyperParams.covFunc;
    
    diRectBounds = [log(bwRange); log(scaleRange)];

    % Learn each GP Separately 
    nlmlF = @(t) normMargLikelihoodStochasticGivenNoise( exp(t(1)), exp(t(2)), xx, yy, ...
                    meanFuncs, noiseVars,covFunc);
    [~, optParams] = diRectWrapMax_module(nlmlF, diRectBounds, diRectOpts);
    %             bounds    - an n x 2 vector of the lower and upper bounds.
    %                         The first column is the lower bounds, and the second
    %                         column contains the upper bounds
    lengthscaleOpt = exp(optParams(1)); % lengthscale
    varGpOpt = exp(optParams(2)); % variance of GPs


    % Optimising for hyper-parameters ends here. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finally fit the GPs and get the Function Handles.
    [funcH, funcCovModifiaiblePost] = gpPosteriorStochasticGivenNoise([], xx, yy, ...
    lengthscaleOpt, varGpOpt, meanFuncs, noiseVars, covFunc);
%     [~, ~, ~, currFuncH] = GPRegression_log_var([], xx, yy, ...
%         bwOpt, scaleOpt, meanFuncs, noiseVars, covFunc);

%     funcHs = currFuncH;

    if(param_BO.isDebug())
%     fprintf('Picked BWs: %s, (%0.4f,%.4f) Scales: %s, (%0.4f,%.4f).\n', ...
        fprintf('> diRect: Selected Lengthscale: %s, (%0.4f,%.4f) GP Variance: %s, (%0.4f,%.4f).\n', ...
        mat2str(round(1e4*lengthscaleOpt)/1e4), bwRange(1), bwRange(2), ...
        mat2str(round(1e4*varGpOpt)/1e4), scaleRange(1), scaleRange(2));
    end
    hyperParams.bw = lengthscaleOpt;
    hyperParams.scale = varGpOpt;
end