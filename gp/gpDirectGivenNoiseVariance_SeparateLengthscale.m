function [infillPostGp, hyperParams, funcCovModifiaiblePost] = gpDirectGivenNoiseVariance_SeparateLengthscale(xx, yy, samplesizes, hyperParams, param_BO, noiseVars)

    meanFuncs = hyperParams.meanFuncs;
%     noiseVars = hyperParams.noiseVars;
%     noiseVars = get_var_logs(samplesizes);
    lengthscaleRange = hyperParams.bwRange;
    gpVarRange = hyperParams.scaleRange;
    diRectOpts = hyperParams.diRectOpts;
    covFunc = hyperParams.covFunc;
    
    diRectBounds = [log(gpVarRange); lengthscaleRange];

      % Learn each GP Separately 
%     nlmlF = @(t) normMargLikelihoodLogSampleVar( lognrndFromMeanStdDev(exp(t(1)),exp(t(1)/1)), exp(t(2)) , xx, yy, ...
    nlmlF = @(t) normMargLikelihoodStochasticGivenNoise_SeparateLengthscale( exp(t(1)), t(2:end), xx, yy, ...
                    meanFuncs, noiseVars,covFunc);
%     diRectBounds = [log(gpVarRange); (lengthscaleRange)];
% 
%       % Learn each GP Separately 
% %     nlmlF = @(t) normMargLikelihoodLogSampleVar( lognrndFromMeanStdDev(exp(t(1)),exp(t(1)/1)), exp(t(2)) , xx, yy, ...
%     nlmlF = @(t) normMargLikelihoodLogSampleVar_SeparateLengthscale( exp(t(1)), (t(2)) , xx, yy, ...
%                     meanFuncs, noiseVars,covFunc);
    [~, optParams, history, queries, queryVals] = diRectWrapMax_module(nlmlF, diRectBounds, diRectOpts);
    %             bounds    - an n x 2 vector of the lower and upper bounds.
    %                         The first column is the lower bounds, and the second
    %                         column contains the upper bounds
    gpVarOpt = exp(optParams(1)); % lengthscale
    lengthscaleOpt = (optParams(2:end)); % variance


    % Optimising for hyper-parameters ends here. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finally fit the GPs and get the Function Handles.
%     [~, ~, ~, currFuncH] = gpPosteriorLogSampleVar([], xx, yy, ...
%     lognrndFromMeanStdDev(lengthscaleOpt, lengthscaleOpt/10), gpVarOpt, meanFuncs, noiseVars, covFunc);
    [infillPostGp, funcCovModifiaiblePost] = gpPosteriorStochasticGivenNoise_SeparateLengthscale([], xx, yy, ...
    lengthscaleOpt, gpVarOpt, meanFuncs, noiseVars, covFunc);
%     [~, ~, ~, currFuncH] = GPRegression_log_var([], xx, yy, ...
%         bwOpt, scaleOpt, meanFuncs, noiseVars, covFunc);

%     infillPostGp = truePostGp;
%     truePostGp = truePostGp;
    

    if(param_BO.isDebug())
%     fprintf('Picked BWs: %s, (%0.4f,%.4f) Scales: %s, (%0.4f,%.4f).\n', ...
%         fprintf('> diRect Selected: ls: ', ...
%         mat2str(round(1e4*lengthscaleOpt)/1e4), lengthscaleRange(1), lengthscaleRange(2), ...
%         );
        fprintf('> diRect Selected: ls: ')
        for i=1:numel(lengthscaleOpt)
            fprintf('%.2g (%.2g,%.2g), ', lengthscaleOpt(i) , lengthscaleRange(i,1), lengthscaleRange(i,2));
        end
        fprintf(' GP Variance: %.2g, (%0.2g,%.2g).\n', gpVarOpt, gpVarRange(1), gpVarRange(2));
    end
    hyperParams.bw = lengthscaleOpt;
    hyperParams.scale = gpVarOpt;
end