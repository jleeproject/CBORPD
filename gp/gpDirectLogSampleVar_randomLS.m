% NOTWORKING
% function [infillPostGp, hyperParams, truePostGp] = gpDirectLogSampleVar_randomLS(xx, yy, samplesizes, hyperParams, param_BO)
% 
%     meanFuncs = hyperParams.meanFuncs;
% %     noiseVars = hyperParams.noiseVars;
%     noiseVars = get_var_logs(samplesizes);
%     lengthscaleRange = hyperParams.bwRange;
%     gpVarRange = hyperParams.scaleRange;
%     diRectOpts = hyperParams.diRectOpts;
%     covFunc = hyperParams.covFunc;
%     
%     diRectBounds = [log(lengthscaleRange); log(gpVarRange)];
% 
%       % Learn each GP Separately 
% %     nlmlF = @(t) normMargLikelihoodLogSampleVar( lognrndFromMeanStdDev(exp(t(1)),exp(t(1)/1)), exp(t(2)) , xx, yy, ...
%     nlmlF = @(t) normMargLikelihoodLogSampleVar( exp(t(1)), exp(t(2)) , xx, yy, ...
%                     meanFuncs, noiseVars,covFunc);
%     [~, optParams] = diRectWrapMax_module(nlmlF, diRectBounds, diRectOpts);
%     %             bounds    - an n x 2 vector of the lower and upper bounds.
%     %                         The first column is the lower bounds, and the second
%     %                         column contains the upper bounds
%     lengthscaleOpt = exp(optParams(1)); % lengthscale
%     gpVarOpt = exp(optParams(2)); % variance
% 
% 
%     % Optimising for hyper-parameters ends here. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Finally fit the GPs and get the Function Handles.
%     [funcH, funcCovModifiaiblePost] = gpPosteriorLogSampleVar([], xx, yy, ...
%     lognrndFromMeanStdDev(lengthscaleOpt, lengthscaleOpt/5), gpVarOpt, meanFuncs, noiseVars, covFunc);
%     [~, ~, ~, truePostGp] = gpPosteriorLogSampleVar([], xx, yy, ...
%     lengthscaleOpt, gpVarOpt, meanFuncs, noiseVars, covFunc);
% %     [~, ~, ~, currFuncH] = GPRegression_log_var([], xx, yy, ...
% %         bwOpt, scaleOpt, meanFuncs, noiseVars, covFunc);
% 
%     infillPostGp = currFuncH;
% %     truePostGp = truePostGp;
%     
% 
%     if(param_BO.isDebug())
% %     fprintf('Picked BWs: %s, (%0.4f,%.4f) Scales: %s, (%0.4f,%.4f).\n', ...
%         fprintf('> diRect: Selected Lengthscale: %s, (%0.4f,%.4f) GP Variance: %s, (%0.4f,%.4f).\n', ...
%         mat2str(round(1e4*lengthscaleOpt)/1e4), lengthscaleRange(1), lengthscaleRange(2), ...
%         mat2str(round(1e4*gpVarOpt)/1e4), gpVarRange(1), gpVarRange(2));
%     end
%     hyperParams.bw = lengthscaleOpt;
%     hyperParams.scale = gpVarOpt;
% end