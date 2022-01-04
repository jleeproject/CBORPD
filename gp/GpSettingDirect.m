classdef GpSettingDirect < matlab.mixin.Copyable
    properties
        gpParam
        typeGpFit
    end
    methods
        function logic = isFitGpml(this, typeGpFit)
            logic = contains(typeGpFit.char,'Gpml');
%             logic = (typeGpFit==TypeEstimatorGp.GpmlLogSampleVar)...
%                 || (typeGpFit==TypeEstimatorGp.GpmlLogSampleVarRandomLS);
        end
        
        function logic = isFitDirect(this, typeGpFit)
            logic = contains(typeGpFit.char,'Direct');
%             logic = (typeGpFit==TypeEstimatorGp.DirectLogSampleVar)...
%                 || (typeGpFit==TypeEstimatorGp.DirectLogSampleVarRandomLS);
        end

        function this = GpSettingDirect(type_gp_kernel, yy_arr, meanFuncConservatives, xDomain, setting, init_x, typeGpFit)
            % -------------------------------------------------------------------------
            %% Gaussian Process Customization: Using initial samples
            if(exist('typeGpFit','var'))
                this.typeGpFit = typeGpFit;
            end

            gpParam = gpDirectInitModule(this, yy_arr, xDomain , meanFuncConservatives, setting, init_x, typeGpFit);
            switch type_gp_kernel 
                case TypeGpKernel.Matern52;
%                     gpParam.covFunc = @(logbw, logscale, x, x2)covMaterniso(1,[logbw, logscale],x, x2);
                    gpParam.covFunc = @gpKernelMatern52;
%                     gpParam.covFunc = @(bw, scale, x, x2)covMaterniso_gpml(5, [bw, scale], x, x2);
%                     error('Not supported kernel');
                case TypeGpKernel.Se;
                    gpParam.covFunc = @gpKernelSqExp;
                otherwise
                    error('Not supported kernel');
            end
            this.gpParam = gpParam;
        end
        
        function gpHPs = gpDirectInitModule(this, init_y, bounds, meanFuncConservative, setting, init_x, typeGpFit)
            params = struct();
            numDims = size(bounds, 1);

            % rangeY
            % stdY = std(init_y);
            meanY = mean(init_y);
            maxY = max(init_y);
            minY = min(init_y);
            rangeY = range(init_y);
            % params.rangeY = rangeY;

              % ============================================================================
            if ~isfield(params, 'bwRange')
%                 params.bwRange = sqrt(numDims) * [1e-2, 10] * mean(bounds(:,2) - bounds(:,1));
                if( contains(typeGpFit.char,'_SeparateLengthscale'));
%                     params.bwRange = [zeros(size(bounds,1),1) 1./numDims./(bounds(:,2) - bounds(:,1)).^2];
                    params.bwRange = [zeros(size(bounds,1),1) 200./numDims./(bounds(:,2) - bounds(:,1)).^2];
                else
                    params.bwRange = sqrt(numDims) * [5e-2, 10] * mean(bounds(:,2) - bounds(:,1));
                end
%                 distancePairs=distSquaredGPModule(init_x, init_x);idx=find(distancePairs);distances = distancePairs(idx);
%                 scaleRange = max(quantile(sqrt(distSquaredGPModule(distances,distances)),.05,'all'), 1e-3);
%                 maxd = max(distances);
%                 scaleRange = max(mean(distances)/50, 1e-3);
%                 params.scaleRange = [1e-6 100] * rangeY^2;
%                 params.bwRange  = [scaleRange maxd];
            end
            if ~isfield(params, 'scaleRange')
            %     params.scaleRange = [1 10] * rangeY;
%                 distancePairs=distSquaredGPModule(init_x, init_x);idx=find(distancePairs);distances = distancePairs(idx);
% %                 scaleRange = max(quantile(sqrt(distSquaredGPModule(distances,distances)),.05,'all'), 1e-3);
%                 maxd = max(distances);
%                 scaleRange = max(mean(distances)/50, 1e-3);
                params.scaleRange = [1e-6 10] * rangeY^2;
%                 params.scaleRange = [scaleRange maxd];
            end


            if ~isfield(params, 'gpMeanFuncs')
%                 switch type_acquisition
%                     case TypeAcquisition.EI
%                         if(conservative)
%                             priorMeanVal = meanFuncConservative(meanY, minY, maxY, rangeY);
%                         else
%                             priorMeanVal = meanY; % Do the usual thing. works best for EI and the rest.
%                         end
%                     case TypeAcquisition.UCB
%                         if(conservative)
%                             priorMeanVal = meanFuncConservative(meanY, minY, maxY, rangeY);
%                         else
%                             priorMeanVal = meanY; % Do the usual thing. works best for EI and the rest.
%                         end
%                     case TypeAcquisition.LCB
%                         if(conservative)
%                             priorMeanVal = meanFuncConservative(meanY, minY, maxY, rangeY);
%                         else
%                             priorMeanVal = meanY; % Do the usual thing. works best for EI and the rest.
%                         end
%                     otherwise
%                         throwError('Unkonwn Acquisition Function');
%                 end
                priorMeanVal = meanFuncConservative(meanY, minY, maxY, rangeY);
                params.gpMeanFuncs = @(t) priorMeanVal * ones(size(t,1), 1);
            end

            if ~isfield(params, 'multGPLearnStrategy'),
            %     params.multGPLearnStrategy = 'jointLearn'; % seems to work best for all.
                params.multGPLearnStrategy = 'sepLearn'; % seems to work best for all.
            end
            if ~isfield(params, 'gpDiRectOpts')
            %     params.gpDirectOpts.maxevals = 200;
                params.gpDirectOpts.maxevals = setting.max_opt_eval;
                params.gpDirectOpts.maxits = setting.max_opt_iter;
            end
              if ~isfield(params, 'diRectParams')
            %     diRectParams.maxevals = ceil(7 * min(5,numDims)^2 * sqrt(min(iter, 1000)));
            %     diRectParams.maxevals = 200;
                diRectParams.maxevals = setting.max_opt_eval;
                diRectParams.maxits = setting.max_opt_iter;
            %     fprintf('t = %d, diREctEvals: %d\n', boIter, diRectParams.maxevals);
              end
              params.diRectParams = diRectParams;


            gpHPs.bwRange = params.bwRange;
            gpHPs.scaleRange = params.scaleRange;
            gpHPs.meanFuncs = params.gpMeanFuncs;
            gpHPs.diRectOpts = params.gpDirectOpts;
            gpHPs.multGPLearnStrategy = params.multGPLearnStrategy;
            gpHPs.zetas  = 0;
            gpHPs.diRectParams = params.diRectParams;
        end        
        
    end
    
end

