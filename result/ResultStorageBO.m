classdef ResultStorageBO < handle
    %_CLASS_TEMPLATE Summary of this class goes here
    %   Detailed explanation goes here
    
%     properties (Constant)
%         Property1
%     end
%     
%     properties (Abstract)
%         Property1
%     end
    
    properties
        param_BO
        objFunc
        niter
%         x_history 
%         y_history 
%         eval_f_at_xnew_history 
%         samplesize_history 
%         mu_min_history 
%         xmin_history 
%         cumElapsedTime 
%         iterElapsedTime 
% 
%         hist_beta_t 
        xHistory 
        yHistory 
        yConsHistory
        samplesizeHistory 

        evalFAtXnewHistory 
        yMuMinHistory
        xMuMinHistory
        modeHistory % 1: opt, 0: feas
        hasMode = false;
        
        evalFAtXnewLiberalHistory 
        feasAtXnewHistory
        yMuMinLiberalHistory
        xMuMinLiberalHistory
        hasLiberal = false;
        evalGXnewHistory;
        feasgapAtXnewHistory;

        obsHistory
%         muMinHistory 
%         xminHistory 
        cumElapsedTime 
        iterElapsedTime 
        histBetaT 
        typeSamplesize
        givenSamplesize
        
        elapsedTotalExp = 0;
        elapsedTotalEval = 0;
        elapsedTotalBoWoEval = 0;
        
        %% For ADMMBO
        firstOpt=true;
        firstFeas=true;
        nEvaluations = 0;
        n_init_sample = 0;
    end
    
%     methods (Abstract)
%         method2(this)
%     end
    methods
        function this = ResultStorageBO(param_BO, objFunc, typeSamplesize, givenSamplesize)
            % Construct an instance of this class
            this.param_BO = param_BO;
            this.objFunc = objFunc;
            init(this, typeSamplesize, givenSamplesize);
            this.typeSamplesize = typeSamplesize;
            this.givenSamplesize = givenSamplesize;
        end
        
        function init(this, typeSamplesize, givenSamplesize)
            switch typeSamplesize
                case TypeSampleSize.Adjust
                    niter = this.param_BO.getEstUpperLimitNumIter(typeSamplesize, 1);
                case TypeSampleSize.Fixed
                    niter = this.param_BO.getEstUpperLimitNumIter(typeSamplesize, givenSamplesize);
                otherwise
                    throwError('[ResultStorageBO] Unknown typeSamplesize');
            end
            if(isa (this.objFunc , 'GrapheneModelSolver'))
                dimEvalFunc = size(this.objFunc.decisionVariables.domains.all,1);
            else
                dimEvalFunc = this.objFunc.getDim();
            end

            this.xHistory  = zeros(niter, dimEvalFunc);
            this.yHistory = zeros(niter, 1);
            this.yConsHistory = zeros(niter, this.param_BO.nConstraints);
            this.evalFAtXnewHistory = zeros(niter,1);
            this.feasAtXnewHistory = zeros(niter,1);
            this.samplesizeHistory = zeros(niter, 1);
            this.yMuMinHistory = zeros(niter, 1);
            this.xMuMinHistory= zeros(niter, dimEvalFunc);
            
            
            
            this.evalGXnewHistory = zeros(niter,this.param_BO.nConstraints);
            this.feasgapAtXnewHistory = zeros(niter,this.param_BO.nConstraints);

            
            this.obsHistory = cell(niter,1);
            
            this.cumElapsedTime = zeros(niter,1);
            this.iterElapsedTime = zeros(niter,1);

            this.histBetaT =  zeros(niter, 1);
            this.niter = niter;
        end
        
        function addElapsedExp(this, time)
            this.elapsedTotalExp = this.elapsedTotalExp + time;
        end
        
        function addElapsedEval(this, time)
            this.elapsedTotalEval = this.elapsedTotalEval + time;
        end
        
        function addElapsedBoWoEval(this, time)
            this.elapsedTotalBoWoEval = this.elapsedTotalBoWoEval + time;
        end
        
        
        %--------------------- [xHistory] ---------------------
        function saveXHistory(this, value, index)
            this.xHistory(index,:) = value;
        end

        function out = getXHistoryTill(this, index)
            out = this.xHistory(1:index,:);
        end

        %--------------------- [yHistory] ---------------------
        function saveYHistory(this, value, index)
            this.yHistory(index,:) = value;
        end
        
        function saveYConsHistory(this, value, index)
            this.yConsHistory(index,:) = value;
        end
    
        function out = getYHistoryTill(this, index)
            out = this.yHistory(1:index,:);
        end
        
        function out = getYConsHistoryTill(this, index)
            out = this.yConsHistory(1:index,:);
        end
        
        function saveObsHistory(this, value, index)
            if(numel(index)>1)
                for i=1:numel(index)
                    this.obsHistory{index(i),1} = value{i};
                end
            else
                this.obsHistory{index,1} = value;
            end
        end
    
        function out = getObsHistoryTill(this, index)
            out = this.obsHistory{1:index,1};
        end
        
        function saveFeasAtXnewHistory(this, value, index)
            this.feasAtXnewHistory(index,:) = value;
        end
        
        function out = getFeasAtXnewHistoryTill(this, index)
            out = this.feasAtXnewHistory(1:index,:);
        end
        
        function saveEvalGXnewHistory( this, value , index, idx_con )
            this.evalGXnewHistory(index,idx_con) = value;
        end
        function saveFeasgapAtXnewHistory( this, value , index, idx_con )
            this.feasgapAtXnewHistory(index,idx_con) = value;
        end

        function out = getEvalGXnewHistoryTill( this, index)
            out = this.evalGXnewHistory(1:index,:);
        end
        function out = getFeasgapAtXnewHistoryTill( this, index)
            out = this.feasgapAtXnewHistory(1:index,:);
        end

        %--------------------- [evalFAtXnewHistory] ---------------------
        function saveEvalFAtXnewHistory(this, value, index)
            this.evalFAtXnewHistory(index,:) = value;
        end


        %--------------------- [samplesizeHistory] ---------------------
        function saveSamplesizeHistory(this, value, index)
            this.samplesizeHistory(index,:) = value;
        end

        function out = getSamplesizeHistoryTill(this, index)
            out = this.samplesizeHistory(1:index,:);
        end

        %--------------------- [muMinHistory] ---------------------
        function saveYMuMinHistory(this, value, index)
            this.yMuMinHistory(index,:) = value;
        end

        function out = getYMuMinHistory(this, index)
            out = this.yMuMinHistory(index,:);
        end


        %--------------------- [xminHistory] ---------------------
        function saveXMuMinHistory(this, value, index)
            this.xMuMinHistory(index,:) = value;
        end
        
        function out = getXMuMinHistory(this, index)
            out = this.xMuMinHistory(index,:);
        end


        %--------------------- [cumElapsedTime] ---------------------
        function saveCumElapsedTime(this, value, index)
            this.cumElapsedTime(index,:) = value;
        end


        %--------------------- [iterElapsedTime] ---------------------
        function saveIterElapsedTime(this, value, index)
            this.iterElapsedTime(index,:) = value;
        end


        %--------------------- [histBetaT] ---------------------
        function saveHistBetaT(this, value, index)
            this.histBetaT(index,:) = value;
        end

        
        %--------------------- [muMinLiberalHistory] ---------------------
        function saveYMuMinLiberalHistory(this, value, index)
            this.yMuMinLiberalHistory(index,:) = value;
            if(this.hasLiberal==false);this.hasLiberal=true;end;            
        end

        function out = getYMuMinLiberalHistory(this, index)
            out = this.yMuMinLiberalHistory(index,:);
        end

        %--------------------- [evalFAtXnewLiberalHistory] ---------------------
        function saveEvalFAtXnewLiberalHistory(this, value, index)
            this.evalFAtXnewLiberalHistory(index,:) = value;
            if(this.hasLiberal==false);this.hasLiberal=true;end;            
        end


        %--------------------- [xminLiberalHistory] ---------------------
        function saveXMuMinLiberalHistory(this, value, index)
            this.xMuMinLiberalHistory(index,:) = value;
            if(this.hasLiberal==false);this.hasLiberal=true;end;            
        end
        
        function out = getXMuMinLiberalHistory(this, index)
            out = this.xMuMinLiberalHistory(index,:);
        end
        %--------------------- [modeHistory] ---------------------
        function setModeHistoryOpt(this, index)
            if(~this.hasMode)
                this.hasMode = true;
            end
            this.modeHistory(index,:) = 1;
        end

        function setModeHistoryFeas(this, index)
            if(~this.hasMode)
                this.hasMode = true;
            end
            this.modeHistory(index,:) = 0;
        end


    end
end

