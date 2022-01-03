%% ----------------------------------------------------------------------------
% noplot_opt = true;
noplot_opt = false;
%% ----------------------------------------------------------------------------
% from_bayesopt = true;         from_feasopt = false;
setting.show_plot_online = true;
setting.grid_bin = 35;
niters = 500;
param_BO.hasConstraints = true;
visualizer = Visualizer_BO();
visualizer.setShowAllGrid(true);
visualizer.setShowSearchTrails(true);


setting.typeProblem = TypeProblem.ConventionalBoStochastic;
show_exact_fmin=true;
optProb.objFunc.getName = 'ADMM';

infillOptimizer.acqFuncName = 'ADMM';

param_BO.info = true;

if(~exist('storage','var'))
    storage = struct();
end

if(numel(fieldnames(storage))==0)
    storage.xMuMinHistory = zeros(niters, opt.dims);
    storage.evalFAtXnewHistory = zeros(niters, 1);
    storage.yMuMinHistory =  zeros(niters, opt.dims);
    storage.samplesizeHistory = zeros(niters, 1);
    if(exist('from_bayesopt','var') && from_bayesopt)
        storage.xHistory = botrace.samples;
        storage.yHistory = botrace.values;
        storage.yConsHistory = con_values;
        nEvaluations = size(botrace.values, 1);
        storage.nEvaluations = nEvaluations;
    end
else
    nEvaluations = storage.nEvaluations+1;
end
if(exist('from_bayesfeas','var') && from_bayesfeas)
    if(~exist('storage_feas','var'))
        storage.xHistory(nEvaluations:nEvaluations+size(samples, 1)) = samples;
    %     storage.yHistory = values;
        storage.yConsHistory(nEvaluations:nEvaluations+size(samples, 1)) = C_values;
        nEvaluations = nEvaluations + size(samples, 1);
    end
end



infillOptimizer.isMaximizeAcq = true;
infillOptimizer.isMaximizeObj = false;

typeSamplesize = TypeSampleSize.Fixed;
givenSamplesize = 1;

type_infill_opt = TypeInfillOptimizer.DirectConstAdmm;

if(exist('from_bayesopt','var') && from_bayesopt && ~noplot_opt)
n_init_sample = numel(init_pt);
infillOptimizer.strStep = sprintf('opt (%d/%d)', i+n_init_sample, opt.max_iters);
infillOptimizer.name = sprintf('ADMM (z=%s)',num2str(zbar));
    if(~exist('storage_opt','var'))
        storage_opt.xMuMinHistory = zeros(niters, opt.dims);
        storage_opt.evalFAtXnewHistory = zeros(niters, 1);
        storage_opt.yMuMinHistory =  zeros(niters, opt.dims);
        storage_opt.samplesizeHistory = zeros(niters, 1);
    end
    % optProb.objFunc.getFnEval = fn_obj_true;
    optProb.objFunc.getXDomain = [opt.mins;opt.maxes]';
%     optProb.con.constraints{1} = @(x) 0;
%     optProb.conFunc{1}.fnEval = @(x) 0;
    optProb.nCon = numel(opt_all.c);
    optProb.con.constraints = cell(optProb.nCon,1);
%     optProb.con.constraints{1}  = problem.C{1};%
    optProb.conFunc = cell(optProb.nCon,1);
    for idx_con=1:optProb.nCon
    	optProb.con.constraints{idx_con} = Constraint(TypeConstraint.UnknownObsConstraint, [], 0);
        optProb.conFunc{idx_con}.evalWithVecX = problem.C{idx_con};
%         optProb.conFunc{idx_con}.getFnEval.evalWithVecX = problem.C{idx_con};
    end
    optProb.conFunc{1}.getName = '';
    optProb.objFunc.getOptSol = zbar; % selected z of ~
    optProb.objFunc.evalWithVecX = fn_obj_true;
    %% Visualizer.
    visualizer.setFnEval(fn_obj_true)
    visualizer.setGrid(setting.grid_bin, optProb.objFunc)

    
    %
    queries_hist.mu = mu;
    queries_hist.acq = ei;
    queries_hist.all = hyper_grid_prev;
    queries_hist.improv = [];
    queries_hist.sigma = sqrt(sigma2);


    iter = n_init_sample+i;
%     iter = nEvaluations;


    storage_opt.xHistory = botrace.samples;
    storage_opt.yHistory = botrace.values;
    storage_opt.yConsHistory = con_values;
    
    
    storage_opt.samplesizeHistory = ones(iter+i,1);
    storage_opt.xMuMinHistory(iter,:) = zbar;
    storage_opt.evalFAtXnewHistory(iter,1) = fn_obj_true(hyper_cand);
    storage_opt.yMuMinHistory(iter,1) = max(mu);
    
		[mv,mi] = min(values);
		minvalue = mv;
		minsample = unscale_point(samples(mi,:),opt.mins,opt.maxes);

    storage.xHistory(nEvaluations,:) = botrace.samples(end,:);
    storage.yHistory(nEvaluations,:) = botrace.values(end,:);
%     storage.yConsHistory(nEvaluations,:) = con_values(end,:);
    storage.samplesizeHistory = ones(nEvaluations,1);
    storage.xMuMinHistory(nEvaluations,:) = minsample;
    storage.evalFAtXnewHistory(nEvaluations,1) = fn_obj_true(minsample);
    storage.yMuMinHistory(nEvaluations,1) = max(mu);

    % %% NEED TO ASSIGN
    funcFeas = @(x,urnd) voidFunc2out(x);
    multiPredictors.cell_predictors{1}.predict = @(x) funcFeas(x,1);
    % 
    % %%
    exactPredictor.predict = @(xx) getScaledMuSigma(xx, out_fn_posterior, opt);
    acquisition = @(xx) getScaledAcq(xx, out_fn_posterior, out_fn_acq, opt);
    
    
    if(numel(gpMeanVar)>0)
        multiPredictors.cell_predictors{1}.predict = gpMeanVar;
        funcFeas = @(xx, urnd)  isEstFeasible(gpMeanVar, xx, urnd);
    end
    
    % 
% 	storage = storage_opt;    %     
    % 
elseif(exist('from_feasopt','var') && from_feasopt)
    % optProb.objFunc.getFnEval 
% infillOptimizer.strStep = 'feas';
n_init_sample = 2;
infillOptimizer.strStep = sprintf('feas (%d/%d)', i+n_init_sample, opt.c{1}.max_iters);
if(~exist('storage_feas','var'))
        storage_feas.xMuMinHistory = zeros(niters, opt.f.dims);
        storage_feas.evalFAtXnewHistory = zeros(niters, 1);
        storage_feas.yMuMinHistory =  zeros(niters, opt.f.dims);
        storage_feas.samplesizeHistory = zeros(niters, 1);
    end
infillOptimizer.name = sprintf('ADMM (z=%s)',num2str(opt.f.x));
fn_obj_true = problem.F;
    
    % optProb.objFunc.getFnEval 
    optProb.objFunc.getXDomain = [opt.f.mins;opt.f.maxes]';;
    
%     optProb.con.constraints{1}  = problem.C{1};%
%     optProb.conFunc{1}.fnEval.evalWithVecX = problem.C{1};
%     optProb.conFunc{1}.getFnEval.evalWithVecX = problem.C{1};
    optProb.objFunc.getOptSol = opt.f.x;%
    optProb.objFunc.evalWithVecX  = problem.F;

    optProb.nCon = numel(opt.c);
    optProb.con.constraints = cell(optProb.nCon,1);
% 	optProb.con.constraints{1} = Constraint(TypeConstraint.UnknownObsConstraint, [], 0);
    
    optProb.conFunc = cell(optProb.nCon,1);
    for idx_con=1:optProb.nCon
    	optProb.con.constraints{idx_con} = Constraint(TypeConstraint.UnknownObsConstraint, [], 0);
        optProb.conFunc{idx_con}.evalWithVecX = problem.C{idx_con};
%         optProb.conFunc{idx_con}.getFnEval.evalWithVecX = problem.C{idx_con};
    end
    optProb.conFunc{1}.getName = '';
%     optProb.conFunc{1}.getFnEval  = problem.F;
    %% Visualizer.
	visualizer.setFnEval(fn_obj_true);
    visualizer.setGrid(setting.grid_bin, optProb.objFunc);


    %
    queries_hist.mu = m;
    queries_hist.acq = EI;
    queries_hist.all = grid_data;
    queries_hist.improv  = [];
    queries_hist.sigma = sqrt(s2);


    storage_feas.xHistory = samples;
    storage_feas.yHistory = [];
    storage_feas.yConsHistory = C_values;
    storage_feas.samplesizeHistory = ones(size(samples,1),1);
%     nEvaluations = size(samples,1); 
    iter = size(samples,1);


    value_opt = opt.f.x;
%     value_opt = z_opt;
    storage_feas.xMuMinHistory(iter,:) = value_opt;
    storage_feas.evalFAtXnewHistory(iter,1) = problem.F(value_opt);
    storage_feas.yMuMinHistory(iter,1) = optProb.objFunc.evalWithVecX(value_opt);
    
    
    storage.xHistory(nEvaluations,:) = samples(end,:);
%     storage.yHistory = [];
    storage.yConsHistory(nEvaluations,:) = C_values(end,:);
    storage.xMuMinHistory(nEvaluations,:) = value_opt;
    storage.evalFAtXnewHistory(nEvaluations,1) = problem.F(value_opt);
    storage.yMuMinHistory(nEvaluations,1) = optProb.objFunc.evalWithVecX(value_opt);
    storage.samplesizeHistory = ones(nEvaluations,1);

    % %% NEED TO ASSIGN
    multiPredictors.cell_predictors{1}.predict = gpMeanVar;
    funcFeas = @(xx, urnd)  isEstFeasible(gpMeanVar, xx, urnd);
    % 
    % %%
	
    out_fn_acq = @(xx,m,s2)EI_Z(m,s2,xx,h_plus,opt,const_num);
%     exactPredictor.predict = @(xx) getMuSigma(xx, gpMeanVar);
    acquisition = @(xx) getAcq(xx, gpMeanVar, out_fn_acq);
	
    exactPredictor.predict = @(xx) getScaledMuSigma(xx, out_fn_posterior, opt.f);
% 	optProb.con.constraints{1} = MultipleConstraints({Constraint(TypeConstraint.UnknownObsConstraint, [], 0)});
% 	maxEIValues

	
% 	storage = storage_feas;    %     
else
end
if ~noplot_opt || from_feasopt
    param_BO.isInfillSearchGrid = false;
    scr_visualize_BO_online_multid_v4;
    
    storage.nEvaluations = nEvaluations;
end
function [out1, out2] = voidFunc2out(in1)
    out1 = zeros(size(in1,1),1);
    out2 = zeros(size(in1,1),1);
end

function [mu, sigma] = getScaledMuSigma(xx, out_fn_posterior, opt)
    xx = scale_point(xx, opt.mins, opt.maxes);
    [mu, sigma2] = out_fn_posterior(xx);
    sigma = sqrt(sigma2);
end
function out = getScaledAcq(xx, out_fn_posterior, out_fn_acq, opt)
    xx = scale_point(xx, opt.mins, opt.maxes);
    [mu, sigma2] = out_fn_posterior(xx);
    out = out_fn_acq(xx,mu, sigma2);
end

function [mu, sigma] = getMuSigma(xx, out_fn_posterior)
    [mu, sigma2] = out_fn_posterior(xx);
    sigma = sqrt(sigma2);
end
function out = getAcq(xx, out_fn_posterior, out_fn_acq)
    [mu, sigma2] = out_fn_posterior(xx);
    out = out_fn_acq(xx,mu, sigma2);
end

function upt = unscale_point(x,mins,maxes)
    if size(x,1) == 1,
        upt = x .* (maxes - mins) + mins;
    else
        upt = bsxfun(@plus,bsxfun(@times,x,(maxes-mins)),mins);
    end
end

function pt = scale_point(x,mins,maxes)
    pt = bsxfun(@rdivide,bsxfun(@minus,x,mins),maxes-mins);
end

function [out_feas,pof] = isEstFeasible(funcPred, xx, urnd)
    [mu, sigma] = funcPred(xx);
    pof = normcdf( -mu./sigma);
%     for i=1:this.nConstraints
%         if(this.constraints{i}.hasLb&&this.constraints{i}.hasUb)
%             lb = this.constraints{i}.lb;
%             ub = this.constraints{i}.ub;
%             ps(:,i) = normcdf( (ub-mus(:,i))./sigmas(:,i)) - normcdf( (lb-mus(:,i))./sigmas(:,i));
%         elseif(this.constraints{i}.hasLb)
%             lb = this.constraints{i}.lb;
%             ps(:,i) = normcdf( 1 - normcdf( (lb-mus(:,i))./sigmas(:,i)) );
%         elseif(this.constraints{i}.hasUb)
%             ub = this.constraints{i}.ub;
%             ps(:,i) = normcdf( (ub-mus(:,i))./sigmas(:,i)  ) ;
%         else
%             throwError('Both Lb and Ub are not specified.');
%         end
%     end
    feas = urnd < pof;
    out_feas = 1 -2.*(feas);
end



function [ EI] = EI_Z(mean_grid,variance_grid,xx,h_plus,opt,const_num)
%% The goal of this function is to evaluate EI(.) for the grid data

    % preprocessing 
    Qz=@(Z) h_plus-(opt.ADMM.rho./(2.*opt.ADMM.M)).*sum((opt.f.x-Z+(opt.f.y{const_num}./opt.ADMM.rho)).^2,2);
%     Qz_val=zeros(opt.c{const_num}.grid_size,1);
%     
%     for j=1:opt.c{const_num}.grid_size
%         Qz_val(j,1)=Qz(grid_data(j,:));
%     end

    Qz_val = Qz(xx);
    
    % finding the EI() regions
    ind_2=find(Qz_val>=0 & Qz_val<1);
    ind_3=find(Qz_val>=1);
    
    % Evaluating EI(z) based on regions
    cdf_at_0=normcdf(0,mean_grid,variance_grid);
    EI=zeros(length(Qz_val),1);
    EI(ind_2)=Qz_val(ind_2).*cdf_at_0(ind_2); 
    EI(ind_3)=Qz_val(ind_3)-(1-cdf_at_0(ind_3));
    
end


