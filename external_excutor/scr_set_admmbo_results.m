% fprintf('# Samples:             fprintf('[[ FINISHED ]]\n');
struct_result.nEvaluations = storage.nEvaluations;
struct_result.storage = storage;
struct_result.setting = setting;
struct_result.x_history  = storage.getXHistoryTill(nEvaluations) ;
struct_result.y_history  = storage.getYHistoryTill(nEvaluations);
struct_result.eval_f_at_xnew_history  = storage.evalFAtXnewHistory;
struct_result.samplesize_history  = storage.samplesizeHistory;
struct_result.xMuMinHistory  = storage.xMuMinHistory;
struct_result.yMuMinHistory  = storage.yMuMinHistory;
struct_result.hist_beta_t  = storage.histBetaT;
struct_result.typeSamplesize  = storage.typeSamplesize;
struct_result.givenSamplesize = storage.givenSamplesize;

struct_result.optProb = optProb;
struct_result.simul = simul;

% struct_result.nEvaluations  = nEvaluations ;        
struct_result.n_init_sample = storage.n_init_sample;
struct_result.param_BO = param_BO;
struct_result.objFunc = optProb.objFunc;

struct_result.type_cell_samplesize = type_cell_samplesize;
struct_result.cell_obj_func = cell_obj_func;
struct_result.type_cell_settings = type_cell_settings;
struct_result.nRepeats = nRepeats;
if(show_plot_online)
    struct_result.visualizer = visualizer;
else
    struct_result.visualizer = [];
end
struct_result.budget = budget;
struct_result.elapsedTime = toc(stt_experiment);
struct_result.cumElapsedTime = storage.cumElapsedTime;
struct_result.iterElapsedTime = storage.iterElapsedTime;
struct_result.opt_val_fn = optProb.objFunc.getOptVal();
struct_result.opt_sol_fn = optProb.objFunc.getOptSol();

if(exist('power','var'))
    struct_result.power = power;
end

if(exist('funcFeas','var'))
    struct_result.funcFeas = funcFeas;
end

struct_result.hasLiberal = storage.hasLiberal;
struct_result.eval_f_at_xnew_liberal_history = storage.evalFAtXnewLiberalHistory;
struct_result.xMuMinLiberalHistory = storage.xMuMinLiberalHistory;
struct_result.yMuMinLiberalHistory = storage.yMuMinLiberalHistory;

struct_result.cell_mean_func = cell_mean_func;
struct_result.nTypeSettings = nTypeSettings;
struct_result.xDomain = xDomain;
struct_result.simulator = simulator;
struct_result.initializer = initializer;
struct_result.acqFunc = acqFunc;
struct_result.infillOptimizer = infillOptimizer;
struct_result.obj = simul.obj;
struct_result.con = simul.con;
struct_result.gpSettingMaker = gpSettingMaker;
struct_result.gpSettingsObj =gpSettingsObj ;
struct_result.gpSettingsCon =gpSettingsCon;
struct_result.estimator =estimator;
if(param_BO.hasConstraints)
    struct_result.estimatorCon = estimatorCon;
end
struct_result.samplesize_cell = samplesize_cell;
struct_result.samplesizeSelector = samplesizeSelector;

results_all{idx_samplesize,repeat, idx_function, idx_setting} = struct_result;

clear funcFeas storage nEvaluations struct_result ;
