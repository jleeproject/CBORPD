

setting = ObjectCopier.copy(setting_default);






% wrong_setting = false;
if(numel(type_cell_settings) ==0)
    ; % DO NOTHING
else
    custom_setting = type_cell_settings{idx_setting};
    for i = 1:numel(custom_setting)
        if(iscell(custom_setting))
            sel_setting = custom_setting{i};
        else
            sel_setting = custom_setting(i);
        end
        switch class(sel_setting)
            case 'TypeInitialSampleSize'
                error('TypeInitialSampleSize cannot be set; because initial samples need to drawn before settings (all the settings have the same initial samples).')
            case 'TypeSimulator'
                error('TypeSimulator cannot be set; because initial samples need to drawn before settings (all the settings have the same initial samples).')
            case 'TypeAcquisition'
                setting.type_acquisition = sel_setting;
            case 'TypeEstimatorGp'
                setting.type_gp_fit = sel_setting;
            case 'TypeGpMean'
                if(sel_setting == TypeGpMean.Conservative)
                    setting.is_conservative = true;
                elseif(sel_setting == TypeGpMean.Liberal)
                    setting.is_conservative = false;
                else
                    throwError('Undefined TypeGpMean');
                end
            case 'TypeInfillOptimizer'
                setting.type_infill_opt = sel_setting;
            otherwise
                error('Undefined Custom Setting.')
        end
    end
end


if(setting.optProb.nCon>0)
    hasConstraints = true;
    if(isa(setting.optProb.con,'MultipleConstraints'))
        if(contains(setting.type_infill_opt.char,'Unconst'))
            str_err = sprintf('[ERROR] Wrong Setting: Unconstraints Infill Optimizer & Constrained Problem\n');
            str_err = sprintf('%s[ERROR] Infill Optimizer = %s \n', str_err, setting.type_infill_opt.char);
            str_err = sprintf('%s[ERROR] Unconstraints Problem  \n', str_err);
            error(str_err);
            wrong_setting = true;
        end
    end
else
%% Unconstraints. Constraints are not set. empty.
    hasConstraints = false;
    if(contains(setting.type_infill_opt.char,'Const'))
            str_err = sprintf('[ERROR] Wrong Setting: Constrained Infill Optimizer & Unconstraints Problem\n');
            str_err = sprintf('%s[ERROR] Infill Optimizer = %s \n', str_err, setting.type_infill_opt.char);
            str_err = sprintf('%s[ERROR] Unconstraints Problem  \n', str_err);
            error(str_err);
        wrong_setting = true;
    end
    if(setting.type_simulator == TypeSimulator.LogSampleVar)
        if(numel(simul.meanFunc)>0)
            meanFunc = [];
            disp('[Warning] In problem of Robust Design Optimization, ');
            disp('[Warning] No constraints, but mean function was satisfied.');
            disp('[Warning] Forced no mean function.');
        end
    end
end



% setting.max_opt_eval= 500;
% setting.max_opt_iter= 100000;

% 
param_BO = Param_BO_v2(setting);



acqFunc = BoFactory.getAcquisitionFunction(setting.type_acquisition, param_BO, setting.optProb.objFunc , setting.isMaximizeObj);


% gpSettingMaker = GpSettingMaker(param_BO.getTypeGpFit);
gpSettingMaker = GpSettingMaker(param_BO.getTypeGpFit, setting.type_acquisition, param_BO.isConservative, setting.type_gp_kernel, yy_arr, yy_cons_arr, setting, xx_arr);

% gpSettingsObj = gpSettingMaker.getGpSetting4ObjFunc(setting.type_acquisition, param_BO.isConservative, setting.type_gp_kernel, yy_arr, setting.optProb.xDomain, setting, xx_arr);
% gpSettingsCon = gpSettingMaker.getGpSetting4Constraints(param_BO.isConservative, setting.type_gp_kernel, yy_cons_arr, setting.optProb.xDomain, setting, xx_arr);
gpSettingsObj = gpSettingMaker.getGpSetting4ObjFunc();
gpSettingsCon = gpSettingMaker.getGpSetting4Constraints();
param_BO.gpSettingMaker = gpSettingMaker;

estimator = BoFactory.getEstimatorGP(param_BO.getTypeGpFit(), param_BO, gpSettingsObj);
if(param_BO.hasConstraints)
    estimatorCon = MultipleConstraintsEstimatorGP(param_BO, gpSettingsCon);
    estimatorCon.addAllConstraints(setting.optProb.con);
end

infillOptimizer = BoFactory.getInfillOptimizer(setting.type_infill_opt, param_BO, acqFunc, setting.optProb);% after gpSettingMaker set in param_BO
