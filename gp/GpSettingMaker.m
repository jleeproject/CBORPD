classdef GpSettingMaker < matlab.mixin.Copyable
    properties
        typeGpFit
        
        type_acquisition
        conservative
        type_gp_kernel
        yy_obj_arr
        yy_con_arr
%         meanFuncConservatives
        xDomain
        setting
        xx_arr
    end
    
    methods
        function this = GpSettingMaker(typeGpFit, type_acquisition, conservative, type_gp_kernel, yy_obj_arr, yy_con_arr, setting, xx_arr)
            this.typeGpFit = typeGpFit;
            this.type_acquisition = type_acquisition;
            this.conservative =conservative;
            this.type_gp_kernel = type_gp_kernel;
            this.yy_obj_arr = yy_obj_arr;
            this.yy_con_arr = yy_con_arr;
%             this.meanFuncConservatives = meanFuncConservatives;
            this.xDomain = setting.optProb.xDomain;
            this.setting = setting;
            this.xx_arr = xx_arr;
        end
        
        function logic = isFitGpml(this)
%             logic = (this.typeGpFit==TypeEstimatorGp.GpmlLogSampleVar)...
%                 || (this.typeGpFit==TypeEstimatorGp.GpmlLogSampleVarRandomLS);
            logic = contains(this.typeGpFit.char, 'Gpml');
        end
        
        function logic = isFitDirect(this)
%             logic = (this.typeGpFit==TypeEstimatorGp.DirectLogSampleVar)...
%                 || (this.typeGpFit==TypeEstimatorGp.DirectLogSampleVarRandomLS);
            logic = contains(this.typeGpFit.char, 'Direct');
        end

%         function setting = getGpSetting(this, conservative, type_gp_kernel, yy_arr, meanFuncConservatives, xDomain, setting)
        function setting = getGpSetting(this)
            % -------------------------------------------------------------------------
            %% Gaussian Process Customization: Using initial samples
            if(this.isFitGpml() )
                setting = GpSettingGpml(this.typeGpFit, this.conservative, this.type_gp_kernel, this.yy_obj_arr, this.meanFuncConservatives, this.xDomain, this.setting.max_opt_iter);
            end
            if(this.isFitDirect())
                setting = GpSettingDirect(this.type_gp_kernel, this.yy_obj_arr, this.meanFuncConservatives, this.xDomain, this.setting, [], this.typeGpFit);
            end
        end
        
%         function setting = getGpSetting4ObjFunc(this, type_acquisition, conservative, type_gp_kernel, yy_obj_arr, xDomain, setting, xx_arr)
        function setting = getGpSetting4ObjFunc(this)
            % -------------------------------------------------------------------------
            meanFuncConservatives = GpMeanFuncFactory.getGpMeanForObjFunc(this.conservative,this.type_acquisition);
            %% Gaussian Process Customization: Using initial samples
            if(this.isFitGpml() )
                setting = GpSettingGpml(this.typeGpFit, this.conservative, this.type_gp_kernel, this.yy_obj_arr, meanFuncConservatives, this.xDomain, this.setting.max_opt_iter);
            elseif(this.isFitDirect())
                setting = GpSettingDirect(this.type_gp_kernel, this.yy_obj_arr, meanFuncConservatives, this.xDomain, this.setting, this.xx_arr, this.typeGpFit);
%                 GpSettingDirect(this.typeGpFit, type_acquisition, conservative, type_gp_kernel, yy_arr, meanFuncConservatives, xDomain, max_opt_iter);
            else
                error('Unkonwn Type of Fitting.');
            end
        end
        
%         function settings = getGpSetting4Constraints(this, conservative, type_gp_kernel, yy_con_arr, xDomain, setting, xx_arr)
        function settings = getGpSetting4Constraints(this)
            % -------------------------------------------------------------------------
            nConst = size(this.yy_con_arr,2);
            settings = cell(nConst,1);
            for i = 1:nConst
                %% Gaussian Process Customization: Using initial samples
                meanFuncConservatives = GpMeanFuncFactory.getGpMeanForConsFunc(this.conservative);
                if(this.isFitGpml() )
                    setting = GpSettingGpml(this.typeGpFit, this.conservative, this.type_gp_kernel, this.yy_con_arr, meanFuncConservatives, this.xDomain, this.setting.max_opt_iter);
%                     (this.typeGpFit, conservative, type_gp_kernel, yy_arr(:,i), meanFuncConservatives, xDomain, max_opt_iter);
                elseif(this.isFitDirect())
                    setting = GpSettingDirect(this.type_gp_kernel, this.yy_con_arr(:,i), meanFuncConservatives, this.xDomain, this.setting, this.xx_arr, this.typeGpFit);
%                     GpSettingDirect(this.typeGpFit, conservative, type_gp_kernel, yy_arr(:,i), meanFuncConservatives, xDomain, max_opt_iter);
                else
                    error('Unkonwn Type of Fitting.');
                end
                settings{i} = setting;
            end
        end
    end
    
end

