% clear all;clc; close all;
nSensor = 10;
%% %
% % filename = '1.Random_W_ratio_1_2';
% % rndm_type = 2; % Random W
% % param_sig_horizontal = 1;
% % sigma2mean_scale = 1/2;
% %% %
% % filename = '2.Random_W_ratio_1_100';
% % rndm_type = 2; % Random W
% % param_sig_horizontal = 10;
% % sigma2mean_scale = 1/100;
% % %
% filename = '3.Random_All_Chng_Params';
% rndm_type = 1; % Random chng Parameters
% param_sig_horizontal = 10;
%
% rndm_type = 1; % Random chng Parameters
% rndm_type = 2; % Random W
% draw_fig = show_progress;
draw_fig = true;
% draw_fig = false;
if(draw_fig)
    fig1 = 1;
    fig2 = 17;
    fig3 = 18;
    figure(fig1);clf;
    figure(fig2);clf;
    figure(fig3);clf;
    weights = [0.8, 0.9, 1, 1.1, 1.2];
    marks2 = {'k--', 'r-', 'b:', 'g--'};
    marks = {'ko', 'ro', 'bo', 'go'};
end
cont = (0:10:20)';
% params_graphene_ext_LW_v3;
Param_used;


load datanew
% params_graphene_ext_LW_modified;
chngParams_all = {chngParams_1, chngParams_2, chngParams_3, chngParams_4};
    true_yys_all = {datanew(1:4:161,3), datanew(1:4:161,4), datanew(1:4:161,5), datanew(1:4:161,2)};
% Vbs= [-80:8:80];
% vg = datanew(1:2:161,1);
% vg = datanew(1:1:161,1);
vg = (-40 :1:80)';
dim_v = max(size(vg));
verbose = false;
nRepeat = nSensor;

nCont = 3;
nSensor = nRepeat;
    obs = zeros(dim_v, nCont, nSensor);
    std_obs = zeros(dim_v, nCont, nSensor);
std_real_obs = zeros(dim_v, nCont);
fit_real_obs = zeros(dim_v, nCont);
real_obs = [true_yys_all{1} true_yys_all{2} true_yys_all{3}];

sig_vb0 = param_sig_horizontal;    
sig_err = .7e-7;

gen_with_extended_regression = true;
gen_with_est_parameter = false;

sig_vb2 = chngParams_all{1}(2)/200;
sig_vb3 = chngParams_all{1}(3)/200;
sig_vb4 = chngParams_all{1}(4)/500;

dim_c = max(size(cont));
params_gen = zeros(dim_c,4);
betas = zeros(2,4);

    params_trend = [chngParams_all{1}; chngParams_all{3}];
for idx = 1:4;
    betas(:,idx) = regress(params_trend(:,idx),[1 1;0 20]');
    params_gen(:,idx) = betas(1,idx) + betas(2,idx)*cont;
end
if(rndm_type==2)
    mean_W = W_orig;
    sig_W = W_orig*sigma2mean_scale;
else
    mean_W = W_orig;
    sig_W = 0;
end
generator = @(contam_lv) generate_test_data_emulator_v3(contam_lv, mean_W, sig_W, vg, fixedParams,betas, constants, sig_vb0, sig_err , obs, rndm_type, sig_vb2, sig_vb3, sig_vb4);
                figure(fig3);

for repeat = 1:nRepeat
    for idx_chng_chngParams = 1:dim_c
        fprintf('[Gen Iter] %d/%d\n',repeat, nRepeat)
        obs(:,idx_chng_chngParams,repeat) = generator(cont(idx_chng_chngParams));
    end
end
% if(gen_with_extended_regression)
% 
% %     cont = (0:10:30)';
%     dim_c = max(size(cont));
% 
% 
%     for repeat = 1:nRepeat
% 
%         fprintf('[Gen Iter] %d/%d\n',repeat, nRepeat)
%         for idx_chng_chngParams = 1:dim_c
%             if(rndm_type==1)
%                 dev_vb0 = normrnd(0,sig_vb0);
%                 dev_param2 = normrnd(0,sig_vb2);
%                 dev_param3 = normrnd(0,sig_vb3);
%                 dev_param4 = normrnd(0,sig_vb4);
%                 W_new = W_orig;
%             elseif(rndm_type==2)
% %                 dev_vb0 = 0;
%                 dev_vb0 = normrnd(0,sig_vb0);
%                 dev_param2 = 0;
%                 dev_param3 = 0;
%                 dev_param4 = 0;
%                 [mm,ss] = fnGetParam4Dist(W_orig, (W_orig*sigma2mean_scale)^2, Prior.LogNormal);
%                 W_new = lognrnd(mm,ss);
%             end
% 
%             L_new = L_orig;
%             fixedParams_Sel = fixedParams;
%             chngParams_Sel = chngParams_all{idx_chng_chngParams};
%             chngParams_Sel_adj = chngParams_Sel;
%             chngParams_Sel_adj(1) = chngParams_Sel(1)+dev_vb0;
%             chngParams_Sel_adj(2) = chngParams_Sel(2)+dev_param2;
%             chngParams_Sel_adj(3) = chngParams_Sel(3)+dev_param3;
%             chngParams_Sel_adj(4) = chngParams_Sel(4)+dev_param4;
% 
%                     yy_base_true = true_yys_all{idx_chng_chngParams};
% 
%                     showplot = @(Vbs,yys)plot(datanew(1:4:161,1),true_yys_all{idx_chng_chngParams},marks{idx_chng_chngParams},Vbs,yys,marks2{idx_chng_chngParams},'linewidth',2);
%                 yys_new = gen_ds_current_LW(dim_v,vg, fixedParams_Sel, chngParams_Sel_adj,constants,showplot, L_orig, W_new, verbose);
%                 obs(:,idx_chng_chngParams,repeat) = yys_new + normrnd(0,sig_err,size(yys_new));
% 
%                 if(draw_fig)
%                     figure(fig1);
%                     hold on;
%                     str_ttl = '';
%                     showplot(vg,obs(:,idx_chng_chngParams,repeat));
% 
%                     hold off;
% 
% 
% 
%                     figure(fig2);
%                     hold on;
% 
%                     plot(vg, obs(:,idx_chng_chngParams,repeat), marks2{idx_chng_chngParams});
% 
%                     hold off;
%                 end
%     %         end
%         end  
%     end
% 
%     
    figure(fig3);clf;
    for j = 1:nRepeat
        for i = 1:dim_c
            obs_ij = obs(:,i,j);
%             if(i == 1)
% %                 fit_w = sum(std_real_obs.*obs_ij)./sum(obs_ij.^2);
%                 fit_w = estimate_W(std_real_obs, obs_ij);
%             end
%             std_obs(:,i,j) = obs_ij*fit_w;
            if(draw_fig)
                figure(fig3);
                hold on;
%                 plot(vg, obs_ij*fit_w, marks2{i});
                plot(vg, obs_ij, marks2{i});
                hold off;
            end
        end
    end
% %     obs = std_obs;
% end

param.generator.mean_W = W_orig;
param.generator.sig_W = W_orig/4;
param.generator.vg = vg;
param.generator.fixedParams = fixedParams;
param.generator.betas = betas;
param.generator.constants = constants;
param.generator.sig_err = sig_err;
param.generator.sig_vb0 = sig_vb0;
% param.generator.W_all = W_all;
param.generator.std_real_obs = std_real_obs;


v_mat3 = zeros(dim_v, nCont, nSensor);
c_mat3 = zeros(dim_v, nCont, nSensor);
j_mat3 = zeros(dim_v, nCont, nSensor);
i_mat3 = zeros(dim_v, nCont, nSensor);
l_mat3 = zeros(dim_v, nCont, nSensor);
for i=1:nCont;
    for j = 1:nSensor;
        for k = 1:dim_v;
            v_mat3(k,i,j) = vg(k);
            j_mat3(k,i,j) = j;
            c_mat3(k,i,j) = cont(i);
            i_mat3(k,i,j) = i;
            l_mat3(k,i,j) = k;
        end

    end;
end;


% vg
% fixedParams
% betas
% constants
% sig_vb0
% sig_err 
% obs
% rndm_type
% sig_vb2
% sig_vb3
% sig_vb4


% %%
% sig_W = 0;
% sig_vb0 =0;
% sig_vb2 =0;
% sig_vb3 =0;
% sig_vb4 =0;
% sig_err =0;
% generator = @(contam_lv) generate_test_data_emulator_v3(contam_lv, mean_W, sig_W, vg, fixedParams,betas, constants, sig_vb0, sig_err , obs, rndm_type, sig_vb2, sig_vb3, sig_vb4)
%     figure(fig3);
% clf;
% %%
% %                 for repeat = 1:10;
% test_cont = [2:2:18]';
% test_cont = [0:10:20]';
% dim_test = max(size(test_cont));
% mark3={'ko','ro','bo'};
% for i=1:dim_test;
%     hold on;
%     figure(fig3);
%     plot(vg,generator(test_cont(i)), mark3{i});
% %     plot(vg,generator(test_cont(i)));
% end
%                 end
%     plot(vg,generator(10), mark3{i});
% 
% save(sprintf('./data_gen/%s.mat',filename));

