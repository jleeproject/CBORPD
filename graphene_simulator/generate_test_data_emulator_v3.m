function [test_obs] = generate_test_data_emulator_v3(contam_lv, mean_W, sig_W, vg, fixedParams,betas, constants, sig_vb0, sig_err , obs, rndm_type, sig_vb2, sig_vb3, sig_vb4)
    if(nargin==1 || nargin==2)
        param = contam_lv;
        if( nargin==2)
            obs = mean_W;
        end
      
        contam_lv = param.cont_new;
        mean_W = param.mean_W;
        sig_W = param.sig_W;
        vg = param.vg;
        fixedParams = param.fixedParams;
        betas = param.betas;
        constants = param.constants;
        sig_err = param.sig_err;
        sig_vb0 = param.sig_vb0;
        std_real_obs = param.std_real_obs;
    end
%     test_obs = zeros(dim_v,1);
% % ======================== MAKE TEST SET ======================
%     test_ind_a = a_line_TRUE(1)+ a_line_TRUE(2)*test_contm + normrnd(0,sig_aa);
%     test_ind_b = b_line_TRUE(1)+ b_line_TRUE(2)*test_contm + normrnd(0,sig_bb);
%     test_ind_c = c_line_TRUE(1)+ c_line_TRUE(2)*test_contm + normrnd(0,sig_cc);
% %     for k = 1:dim_v;
% %         test_obs(k) = test_ind_a .*( xx(k)-test_ind_b).^2 + test_ind_c + normrnd(0,sigerr);
% %     end
%     test_obs = test_ind_a .*( xx-test_ind_b).^2 + test_ind_c + normrnd(0,sigerr, dim_v, 1);
    params_gen = zeros(1,4);
    params_gen_0 = zeros(1,4);
    if(size(betas,1)==2)
        for idx = 1:4;
            params_gen(:,idx) = betas(1,idx) + betas(2,idx)*contam_lv;
            params_gen_0(:,idx) = betas(1,idx) ;
        end
    elseif(size(betas,1)==3)
        for idx = 1:4;
            params_gen(:,idx) = betas(1,idx) + betas(2,idx)*contam_lv + betas(3,idx)*contam_lv.^2;
            params_gen_0(:,idx) = betas(1,idx) ;
        end
    end
        chngParams_Sel = params_gen(1 ,:);
        chngParams_Sel_adj = chngParams_Sel;
%%
            if(rndm_type==1)
                dev_vb0 = normrnd(0,sig_vb0);
                dev_param2 = normrnd(0,sig_vb2);
                dev_param3 = normrnd(0,sig_vb3);
                dev_param4 = normrnd(0,sig_vb4);
                W_new = mean_W;
            elseif(rndm_type==2)
%                 dev_vb0 = 0;
                dev_vb0 = normrnd(0,sig_vb0);
                dev_param2 = 0;
                dev_param3 = 0;
                dev_param4 = 0;
                [mm,ss] = fnGetParam4Dist(mean_W, (sig_W)^2, Prior.LogNormal);
                W_new = lognrnd(mm,ss);
            elseif(rndm_type==3)
                dev_vb0 = normrnd(0,sig_vb0);
                dev_param2 = normrnd(0,sig_vb2);
                dev_param3 = normrnd(0,sig_vb3);
                dev_param4 = normrnd(0,sig_vb4);
                [mm,ss] = fnGetParam4Dist(mean_W, (sig_W)^2, Prior.LogNormal);
                W_new = lognrnd(mm,ss);
            end

%             L_new = L_orig;
%             fixedParams_Sel = fixedParams;
%             chngParams_Sel = chngParams_all{idx_chng_chngParams};
            chngParams_Sel_adj = chngParams_Sel;
            chngParams_Sel_adj(1) = chngParams_Sel(1)+dev_vb0;
            chngParams_Sel_adj(2) = chngParams_Sel(2)+dev_param2;
            chngParams_Sel_adj(3) = chngParams_Sel(3)+dev_param3;
            chngParams_Sel_adj(4) = chngParams_Sel(4)+dev_param4;

%%
%         dev_vb0 = normrnd(0,sig_vb0);
%         chngParams_Sel_adj(1) = chngParams_Sel(1)+dev_vb0;
%         [mm,ss] = fnGetParam4Dist(mean_W, (sig_W)^2, Prior.LogNormal);
%         W_new = lognrnd(mm,ss);

        
        nVbs = max(size(vg));
        yys_new = zeros(nVbs,1);
        verbose = false;
        for i=1:nVbs
            yys_new(i) = fn_solve_eq_graphene_LW(vg(i),fixedParams,chngParams_Sel_adj, constants, constants.L, W_new, verbose);
        end
        
%         chngParams_Sel = params_gen_0(1 ,:); %% When c=0
%         chngParams_Sel_adj = chngParams_Sel;
%         dev_vb0 = normrnd(0,sig_vb0);
%         chngParams_Sel_adj(1) = chngParams_Sel(1)+dev_vb0;

%         yys_new_baseline = zeros(nVbs,1);
%         for i=1:nVbs
%             yys_new_baseline(i) = fn_solve_eq_graphene_LW(vg(i),fixedParams,chngParams_Sel_adj, constants, constants.L, W_new, verbose);
%         end
        
        %       yys_new = fn_calc_draw_plots_LW(dim_v,vg, fixedParams, chngParams_Sel_adj,constants,showplot, L_new, W_new, verbose);title(str_ttl);
%         test_baseline = yys_new_baseline + normrnd(0,sig_err,size(yys_new));
%         test_baseline = obs(:,1,1);
        test_obs = yys_new + normrnd(0,sig_err,size(yys_new));
%         fit_w = estimate_W(obs(:,1,1), test_baseline);
%         std_test_obs = test_obs*fit_w;
%         figure(17);clf;plot(vg,reshape(obs,nVbs,[]),'-k', vg,std_test_obs,'y-', vg, test_baseline*fit_w, '-r');
    end