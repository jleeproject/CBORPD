clear;clc;
% clearProblems;
% % setting =struct();
% setting.typeProblem = TypeProblem.RobustDesignOptimization;
% % -----------------------------------------------------------------------
% % import_problem_RDO_gramacy_easier_43;     addProblems;
% import_problem_RDO_gardener14_1_medium;    addProblems;
% % import_problem_RDO_gardener14_2_41;        addProblems;
% % import_problem_RDO_gelbart14_hard;       addProblems;
% % import_problem_RDO_ariafar19_easier_44;     addProblems;
% % 
% % 
% 
% obj = cell_obj_func{1};
% dom = obj.xDomain;
% numBin = 100;
% numBin2 = numBin;
% [xx,yy, ranges, ranges2]= make2dRangesForNumericalStudy(dom(1,1),dom(1,2), numBin, dom(2,1),dom(2,2), numBin2);
% 
% fcon = cell_type_mean_func{1}{1};
% 
% con = cell_const_func{1}.constraints{1};
% 
% f = figure(55);clf;
% flippedImagesc(ranges, ranges2, obj.fnEval(xx,yy));

resume = false;
% resume = true;
if(resume); if (~askYesNo('Are you sure you want RESUME?'));return;end;end
clc;fprintf('BO starts at %s\n', showPrettyDateTime(now()));setVersion; fprintf('%s\n',versionBO);if(~resume);clear all; resume=false;end

% scr_run_on_server
scr_run_on_my_pc
% scr_setting_0_preprocessing;  
% scr_RDO_selected;
% scr_RDO_easier_gramacy
% scr_RDO_easier_gardener2
% scr_RDO_easier_ariafar
scr_RDO_selected2;

% scr_RDO_easier;
% scr_RDO_easy;
% scr_RDO;
% scr_RDO_medium;
% scr_RDO_hard;
% scr_DET;
% scr_STO;
% scr_STO_letham;
% scr_STO_letham_hard
% scr_STO_letham_medium
% scr_Casestudy;
% scr_common;
% scr_common_sel2;
scr_test_new;

create_pdf = true;
bin = 20;
% textstep =0;
textstep = 1;
fontsize = 35;
markerwidth = 20;
dec_roundup = 1;

if (create_pdf)
%     mysubplot = @(m,n,p) figure(p);
    mysubplot = @(m,n,p) subtightplot(1,1,1);
else
    mysubplot = @(m,n,p) subtightplot(m,n,p);
end

stt_all = tic();cnt_exp = 1;; time_exp_start = now;
elapsedTime_settings = zeros(nTypeSettings,1);
tot_exp = numel(cell_obj_func)* nRepeats* nTypeSettings* numel(type_cell_samplesize);
stt_periodic = tic();
fig = figure(55);clf;
% for repeat = 1:nRepeats;
    stt_repeat = tic();
    for idx_function = 1:numel(cell_obj_func)
        scr_setting_1_in_functions_sameObj;

        stt_function = tic();
        warning('off');
        if(~resume);scr_setting_2_in_repeats_same_initPts;end;
%         for idx_setting = 1:nTypeSettings
%             stt_setting = tic;
%             scr_setting_3_in_settings; 
%     
%             for idx_samplesize = 1:numel(type_cell_samplesize)
%                 fprintf('[Exp] func:%d/%d, rep:%d/%d, set:%d/%d, samp:%d/%d (%d/%d)\n', ...
%                     idx_function, numel(cell_obj_func), repeat, nRepeats, idx_setting, nTypeSettings,...
%                     idx_samplesize, numel(type_cell_samplesize),...
%                     cnt_exp, tot_exp);
                
%                 scr_setting_4_in_samples;
                
                scr_set_settings_pre_run;
                %%
                if setting.typeProblem == TypeProblem.RobustDesignOptimization
                    fcon = simulator.obj.meanFunc;
                    fobj = simulator.obj.sigmaFunc;
                dom = fobj.xDomain;
                else
                    fcon = simul.con{1}.meanFunc;                    
                    fobj = optProb.objFunc;
                    dom = optProb.xDomain;
                end
                numBin = 100;
                numBin2 = numBin;
                [xx,yy, ranges, ranges2]= make2dRangesForNumericalStudy(dom(1,1),dom(1,2), numBin, dom(2,1),dom(2,2), numBin2);

%                 fcon = cell_type_mean_func{1}{1};

%                 con = cell_const_func{1}.constraints{1};
                
                xeval = reshape(xx,[],1);
                yeval = reshape(yy,[],1);
                if(size(dom,1)>2)
                    x3 = .1* ones(size(xeval));
                    x4 = .1* ones(size(xeval));
                    x5 = .1* ones(size(xeval));
                    x6 = .1* ones(size(xeval));
                    geval = reshape( fcon.fnEval(xeval,yeval, x3, x4, x5, x6), size(xx,1),[]);
                    feval = reshape( log(fobj.fnEval(xeval,yeval, x3, x4, x5, x6)), size(xx,1),[]);
                    
                    
                    v = @(x) reshape(x,[],1);
                    b = @(x,X) reshape(x,size(X,1),size(X,2));
                    [X,Y] = meshgrid(ranges,ranges2);
                    X1 = v(X);
                    X2 = v(Y);
                    X3 = v( .1* ones(size(X)) );
                    X4 = v( .1* ones(size(X)) );
                    X5 = v( .1* ones(size(X)) );
                    X6 = v( .1* ones(size(X)) );
                    geval_mesh = b(fcon.fnEval(X1,X2,X3,X4,X5,X6), X);
%                     feval = reshape( log(fobj.fnEval(xeval,yeval)), size(xx,1),[]);
                    feval_mesh = b(fobj.fnEval(X1,X2,X3,X4,X5,X6), X);
                else
                    geval = reshape( fcon.fnEval(xeval,yeval), size(xx,1),[]);
%                     feval = reshape( log(fobj.fnEval(xeval,yeval)), size(xx,1),[]);
                    feval = reshape( (fobj.fnEval(xeval,yeval)), size(xx,1),[]);
                    

%                     x = -2:0.2:2;
%                     y = -2:0.2:3;
                    [X,Y] = meshgrid(ranges,ranges2);
                    geval_mesh = fcon.fnEval(X,Y);
%                     feval = reshape( log(fobj.fnEval(xeval,yeval)), size(xx,1),[]);
                    feval_mesh = fobj.fnEval(X,Y);
                    
                end
                    con = cell_const_func{idx_function}.constraints{1};
                
                    if con.hasUb && ~con.hasLb
                        feasible = geval < con.ub;
                    elseif ~con.hasUb && con.hasLb
                        feasible = geval > con.lb;
                    elseif con.hasUb && con.hasLb
                        feasible = con.lb < geval && geval < con.ub;
                    else
                        disp('error');
                        return;
                    end
                    vec_feval = reshape(feval,[],1);
                    feasobj = vec_feval(reshape(feasible,[],1));
                    [vmin, imin] = min(feasobj);
                    
                    x2 = xeval(feasible);
                    y2 = yeval(feasible);
                
                    
                    sp = mysubplot(2, numel(cell_obj_func), idx_function);
%                     subplot(numel(cell_obj_func), 2, (idx_function-1)*2+ 1);
%                     flippedImagesc(ranges, ranges2, feval);
                    %%
                    [C,h] = contour(X,Y,feval_mesh,bin,'ShowText','on', 'TextStep',textstep);
%                     [C,h] = contour(X,Y,feval_mesh,'ShowText','on');
%                     [C,h] = contour(X,Y,feval_mesh,bin, 'ShowText','on');
                    h.LevelList=round(h.LevelList,dec_roundup)
                    clabel(C,h,'FontSize',fontsize)

                    hold on;
                    plot(x2(imin),y2(imin),'rx', 'LineWidth',markerwidth)
%                     colorbar()
                    if(create_pdf)
                        hold off;
                        fig.Position = [1, 50, 400 ,400];
                        setFontSize(fig,fontsize);
%                         saveas(sp, sprintf('./figs/%s_o.pdf',strrep( cell_prob_names{idx_function},':','-' ) ) , 'pdf'  );
                        print(fig, sprintf('./figs/%s_o.pdf',regexprep( cell_prob_names{idx_function},'\s?\(.*\)\s?','' ) ), '-dpdf', '-bestfit')
                        clf();cla();
                    else
                        title(cell_prob_names{idx_function})
                    end

                    sp = mysubplot(2, numel(cell_obj_func), numel(cell_obj_func) + idx_function);
%                     subplot(numel(cell_obj_func), 2, (idx_function-1)*2+ 2);

%                     flippedImagesc(ranges, ranges2, geval);
                    %%
                    [C,h] = contour(X,Y,geval_mesh,bin,'ShowText','on', 'TextStep',textstep);
%                     [C,h] = contour(X,Y,geval_mesh,'ShowText','on');
% %                     [C,h] = contour(X,Y,geval_mesh, bin, 'ShowText','on');
                    h.LevelList=round(h.LevelList,dec_roundup)
                    clabel(C,h,'FontSize',fontsize)

                    hold on;
                    boundaries = bwboundaries(feasible);
                    drawBoundaries(boundaries, ranges, ranges2, 'm', 'k', 0.3, 0.3)
%                     colorbar()
                    disp('');
                    
                    plot(x2(imin),y2(imin),'rx', 'LineWidth',markerwidth)
                    
                    if(create_pdf)
                        hold off;
                        fig.Position = [1, 50, 400 ,400];
%                         fig.PaperPositionMode = 'auto';

                        setFontSize(fig,fontsize);
%                         setFontSize(sp2,18);
%                         saveas(sp, sprintf('./figs/%s_c.pdf',strrep( cell_prob_names{idx_function},':','-' ) ) , 'pdf'  );
                        print(fig, sprintf('./figs/%s_c.pdf',regexprep( cell_prob_names{idx_function},'\s?\(.*\)\s?','' ) ), '-dpdf', '-bestfit')
                        clf();cla();
%                     else
%                         title(cell_prob_names{idx_function})
                    end
                
%                 return;
%             end
%         end
    end
% end