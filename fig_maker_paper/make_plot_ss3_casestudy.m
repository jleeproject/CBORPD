clc;
clear all;
% names_problems = {'};
% filenames ={'results_210722_v1_8_01_RDOs(sel)_ss2_ninit2_bx100_intermediate','results_210719_v1_8_01_RDOs(gardener2)_ss2_ninit2_bx100_intermediate'};
% filenames ={'results_210719_v1_8_01_RDOs(sel)_ss3_ninit2_bx100_intermediate','results_210718_v1_8_01_RDOs(gardener2)_ss3_ninit2_bx100_intermediate'};
% filenames ={'',''};
% filenames ={'','results_210718_v1_8_01_RDOs(gardener2)_ss5_ninit2_bx100_intermediate'};
% prob_indices={[1, 3, 4], 1};
% plot_indices={[1, 0, 3, 4], 2};



filenames ={'FINAL_results_210724_v1_8_01_Casestudy_ss3'};
% filenames ={'New Folder/results_210723_v1_8_01_Casestudy_ss2_intermediate'};
% disp(repeat);


prob_indices={1};
plot_indices={1};
create_plot_4_paper;






select_sp = @(idx) subtightplot(1,1,idx, [.06, 0.02]);



f = figure(4);

sp = select_sp(1);
sp.YLim = [0.06, 0.12];
sp.XLim = [0, 500];

f = figure(idxSumCon3);
sp = select_sp(1);
sp.YLim = [0.025, .07];
sp.XLim = [0, 500];

save_plots_in_pdf