clear all;
type = 2 ;


% names_problems = {'};
filenames ={'save_recent3/Final/Penalty/results_210802_v1_8_01_RDOs(sel)_ss2_ninit2_bx100_penalty10_intermediate',''};
prob_indices={[1, 2, 3], 1};
plot_indices={[1, 2, 3], 2};



% filenames ={'results_210723_v1_8_01_Casestudy_ss2_intermediate'};
% filenames ={'New Folder/results_210723_v1_8_01_Casestudy_ss2_intermediate'};
% prob_indices={1};
% plot_indices={1};
if type==2
    create_plot_4_paper_diff_arr
    select_sp = @(idx) subtightplot(2,2,idx, [.06, 0.02]);
else
    create_plot_4_paper;
    select_sp = @(idx) subtightplot(4,1,idx, [.06, 0.02]);
end

create_plot_4_paper;


figure(4);

for i=1:4
    sp = select_sp(i);
    disp(sp.YLim);
end

sp = select_sp(1);
sp.YLim = [-.01, 2.2];

sp = select_sp(2);
sp.YLim = [1.25, 2];

sp = select_sp(3);
sp.YLim = [0, 20];

% sp = select_sp(4);
% % sp.YLim = [3.5, 4.3];
% sp.YLim = [3.5, 10];



figure(idxSumCon3);

for i=1:4
    sp = select_sp(i);
    disp(sp.YLim);
end

sp = select_sp(1);
sp.YLim = [-.01, 1.5];

sp = select_sp(2);
sp.YLim = [-.01, 100];
sp.XLim = [0,200]
sp = select_sp(3);
sp.YLim = [-.01, 250];

sp = select_sp(4);
% sp.YLim = [3.5, 4.3];
sp.YLim = [-.01, 30];


save_plots_in_pdf