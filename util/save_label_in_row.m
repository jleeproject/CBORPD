function save_label_in_row(usedsubplots, linewidth, markersize, fontsize, titles, legends, filename_prefix)
%     settingfig = figure(999);ssp = subtightplot(1,1,1);;settingfig.Visible='off';
% 
%     for idx_subplot=1:numel(usedsubplots)
%         components=get(usedsubplots(idx_subplot),'children');
%         if(numel(components)>0)
% 
%             outfig = figure(2);clf;
%             copyobj(usedsubplots(idx_subplot),outfig)
%             
%             ax = findobj(outfig,'Type','Axes');
%             ax.Position = ssp.Position;
%             
%             offset = .24;
%             ax.Position(2)=ssp.Position(2) + offset
%             ax.Position(4)=ssp.Position(4) - offset
% 
%     %         figure(outfig)
% 
%             outfig.Position = [1, 50, 400 ,300];
%             % lw = findobj(ax,'Type','LineWidth');
%             set(ax.Children,'LineWidth',linewidth);
%             set(ax.Children,'MarkerSize',markersize);
%             %                         fig.PaperPositionMode = 'auto';
%             % colormap('gray');
% 
%             setFontSize(outfig,fontsize);
%             %                         setFontSize(sp2,18);
%             %                         saveas(sp, sprintf('./figs/%s_c.pdf',strrep( cell_prob_names{idx_function},':','-' ) ) , 'pdf'  );
% %             out = regexp(ax.Title.String, '(.*) \(','tokens')
%             ax.Title.String='';
% %             print(outfig, sprintf('figs/%s_%s',filename_prefix, titles{idx_subplot}), '-dpdf', '-bestfit')
%             
%             if (numel(usedsubplots(idx_subplot).Legend)>0)
%                 legend();
%                 leg = findobj(outfig,'Type','Legend')
%                 leg_comp = get(leg,'children');
%                 copyobj(legends{idx_subplot}, leg_comp);
%                 leg.NumColumns = ceil(numel(leg.String)/2);
%             %     leg.Location = 'northoutside';
% %                 leg.Location = 'northoutside';
% %                 leg.Position(1)=0.3;
% %                 leg.Position(2)=0.9;
%                 leg.Location = 'southoutside';
%                 
%                 ax.Position(2)=ssp.Position(2) + offset
%                 ax.Position(4)=ssp.Position(4) - offset
% %                 leg.Position(1)=0.3;
%                 pos = leg.Position;
%                 pos(2) = pos(2) + 0.013;
%                 leg.Position=pos;
% 
%                 leg.FontSize = 11;
% %                 leg.Position(4)=.06;
%             end
%             
% 
%             outfig.PaperPositionMode = 'auto';
%             fig_pos = outfig.PaperPosition;
%             outfig.PaperSize = [fig_pos(3) fig_pos(4)];
%             print(outfig,'-dpdf','-painters','-r600','-bestfit',sprintf('figs/%s_%s',filename_prefix, titles{idx_subplot}));
% 
%         end
%     end
end