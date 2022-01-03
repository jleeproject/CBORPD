function pd1 = fitDistributionPlot(distname, data, strLegend, strTitle, params)
hold on;
LegHandles = []; LegText = {};
% --- Plot data originally in dataset "data data"
[CdfF,CdfX] = ecdf(data,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 1;
[~,BinEdge] = internal.stats.histbins(data,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
% hLine = bar(BinCenter,BinHeight,'hist');
hLine = histogram(data,'Normalization','pdf'); %title(strTitle);
% set(hLine,'FaceColor','none','EdgeColor',[0.333333 0 0.666667],...
set(hLine,'EdgeColor',[0.333333 0 0.666667],...
    'LineStyle','-', 'LineWidth',1);
% xlabel('Qm');
% ylabel('Density')
LegHandles(end+1) = hLine;
% LegText{end+1} = strLegend;

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);
if(strcmp(distname,'chi2'))
    if(exist('params','var'))
        df = params;
    else
        df = mean(data,'all');
    end
    YPlot = pdf('chi2',XGrid,df);
    LegText{end+1} = sprintf('Chi^2(%d)',df);
    temp_axis = axis; axis([0 max(data) temp_axis(3:4)]);

else
    pd1 = fitdist(data, distname);
    YPlot = pdf(pd1,XGrid);
    LegText{end+1} = distname;
end
hLine = plot(XGrid,YPlot,'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
LegHandles(end+1) = hLine;
hold off;hLine.Parent.FontSize=14;hLine.Parent.XAxis.Limits=[0 , 15];hLine.Parent.Position = [0.1300    0.1285    0.7750    0.7965];
% hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northeast');
% set(hLegend,'Interpreter','none');

