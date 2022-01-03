function pd1 = fitNormalPlot(data, strLegend, strTitle, mu, std)
hold on;
LegHandles = []; LegText = {};
% --- Plot data originally in dataset "data data"
[CdfF,CdfX] = ecdf(data,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 1;
[~,BinEdge] = internal.stats.histbins(data,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
% hLine = bar(BinCenter,BinHeight,'hist');
hLine = histogram(data,'Normalization','pdf'); 
% title(strTitle);
% set(hLine,'FaceColor','none','EdgeColor',[0.333333 0 0.666667],...
set(hLine,'EdgeColor',[0.333333 0 0.666667],...
    'LineStyle','-', 'LineWidth',1);
% xlabel('Data');
% ylabel('Density')
% LegHandles(end+1) = hLine;
% LegText{end+1} = strLegend;

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);
if(exist('mu','var'))
    YPlot = normpdf(XGrid,mu,std);
else
    pd1 = fitdist(data, 'normal');
    YPlot = pdf(pd1,XGrid);
end
hLine = plot(XGrid,YPlot,'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
% LegHandles(end+1) = hLine;
% LegText{end+1} = 'normal';
hold off;  hLine.Parent.FontSize=13;
% hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northeast');
% set(hLegend,'Interpreter','none');

