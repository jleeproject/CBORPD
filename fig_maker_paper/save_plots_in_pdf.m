

fontsize = 13;
linewidth = 2;
markersize = 4;
% haveLegend = true;
% infig=figure(idxSummary);
settingfig = figure(3);ssp = subtightplot(1,1,1);;settingfig.Visible='off';


infig=figure(idxSumCon3);

usedsubplots = findobj(infig, 'Type','axes');

titles = cell(numel(usedsubplots),1);
for idx_subplot=1:numel(usedsubplots)
    out = regexp(usedsubplots(idx_subplot).Title.String, '(.*) \(','tokens');
    if numel(out)>0
        titles{idx_subplot} = out{1}{1};
    end
end


legends = cell(numel(usedsubplots),1);
for idx_subplot=1:numel(usedsubplots)
    out= usedsubplots(idx_subplot).Legend;
%     out = regexp(usedsubplots(idx_subplot).Title.String, '(.*) \(','tokens');
    if numel(out)>0
        legends{idx_subplot} = out;
    end
end

printSubplotsSeparate(usedsubplots, linewidth, markersize, fontsize, titles, legends, sprintf('d%d_fg',setting.samplesize_for_variance_init))


infig=figure(idxSummary);
usedsubplots = findobj(infig, 'Type','axes');
printSubplotsSeparate(usedsubplots, linewidth, markersize, fontsize, titles, legends, sprintf('d%d_fo',setting.samplesize_for_variance_init))
