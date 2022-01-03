function hideLegend(fig_handler, str_legend_no_use)
    for i=1:numel(fig_handler)
        fig_handler(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end

end