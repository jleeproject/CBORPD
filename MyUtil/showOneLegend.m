function showOneLegend(fig_handler, str_legend)
    if(numel(fig_handler)>0)
        fig_handler(1).Annotation.LegendInformation.IconDisplayStyle = 'on';
        fig_handler(1).DisplayName = str_legend;
        for i=2:numel(fig_handler)
            fig_handler(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end
end