function savePng(gcf, filename,position)
    if(exist('position','var') )
        if(isnan(position))
            set(gcf, 'Units', 'inch', 'position', position);
        else
            set(gcf, 'Units', 'inch', 'position', [0, 0,5, 5]);
        end
    else
        set(gcf, 'Units', 'inch', 'position', [0, 0,5, 5]);
%         set(gcf, 'Units', 'inch', 'position', [-6.5000,3.8750,5, 5]);
    end
%     set(gcf,'Resize','off')
    set(gca, 'LooseInset', get(gca,'TightInset'))

    saveas(gcf,sprintf('%s.fig',filename),'fig');
    saveas(gcf,sprintf('%s.png',filename),'png');
end

