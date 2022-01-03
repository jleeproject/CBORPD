function save_cropped_pdf(f1,path)
    figure(f1);
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters epsFig
%     saveas(f1,'./figs/f1_s_all.pdf');
    saveas(f1,path);
end