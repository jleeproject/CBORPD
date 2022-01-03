function mydrawpdf(fig, fontsize, filename)
%     colormap('gray');
    setFontSize(fig,fontsize);
%     print(fig, filename, '-dpdf', '-bestfit');

    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3)*1.15 fig_pos(4)];
    print(fig,'-dpdf','-painters','-r600','-bestfit',filename);

    clf();cla();
end