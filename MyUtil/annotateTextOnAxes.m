function annotateTextOnAxes(fig, nRow, nCol, strRows, strCols, strTitle)
    if(numel(fig)>0)
    %     fig4=figure(4);clf;
    %     nRow = 4;
    %     nCol = 3;
    %     strRows = {'r1','r2','r3','r4'};
    %     strCols = {'c1','c2','c3'};
    %     strTitle = 'ABC';

        left = 0.00;
        top = 0.01;
        width_row = 0.1;
        width_col = 0.21;
        width_title = 0.8;
        height = 0.06;
        left_margin = 0.085;
        top_margin = 0.05;
        adjust_vertical = 0.1;
        adjust_horizontal = 0.05;


        annotateText(fig,left,top,width_title,height,strTitle);

        delta = (1-2*top_margin)/(nRow);
        for i=1:nRow
            sel_top = (i-1)*delta +  top_margin + adjust_vertical;
            annotateText(fig,left,sel_top,width_row,height,strRows{i}, 'HorizontalAlignment','center');
        end
        delta = (1-2*left_margin)/nCol;
        for i=1:nCol
            sel_left = (i-1)*delta + left_margin + adjust_horizontal;
            annotateText(fig,sel_left,top,width_col,height,strCols{i}, 'HorizontalAlignment','center');
        end

    %     for i=1:nRow*nCol
    %         subplot(nCol,nRow,i);
    %     end
    end
end