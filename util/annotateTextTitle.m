function annotateTextTitle(fig, strTitle)
    if(numel(fig)>0)

        left = 0.00;
        top = 0.01;
        width_title = 0.8;
        height = 0.06;



        annotateText(fig,left,top,width_title,height,strTitle);

    end
end