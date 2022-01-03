function f1=setFontSize(f1, fontSize)
    if(nargin==1)
        fontSize=15;
    end
    if(max(size(f1.Children))>1)
        nn = max(size(f1.Children));
        for ii= 1: nn;
            if(isa(f1.Children(ii),'matlab.graphics.axis.Axes'))
                f1.Children(ii).FontSize=fontSize;
            end
        end
    else
        f1.Children.FontSize=fontSize;
    end

end