function f1=setAxis(f1, xlim, ylim)
    if(nargin==1)
        fontSize=15;
    end
    if(max(size(f1.Children))>1)
        nn = max(size(f1.Children));
        for ii= 1: nn;
            if(isa(f1.Children(ii),'matlab.graphics.axis.Axes'))
%                 f1.Children(ii).FontSize=fontSize;
                if(min(size(xlim))>0)
                    f1.Children(ii).XLim = xlim;
                end
                if(min(size(ylim))>0)
                    f1.Children(ii).YLim = ylim;
                end
            end
        end
    else
                if(min(size(xlim))>0)
                    f1.Children.XLim = xlim;
                end
                if(min(size(ylim))>0)
                    f1.Children.YLim = ylim;
                end
    end

end