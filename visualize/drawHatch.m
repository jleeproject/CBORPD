function drawHatch(boundaries, gridDomainD1, gridDomainD2, edgeColor, faceColor, edgeAlpha, faceAlpha, hatchColor)
    if(nargin<4)
        edgeColor = 'r';
        faceColor = 'w';
    end
    if(nargin<6)
        faceAlpha = 0.2;
        edgeAlpha = 0.5;
    end
    numberOfBoundaries = size(boundaries,1);
    for k = 1 : numberOfBoundaries
        thisBoundary = boundaries{k};
        h = fill(gridDomainD1(thisBoundary(:,2)), gridDomainD2(thisBoundary(:,1)), faceColor, 'LineWidth', 2, 'EdgeColor', edgeColor);
        h1 = hatchfill2(h,'HatchColor',hatchColor, 'HatchDensity', 30);

        hideLegend(h);
        h.FaceAlpha = faceAlpha;
        h.EdgeAlpha = edgeAlpha;
%         set(h,'FaceAlpha',alpha)
    end
end






% 
% hold off;                                 % if you want to have more control
% ax1 = gca;
% ax2 = copyobj(ax1,figure);
% 
% % Example 1: Default hatching
% hp = findobj(ax1,'Tag','HatchingRegion');
% hh = hatchfill2(hp,'cross','LineWidth',1,'Fill','off');
% title('Example 1: hatchfill2(hp,''HatchColor'',''w'',''FaceColor'',''none'')');
% 
% % Example 2: Set logarithmic yscale and reverse yaxis & speckle
% set(ax2,'ylim',[50 700],'yscale','log','ydir','reverse');
% hp = findobj(ax2,'Tag','HatchingRegion');
% % h1 = hatchfill2(hp,'speckle','HatchColor','w');
% hp = findobj(ax2,'Tag','HatchingRegion');
% 
