function drawBoundaries(boundaries, gridDomainD1, gridDomainD2, edgeColor, faceColor, edgeAlpha, faceAlpha)
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
        hideLegend(h);
        h.FaceAlpha = faceAlpha;
        h.EdgeAlpha = edgeAlpha;
%         set(h,'FaceAlpha',alpha)
    end
end