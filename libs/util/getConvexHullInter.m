function [stPnt, edPnt] = getConvexHullInter(pnts,stPntIN,edPntIN)
    
chInd = convhull(pnts(:,1),pnts(:,2));
ch = pnts(chInd,:);
[xi,yi] = intersections(ch(:,1),ch(:,2),...
    [stPntIN(1),edPntIN(1)],[stPntIN(2),edPntIN(2)]);
%[xi,yi] = polyxpoly(ch(:,1),ch(:,2),...
%    [stPntIN(1),edPntIN(1)],[stPntIN(2),edPntIN(2)]);
if(~isempty(xi) && numel(xi)>1)
    stPnt = [xi(1),yi(1)];
    edPnt = [xi(2),yi(2)];
else
    stPnt = stPntIN;
    edPnt = edPntIN;
end