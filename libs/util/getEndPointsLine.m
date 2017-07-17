function [stPnt,edPnt] = getEndPointsLine(X,Y)

[X,idx] = sort(X);
Y = Y(idx);
coeffsVer = polyfit(X,Y, 1);
fittedxVer = [min(X), max(X)];
fittedyVer = polyval(coeffsVer, fittedxVer);
stPntVer = [fittedxVer(1),fittedyVer(1)];
edPntVer = [fittedxVer(2),fittedyVer(2)];

[Y,idx] = sort(Y);
X = X(idx);
coeffsHor = polyfit(Y,X, 1);
fittedyHor = [min(Y), max(Y)];
fittedxHor = polyval(coeffsHor, fittedyHor);
stPntHor = [fittedxHor(1),fittedyHor(1)];
edPntHor = [fittedxHor(2),fittedyHor(2)];

if(sum(isnan(stPntHor))>0 || sum(isnan(edPntHor))>0)
    stPnt = stPntVer;
    edPnt = edPntVer;
    return;
end
if(sum(isnan(stPntVer))>0 || sum(isnan(edPntVer))>0)
    stPnt = stPntHor;
    edPnt = edPntHor;
    return;
end
if(pdist([stPntVer;edPntVer],'euclidean') > ...
        pdist([stPntHor;edPntHor],'euclidean'))
    stPnt = stPntVer;
    edPnt = edPntVer;
else
    stPnt = stPntHor;
    edPnt = edPntHor;
end