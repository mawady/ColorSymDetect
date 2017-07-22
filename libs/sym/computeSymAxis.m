function symData = computeSymAxis(Image,voteData,maxData,voteParam,maxParam)

axisLoc = cell(size(maxData.locs,2),1);
axisPnts = cell(size(maxData.locs,2),1);
scores = zeros(size(maxData.locs,2),1);
wd = maxParam.halfwindow;
for i=1:size(maxData.locs,2)
    [row,col] = getSegment(Image,maxData.locs(:,i),voteParam);
    axisLoc{i} = [row; col]';
    axisRows = maxData.locs(1,i)-wd:...
        maxData.locs(1,i)+wd;
    axisRows(axisRows<1) = 1;
    axisRows(axisRows>size(voteData.pnts,1)) = size(voteData.pnts,1);
    axisRows = unique(axisRows);
    axisCols = maxData.locs(2,i)-wd:...
        maxData.locs(2,i)+wd;
    axisCols(axisCols<1) = 1;
    axisCols(axisCols>size(voteData.pnts,2)) = size(voteData.pnts,2);
    axisCols = unique(axisCols);
    axisPnts{i} = voteData.pnts(axisRows,axisCols);
%     scores(i) = voteData.voteMap(maxData.locs(1,i),maxData.locs(2,i));
    scores(i) = maxData.voteMapBlur(maxData.locs(1,i),maxData.locs(2,i));
end

axsSt = {};
axsEd = {};
for i=1:size(maxData.locs,2)
    bilSymC_Our = axisLoc{i};
    if(~isempty(bilSymC_Our))
        [stPntC_Our,edPntC_Our] = getEndPointsLine(bilSymC_Our(:,1),...
            bilSymC_Our(:,2));
        pntsC_Our = axisPnts{i};
        pntsC_Our = pntsC_Our(:);
        tmpVotePnts = [];
        for k=1:numel(pntsC_Our)
            tmpVotePnts = [tmpVotePnts,pntsC_Our{k}];
        end
        pntsC_Our = cell2mat(tmpVotePnts)';
        pntsC_Our = unique(pntsC_Our,'rows');
        if(size(pntsC_Our,1)>2)
            pntsC_Our(:,1) = pntsC_Our(:,1)*max(size(Image))+size(Image,1)/2;
            pntsC_Our(:,2) = pntsC_Our(:,2)*max(size(Image))+size(Image,2)/2;
            try
            [stPntC_Our,edPntC_Our] = getConvexHullInter(pntsC_Our,...
                stPntC_Our,edPntC_Our);
            catch
            end
        end
        axsSt{i} = fliplr(stPntC_Our);
        axsEd{i} = fliplr(edPntC_Our);
    end
end

symData.axisLoc = axisLoc;
symData.axisPnts = axisPnts;
symData.scores = scores;
symData.axsSt = axsSt;
symData.axsEd = axsEd;

end

function [row,col] = getSegment(Image,loc,voteParam)

ro = (-sqrt(2)/2) + (0:(voteParam.accwidth-1)) ...
    * sqrt(2)/(voteParam.accwidth-1);
ag = (pi*(0:(voteParam.accheight-1))/(voteParam.accheight-1));
X_d = ro(loc(2));
X_g = ag(loc(1));

x = -round(sqrt(size(Image,1)^2+size(Image,2)^2)):...
    round(sqrt(size(Image,1)^2+size(Image,2)^2));
y = (X_d-x*cos(X_g))/sin(X_g);
xx = x*max(size(Image))+(size(Image,1)/2);
yy = y*max(size(Image))+(size(Image,2)/2);

y1 = -round(sqrt(size(Image,1)^2+size(Image,2)^2)):...
    round(sqrt(size(Image,1)^2+size(Image,2)^2));
x1 = (X_d-y1*sin(X_g))/cos(X_g);
xx1 = x1*max(size(Image))+(size(Image,1)/2);
yy1 = y1*max(size(Image))+(size(Image,2)/2);

if(var(xx)<var(xx1))
    row = xx;
    col = yy;
else
    row = xx1;
    col = yy1;
end

end

