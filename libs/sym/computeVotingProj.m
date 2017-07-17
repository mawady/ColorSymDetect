function voteData = computeVotingProj(triData,voteParam)

voteMap = zeros(voteParam.accheight,voteParam.accwidth);
countMap = zeros(voteParam.accheight,voteParam.accwidth);
pnts = cell(voteParam.accheight,voteParam.accwidth);

pntX = round(triData.gamma/pi*voteParam.accheight);
pntY = round(voteParam.accwidth*((triData.displacement+(sqrt(2)/2))/sqrt(2)));
for i=1:numel(triData.gamma)
    x = pntX(i);
    y = pntY(i);
    if(x==0)
        x = voteParam.accheight;
    end
    voteMap(x,y) = voteMap(x,y)+ ...
        triData.sym_clr(i)*triData.sym_wmp(i)*triData.sym_hst(i);
    countMap(x,y) = countMap(x,y)+1;
    if(~isempty(pnts{x,y}))
        pnts{x,y} = [pnts{x,y},{triData.p{i}},{triData.q{i}}];
    else
        pnts{x,y} = [{triData.p{i}},{triData.q{i}}];
    end
end
voteData.voteMap = voteMap;
% voteData.voteMap = voteMap./sum(sum(voteMap));
voteData.countMap = countMap;
voteData.pnts = pnts;

end

