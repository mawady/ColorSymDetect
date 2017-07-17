function [symRes,voteMapBlur] = symBilOurCentLogGaborHSV(Image)
%% Wavelet feature extraction:
%   - Compute LogGabor response
%   - Define characteristics of edge features
%   - Declare textural and color histograms
wavParam.nAngs = 32;
wavParam.nScls = 12;
wavParam.minWaveLength = 6; 
wavParam.mult = 1.2;
wavParam.radSigma = 0.55;
wavParam.angSigma = 0.2;

wavParam.histBinNumT = 32;
wavParam.halfWindowSize = round(max(size(Image))/50);
wavParam.hopSize = 2*wavParam.halfWindowSize+1;
wavParam.magThreshold = 0.01;

wavData = computeWaveletCoeff_LogGabor(Image, wavParam);

%% Pairwise symmetry triangulation

triData = computeTriangles(wavData);

idx = [];
idx = [idx, find(isnan(triData.sym_wmp)==1)];
idx = [idx, find(isnan(triData.sym_hst)==1)];
idx = [idx, find(isnan(triData.sym_wgt)==1)];
idx = [idx, find(isnan(triData.displacement)==1)];
idx = [idx, find(isnan(triData.gamma)==1)];
idx = unique(idx);
triData.sym_wmp(idx) = [];
triData.sym_hst(idx) = [];
triData.sym_wgt(idx) = [];
triData.displacement(idx) = [];
triData.gamma(idx) = [];
triData.p(idx) = [];
triData.q(idx) = [];

%% 2D Voting projection based on symmetry weights

voteParam.accwidth = 2*ceil(sqrt(size(Image,1)^2+size(Image,2)^2))+1; % displacement
voteParam.accheight = 360; % angle
% voteParam.accwidth = 800; % displacement
% voteParam.accheight = 180; % angle

voteData = computeVotingProj(triData,voteParam);

%% Local maxima selection

maxParam.hsize = 20;
maxParam.halfwindow = 10;
maxParam.mindistbetcent = 10;
maxParam.lowerbound = 0.01;
maxParam.minarea = 0.1;

maxData = computeVotingMax(voteData,maxParam);

%% Symmetry axis computation

symData = computeSymAxis(Image,voteData,maxData,voteParam,maxParam);

%%

voteMapBlur = maxData.voteMapBlur;
symRes = zeros(size(symData.scores,1),5);
for i=1:size(symData.scores,1)
    symRes(i,1:2) = symData.axsSt{i};
    symRes(i,3:4) = symData.axsEd{i};
    symRes(i,5) = symData.scores(i);
end

if(size(symRes,1)>1)
    scr = symRes(:,5);
    [~,idx] = sort(scr,'descend');
    symRes = symRes(idx,:);
%     symRes(:,5) = symRes(:,5) ./ symRes(1,5);
end

end

