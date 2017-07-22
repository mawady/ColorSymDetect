function wavData = computeWaveletCoeff_LogGabor( img, wavParam)
%%
nAngs = wavParam.nAngs;
nScls = wavParam.nScls;
magThreshold = wavParam.magThreshold;
halfWindowSize = wavParam.halfWindowSize;
hopSize = wavParam.hopSize;
histBinNumT = wavParam.histBinNumT;
%%
binEdgesT = linspace(1, nAngs, histBinNumT);
%% Image pre-processing: 
%   - RGB-to-GrayScale 
%   - Rescale to [0,1]
imgRGB = img;
if(size(imgRGB,3)>1)
    img = rgb2gray(imgRGB);
end
img = double(img)/255;
%%
EO = computeLogGaborResponse(img, wavParam);
%%
imRs = zeros(size(img,1),size(img,2),nScls,nAngs);
for i = 1:nScls
    for j = 1:nAngs
        gaborAbs = EO{i,j};
        R = real(gaborAbs);
        Z = imag(gaborAbs);
        rsp = sqrt(R.*R+Z.*Z);
        imRs(:,:,i,j) = rsp;
    end
end

[M,IM] = max(imRs,[],4);
[ampM,sclM] = max(M,[],3);
angM = zeros(size(sclM));
for i=1:size(sclM,1)
    for j=1:size(sclM,2)
        angM(i,j) = IM(i,j,sclM(i,j));
    end
end
%%
nr = min([floor((size(img,1)-2*halfWindowSize)/hopSize)+1 size(img,1)]);
nc = min([floor((size(img,2)-2*halfWindowSize)/hopSize)+1 size(img,2)]);
C = zeros(nr,nc);
A = zeros(nr,nc);
X = zeros(nr,nc);
Y = zeros(nr,nc);
V = zeros(nr,nc,histBinNumT);
VI = zeros(nr,nc);
SC = zeros(nr,nc);
rows = round(linspace(halfWindowSize+1,size(img,1)-halfWindowSize,nr));
cols = round(linspace(halfWindowSize+1,size(img,2)-halfWindowSize,nc));
for k = 1:nr
    for l = 1:nc
        row = rows(k);
        col = cols(l);
        if(row <= halfWindowSize || col <= halfWindowSize )
            continue;
        end
        ampMs = ampM(row-halfWindowSize:row+halfWindowSize,col-halfWindowSize:col+halfWindowSize);
        angMs = angM(row-halfWindowSize:row+halfWindowSize,col-halfWindowSize:col+halfWindowSize);
        sclMs = sclM(row-halfWindowSize:row+halfWindowSize,col-halfWindowSize:col+halfWindowSize);
        C(k,l) = max(max(ampMs));
        [rm,cm] = find(ampMs == C(k,l));

%         A(k,l) = (angMs(rm,cm)-1)*pi/nAngs+pi/2;
        A(k,l) = (angMs(rm,cm)-1)*pi/nAngs-pi/2;
        X(k,l) = row+rm-(halfWindowSize+1);
        Y(k,l) = col+cm-(halfWindowSize+1);
        SC(k,l) = sclMs(rm,cm);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %--Calculation of histogram of magnitude orientations
        orVec = angMs(:);
        mgVec = ampMs(:);
                
        [~,~,whichBin] = histcounts(orVec,binEdgesT);
        whichBin = whichBin+1;
        histOrMag = zeros(1,histBinNumT);
        for z=1:numel(whichBin)
            histOrMag(whichBin(z)) = histOrMag(whichBin(z)) + mgVec(z);
        end
        indMaxOr = whichBin(find(mgVec == max(mgVec),1));
        histOrMag = histOrMag ./ sum(histOrMag);
        V(k,l,:) = histOrMag;
        VI(k,l) = indMaxOr;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%
M = zeros(size(C));
t = magThreshold*max(max(C));
for k = 2:nr-1
    for l = 2:nc-1
        bM = C(k-1:k+1,l-1:l+1);
        if C(k,l) > 0.5*max(max(bM)) && min(min(bM)) > t
            M(k,l) = C(k,l);
        end
    end
end
%%
s = sum(sum(M > 0));
m = zeros(1,s);
a = zeros(1,s);
x = zeros(1,s);
y = zeros(1,s);
sc = zeros(1,s);
v = cell(1,s);
vi = zeros(1,s);
[nr,nc] = size(M);
index = 0;
for i = 1:nr
    for j = 1:nc
        if M(i,j) > 0
            index = index+1;
            m(index) = M(i,j);
            a(index) = A(i,j);
            x(index) = X(i,j);
            y(index) = Y(i,j);
            sc(index) = SC(i,j);
            v{index} = V(i,j,:);
            vi(index) = VI(i,j);
        end
    end
end
%%
hC = cell(1,s);
for i=1:s
    rngR = x(i)-halfWindowSize:x(i)+halfWindowSize;
    rngC = y(i)-halfWindowSize:y(i)+halfWindowSize;
    if(size(imgRGB,3)==1)
        imP = im2double(imgRGB(rngR,rngC));
        imP = (imP - mean(imP(:)))/std(imP(:));
        imP(isnan(imP)) = 0;
        RGBq = 32;
        [histCF,~,~] = histcounts(imP(:),(0:(RGBq))/(RGBq));
    else
        imP = imgRGB(rngR,rngC,:);
        imHSV = rgb2hsv(imP);
        imH = imHSV(:,:,1);
        imS = imHSV(:,:,2);
        imV = imHSV(:,:,3);
        imH(imH==0)=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Hq = 8; %8
        Sq = 2; %2
        Vq = 2; %2
        [~,~,binH] = histcounts(imH(:),(0:(Hq))/(Hq));
        [~,~,binS] = histcounts(imS(:),(0:(Sq))/(Sq));
        [~,~,binV] = histcounts(imV(:),(0:(Vq))/(Vq));
        histC = zeros(Hq,Sq,Vq);
        for b=1:numel(binH)
            histC(binH(b),binS(b),binV(b)) = ...
                histC(binH(b),binS(b),binV(b))+1;
        end
        histCF = histC(:);
    end
    histCF = histCF ./ sum(histCF);
    hC{i} = histCF;
end
%%
wavData.c = sc;
wavData.m = m;
wavData.a = a;
wavData.x = (x-(size(img,1)/2))/max(size(img));
wavData.y = (y-(size(img,2)/2))/max(size(img));
wavData.xS = x;
wavData.yS = y;
wavData.v = v;
wavData.vi = vi;
wavData.hC = hC;
end

