function maxData = computeVotingMax(voteData,maxParam)

d = maxParam.halfwindow;
voteMapEx = [voteData.voteMap(end-d+1:end, :); ...
       voteData.voteMap;...
       voteData.voteMap(1:d, :)];
   
H = fspecial('gaussian',maxParam.hsize,maxParam.hsize/4);
I = conv2(double(voteMapEx),H,'same');
voteMapBlur = I(d+1:end-d,:);
I = I/max(max(I));

[nr,nc] = size(I);
J = zeros(nr,nc);
n = 8;
for i = d+1:nr-d
    for j = d+1:nc-d
        values1 = zeros(1,n);
        values2 = zeros(1,n);
        values3 = zeros(1,n);
        for k = 1:n
            ag = (k-1)/n*2*pi;
            vector = round(d*[cos(ag) sin(ag)]);
            values3(k) = I(i+vector(1),j+vector(2));
            vector = round(d/4*[cos(ag) sin(ag)]);
            values2(k) = I(i+vector(1),j+vector(2));
            vector = round(d/8*[cos(ag) sin(ag)]);
            values1(k) = I(i+vector(1),j+vector(2));
        end
        center = I(i,j);
        if (center > max([values1 values2 values3]) && ...
                center > maxParam.lowerbound)
            J(i,j) = 1;
        end
    end
end

cc = bwconncomp(J,8);
stats = regionprops(cc,'Area','Centroid');
l = length(stats);
areas = zeros(1,l);
centers = zeros(2,l);
for i = 1:l
    areas(i) = stats(i).Area;
    centers(:,i) = round([stats(i).Centroid(2) stats(i).Centroid(1)]);
end

% get rid of blobs that are too small
mx = max(areas);
[mn,imn] = min(areas);
while mn < maxParam.minarea*mx % default: 0.25
    areas(imn) = [];
    centers(:,imn) = [];
    [mn,imn] = min(areas);
end

% get rid of blobs for which there's a strong blob nearby
D = distbetcent(centers);
[mn,imn] = min(D);
[mn1,imn1] = min(mn);
while mn1 < maxParam.mindistbetcent
    c = imn1;
    r = imn(c);
    if areas(c) < areas(r)
        idx = c;
    else
        idx = r;
    end
    areas(idx) = [];
    centers(:,idx) = [];
    
    D = distbetcent(centers);
    [mn,imn] = min(D);
    [mn1,imn1] = min(mn);
end

locs = centers;

if(numel(locs)~=0)
    locs(1,:) = locs(1,:) - d;
end
% voteMapBlur = I(d+1:end-d,:);

maxData.locs = locs;
maxData.voteMapBlur = voteMapBlur;

end

function D = distbetcent(centers)
    d = size(centers,2);
    D = zeros(d,d);
    for i = 1:d
        for j = 1:d
            if i == j
                D(i,j) = Inf;
            else
                D(i,j) = norm(centers(:,i)-centers(:,j));
            end
        end
    end
end

