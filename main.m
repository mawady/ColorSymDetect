clc; clear all; close all; warning off;
%%
addpath(genpath(fullfile('.','libs')));
%%
% srcDir = '../data/input/dataset_AVA';
% SymGT = load('../data/GT/GT_AVA/AVA_GT.mat');
% FileNames = SymGT.new_cell_GT(:,1);
% SymEll = SymGT.new_cell_GT(:,2);
% i = 221;
% file = strtrim(FileNames{i});
%%
srcDir = './input';
file = '16057.jpg';
%%
[~,name,~] = fileparts(file);
img = imread([srcDir '/' file]);
disp(['Processing : ' name]);
tic;
SymOcLgHSV = symBilOurCentLogGaborHSV(img);
toc;
%%
Num = 10;

MarkerEdgeColors=hsv(Num);
str = {};
mkSize = 10;
lnWidth = 2;

figure;
imshow(img); hold on;
for j=1:min(Num,size(SymOcLgHSV,1))
    X = [SymOcLgHSV(j,1) SymOcLgHSV(j,3)];
    Y = [SymOcLgHSV(j,2) SymOcLgHSV(j,4)];
    plot(X,Y, '-s','Color',MarkerEdgeColors(j,:),...
              'MarkerSize',mkSize,'LineWidth',lnWidth);
    str{j} = [num2str(j) ' - ' num2str(SymOcLgHSV(j,5))];
end
legend(str,'Location','BestOutside');
hold off;