%% Initialization
clc;
clear;
close all;
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();

% Paths and software
csvFile = 'AgesVisits.csv';
nccam2Root = '/study/rinpoche/MRAnalysis/NCCAM2Data';
maskFile = sprintf('%s/MaskForBrainAGE.nii', nccam2Root);
myrRoot = '/study/rinpoche/Processed';
myr2005Root = '/study/rinpoche/2005Data/2005';

addpath(genpath('drtoolbox'));
addpath(genpath('spider'));

%% Control data
csv = table2cell(readtable(csvFile, 'Delimiter', ','));
mask = load_nifti(maskFile);
maskIdx = find(mask.vol(:) > 0);
numVoxels = length(maskIdx);

imgFiles = cellfun(@(id, visit) sprintf('%s/TimePoint%s/swmwc1o%s_%s_T1High+orig.nii', nccam2Root, visit, id, visit), csv(1:239, 1), csv(1:239, 4), 'UniformOutput', false);
gmVols = cellfun(@(img) load_nifti(img).vol, imgFiles, 'UniformOutput', false);
CtrlX = cell2mat(cellfun(@(gm) gm(maskIdx)', gmVols, 'UniformOutput', false));
CtrlY = cell2mat(csv(1:239, 5));

%% MYR data
MYRX = zeros(4, numVoxels);

% TP1 -- 2002
img2002 = zeros(3, numVoxels);
% Day 1
imgFile = sprintf('%s/swmwc1o_SP3_day1_orig_2002.nii', myrRoot);
img = load_nifti(imgFile);
img2002(1, :) = img.vol(maskIdx);
% Day 2
imgFile = sprintf('%s/swmwc1o_SP3_day2_orig_2002.nii', myrRoot);
img = load_nifti(imgFile);
img2002(2, :) = img.vol(maskIdx);
% Day 3
imgFile = sprintf('%s/swmwc1o_SP3_day3_orig_2002.nii', myrRoot);
img = load_nifti(imgFile);
img2002(3, :) = img.vol(maskIdx);
% Mean
MYRX(1, :) = mean(img2002, 1);

% TP2 -- 2005
img2005 = zeros(2, numVoxels);
% Day 1
imgFile = sprintf('%s/swmwc1oSP17_day1_orig.nii', myr2005Root);
img = load_nifti(imgFile);
img2005(1, :) = img.vol(maskIdx);
% Day 2
imgFile = sprintf('%s/swmwc1oSP17_day1_orig.nii', myr2005Root);
img = load_nifti(imgFile);
img2005(2, :) = img.vol(maskIdx);
% Mean
MYRX(2, :) = mean(img2005, 1);

% TP3 -- 2007
imgFile = sprintf('%s/2007/swmwc1o_EFGRE3D.nii', myrRoot);
img = load_nifti(imgFile);
MYRX(3, :) = img.vol(maskIdx);

%TP4 -- 2016
imgFile = sprintf('%s/2016/swmwc1o_MR_PU_2016.nii', myrRoot);
img = load_nifti(imgFile);
MYRX(4, :) = img.vol(maskIdx);
MYRY = [27; 30; 32; 41];

%% Saving the data
CtrlData = data(CtrlX, CtrlY);
MYRData = data(MYRX, MYRY);
save('ControlData.mat', 'CtrlData')
save('MYRData.mat', 'MYRData')

%% Loading the XY
load('ControlData.mat');
load('MYRData.mat');
global allX
global allY
allX = [CtrlData.X; MYRData.X];
allY = [CtrlData.Y; MYRData.Y];
clear ControData MYRData;
numVols = size(allX, 1);

%% RVM
BrainAGEExp = cell2mat(arrayfun(@(i) TrainTest(i, numVols), 1:numVols, 'UniformOutput', false)');

% Sanity check
iif(sum(allY-(BrainAGEExp(:, 2)-BrainAGEExp(:, 1))) == 0, 'OK',...
    true, 'Check again')

% Save results
BrainAGEExp = [BrainAGEExp, allY];
save('BrainAGE_Results.csv', 'BrainAGEExp', '-ascii');

%% function TrainTest
function [BrainAGEExp] = TrainTest(i, numVols)
tic;
global allX
global allY
% Training
trainX = allX;
trainY = allY;
trainX(i, :) = [];
trainY(i) = [];
trainD = data(trainX, trainY);
[~, trainedModel] = train(relvm_r(kernel('poly', 1)), trainD);
% Testing
testD = data(allX(i, :), allY(i));
predTest = test(trainedModel, testD);
% BrainAGE (\Delta)
BrainAGE = predTest.X - predTest.Y;
% Estimated brain age
estAge = predTest.X;
BrainAGEExp = [BrainAGE estAge];
fprintf('BrainAGE[%d/%d]=%f, EstAge[%d/%d]=%f, CA[%d/%d]=%f in %f sec.\n', i, numVols, BrainAGE,...
    i, numVols, estAge, i, numVols, allY(i), toc);
end