% This is the C implementation of the paper "Color-guided Depth Recovery from RGB-D Data 
% Using an Adaptive Auto-regressive Model" in TIP 2014. This implementation
% is the same as the original MATLAB source code which can be download
% here:http://cs.tju.edu.cn/faculty/likun/projects/depth_recovery/index.htm.
% The difference is that our implementation uses C and sparse matrix
% which can be much faster and can handle much larger depth maps..


clear
clc
close all


rWin = 4; % the radius of the AR neighborhood system N(x) in Eq.(5)
rPatch = 3; %the radius of the patch to compute the shape-based color term of the AR coefficient.

%%%%%% the original parameter setting in the source code provided by the author
sigma1=4; %depth weighted
sigma_c=0.35;
sigma2=sigma_c*1.732*11; %color weighted
% sigma3=3; %gausian weighted, however, this parameter is not used in their MATLAB
% code
sigma4=0.1;%bi weighted

mu = 200; % default value, try 0.1 for noisy data to see better results. we use 0.1 for noisy data in our experiments.

Data = load('ARProvideSample.mat');
% Data = load('ToFSimulated_2.mat');
% Data = load('ToFReal.mat');

DepthGT = Data.DepthGT;
Color = Data.Color;
DepthGuide = Data.DepthGuide;
DepthSample = Data.DepthSample;

%%%%%%%%%%
% dMax = max(DepthGuide(:));
% DepthGuide = DepthGuide/dMax*255; 
% %%%%% uncomment these two columns if Data = load('ToFReal.mat');
%%%%%%%%%%%%

Color = double(image_rgb2yuv(Color));
Color = padarray(Color, [rPatch, rPatch], 'symmetric');


[m, n] = size(DepthGT);

Index = zeros(m, n);
Index(DepthSample>0) = 1;
PtP = spdiags(Index(:), 0, m*n, m*n);

B = DepthSample(:);

t11 = tic;
Q = mexGetARWeight(Color(:,:,1), Color(:,:,2), Color(:,:,3), DepthGuide, rWin, rPatch, sigma1, sigma2, sigma4);
QtQ = Q*Q';
t12 = toc(t11);
fprintf('Computing the AR weights costs %f seconds\n', t12);
clear Q;


t11 = tic;
X=(PtP+(1/mu)*QtQ)\B; 
t12 = toc(t11);
fprintf('Solving the system costs %f seconds\n', t12);
Res = reshape(X, m, n);

Mask = zeros(m, n);
Mask(DepthGT>0) = 1;
DepthGT(Mask<1) = 0;
Diff = abs(DepthGT - Res).*Mask;
MAE = sum(Diff(:))/sum(Mask(:));

fprintf('MAE is %f\n', MAE);

dMax = max(Res(:));
figure
imagesc([DepthGT, Res])
colormap gray
% imshow([DepthGT, Res]./dMax)

