function [results1, weight] = highf_sim(RefName,SRName)
% RefName: 比如'img19'
% SRName: 比如'img19_3_12' or 'img19_2_01'
%% Read
% RefName = 'img19';%'img08'; ppt & ppt 02_01
% SRName1 = 'img19_3_12';% 'img08_3_09';% 0.61100
% SRName2 = 'img19_2_01';%'img08_3_10'; % 0.38300

% refpath = 'G:\Database_IQA\database\SRD\reference_images\';
% srpath = 'G:\Database_IQA\database\SRD\super-resolution_images\';
% I = imread([refpath RefName '.bmp']);
% I1 = imread([srpath SRName '.bmp']);
% I2 = imread([srpath SRName2 '.bmp']);
% Y = rgb2ycbcr(I);
% Y = double(Y(:,:,1))/255;
% Y1 = rgb2ycbcr(I1);
% Y1 = double(Y1(:,:,1))/255;
% Y2 = rgb2ycbcr(I2);
% Y2 = double(Y2(:,:,1))/255;
% Y1 = Y1 - mean(Y1(:)) + mean(Y(:));
% Y2 = Y2 - mean(Y2(:)) + mean(Y(:));
refpath = 'G:\Database_IQA\STD\refImg_s\';
refpath_v = 'G:\Database_IQA\STD\refImg_t\';
% srpath = 'G:\Database_IQA\STD\srimg_s\';
srpath_I = 'G:\Database_IQA\database\SRD\super-resolution_images\';
load([refpath RefName '_u.mat']);
load([refpath_v RefName '_v.mat']);
I = imread([srpath_I  SRName '.bmp']);
% load([srpath SRName '_u.mat'])
% u1 = sr_u;
if ndims(I) ==3
    YCbCr = rgb2ycbcr(I);
    u1 = double(YCbCr(:, :, 1))-v;
end  
Y = u;
Y1 = u1 - mean(u1(:)) + mean(u(:));

%% Config
C = 2; %1e-5;
sigma = 5;
window = sigma*2+1;
% patch_size = window;
% length_sigma = length(sigma);

%% 准备滤波器和滤波器响应
% mm = 1;
LF  =fspecial('gaussian', window, sigma);

L =imfilter(Y,LF,'replicate'); % 
L1=imfilter(Y1,LF,'replicate'); % psnr(out,out1) = 56.8974;
% L2=imfilter(Y2,LF,'replicate'); % psnr(out,out2) = 54.5598;
H = Y - L;
H1 = Y1 - L1;
% H2 = Y2 - L2;

EF = fspecial('disk', 3);
% E = imfilter(H.^2,EF,'replicate');
% E1 = imfilter(H1.^2,EF,'replicate');
E = imfilter(abs(H),EF,'replicate');
E1 = imfilter(abs(H1),EF,'replicate');
% E2 = imfilter(H2.^2,EF,'replicate');
% 

% E = H.^2;
% E1 = H1.^2;
% E2 = H2.^2;
E1 = min(E,E1);
results1 = (2.*E.*E1+C)./(E.^2+E1.^2+C); 
weight = E/sum(E(:));

% results1 = exp((E1-E)./(0.5*E1+0.5+E));%'img04_2_06'
% results2 = (2.*E.*E2+C)./(E.^2+E2.^2+C);


% for kk = length_sigma:-1:3
%     mm = mm+1;
%     F{mm} = fspecial('gaussian', window, sigma(kk-1)) - fspecial('gaussian', window, sigma(kk));
%     L{mm}=imfilter(Y,F{mm},'replicate');
%     L1{mm}=imfilter(Y1,F{mm},'replicate');
%     L2{mm}=imfilter(Y2,F{mm},'replicate');   
% end
% if mm~= length_sigma-1
%     print('ZF: Something wrong in the number of filters!');
% end
% F{length_sigma} = fspecial('gaussian', window, sigma(2));
% L{length_sigma}= Y - imfilter(Y,F{length_sigma},'replicate');
% L1{length_sigma}= Y1 - imfilter(Y1,F{length_sigma},'replicate');
% L2{length_sigma}= Y2 - imfilter(Y2,F{length_sigma},'replicate');
%% 
if numel(find(isnan(results1)))
    error('ZF2:Something wrong in high-freq part!!!!');
end
% mean(results1(:))
% mean(results2(:))
% imtool(results1,[])
% imtool(results2,[])

%当前单独结果SROCC，mean pooling是0.8240；std是0.8590;
% 过好的估计:效果比较差的是锯齿效应比较明显的（因为锯齿效应也有很丰富的高频成分）。
% 对纯模糊的，较低的估计了
