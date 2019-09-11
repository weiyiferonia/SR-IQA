function [results1,array_weight] = texture_sim(RefName,SRName)
% RefName: 比如'img16'
% SRName: 比如'img16_2_20'
%% Config

binsize = 4; % SFIT特征的bin size
G_N = 5; % 高斯窗尺寸
sigma = 1; % 高斯窗方差
refpath = 'G:\Database_IQA\STD\refImg_t\';
srpath = 'G:\Database_IQA\STD\srimg_t\';

firstcenter = 1+1.5*binsize;
padding_size = firstcenter-1;
%pathsize = padding_size*2+1;
%% data pre-processing
% srpath = 'G:\Database_IQA\database\SRD\super-resolution_images\';
% refpath = 'G:\Database_IQA\database\SRD\reference_images\';
% I = rgb2gray(imread([refpath 'img16.bmp']));
load([refpath RefName '_v.mat']);
% I1 = rgb2gray(imread([srpath 'img16_2_20.bmp'])); %MOS 0.90996
load([srpath SRName '_v.mat'])%;'G:\Database_IQA\STD\srimg_t\img16_2_20_v.mat')
v1 = sr_v;%_16_2_20
% I2 = rgb2gray(imread([srpath 'img16_4_21.bmp']));%MOS 0.70047
% load('G:\Database_IQA\STD\srimg_t\img16_4_21_v.mat')
% v2 = sr_v;%_16_4_21
[M, N] = size(v);
% F = extractLBPFeatures(v(230:236,40:46),'Upright',false);
% Padding the input
v_pad = padarray(v, [padding_size padding_size], 'replicate');
v1_pad = padarray(v1, [padding_size padding_size], 'replicate');
% v2_pad = padarray(v2, [padding_size padding_size], 'replicate');

% Adding Gaussian before SIFT
gausFilter = fspecial('gaussian',[G_N G_N],sigma); 
v_pad = imfilter(v_pad,gausFilter,'same');  
v1_pad = imfilter(v1_pad,gausFilter,'same'); 
% v2_pad = imfilter(v2_pad,gausFilter,'same'); 

if ~isequal(size(v_pad),size(v1_pad))%||~isequal(size(v_pad),size(v2_pad))
    error('ZF0:Size does not match!!!!');
end
% v1_results = zeros([M N]);
% v2_results = zeros([M N]);
%% 计算相似性
% LBP test：
% for ii = 1:M
%     for jj = 1:N      
%         rang_x = ii:ii+pathsize-1;
%         rang_y = jj:jj+pathsize-1;
% %         F = extractLBPFeatures(v_pad(rang_x,rang_y)); %LBP直方图 ,'Upright',false
% %         F1 = extractLBPFeatures(v1_pad(rang_x,rang_y)); 
% %         F2 = extractLBPFeatures(v2_pad(rang_x,rang_y)); 
% 
%         v1_results(ii,jj)=F*F1';%sum(abs(F-F1));
%         v2_results(ii,jj)=F*F2';%sum(abs(F-F2));
%     end
% end
% SIFT test：
[F,D] = vl_dsift(single(v_pad),'size',binsize,'Fast','FloatDescriptors');
[~,D1] = vl_dsift(single(v1_pad),'size',binsize,'Fast','FloatDescriptors');
% [F2,D2] = vl_dsift(single(v2_pad),'size',binsize,'Fast','FloatDescriptors');
K = max(F(1,:))-min(F(1,:))+1;% 最大的维数
L = max(F(2,:))-min(F(2,:))+1;
if K~=N||L~=M
    error('ZF1:Size does not match!!!!');
end
results1 = zeros(L,K);
array_weight = zeros(L,K);
% results2 = zeros(L,K);
% shift = padding_size;%min(F(1,:))-1;
% Cmap1 = zeros(L,K);
% Cmap2 = zeros(L,K);
for ii = 1:L % L=380
    for jj = 1:K %K=500
        index = (jj-1)*L +ii;        
        %
        local_patch_ii = ii:ii+padding_size-1;
        local_patch_jj = jj:jj+padding_size-1;
        local_patch = v_pad(local_patch_ii,local_patch_jj);
        local_patch1 = v1_pad(local_patch_ii,local_patch_jj);
%         local_patch2 = v2_pad(local_patch_ii,local_patch_jj);
        var0 = var(local_patch(:))+1e-3;
        var1 = var(local_patch1(:))+1e-3;
%         var2 = var(local_patch2(:))+1e-3; % var 范围从10-3到10+3
        
        C1 = 1/max(var0,var1);
%         C2 = 1/max(var0,var2);
%         Cmap1(ii,jj) = C1;
%         Cmap2(ii,jj) = C2;
        fea = D(:,index);% norm0 = norm(D(:,index));
        fea_n = fea/(norm(fea)+1e-3);
        fea1 = D1(:,index); 
        fea_n1 = fea1/(norm(fea1)+1e-3);
%         fea2 = D2(:,index); 
%         fea_n2 = fea2/norm(fea2);
        results1(ii,jj) = (fea_n'*fea_n1+C1)/(1+C1);
        array_weight(ii,jj) = max(var0,var1);
%         results2(ii,jj) = (fea_n'*fea_n2+C2)/(1+C2);
    end
end

array_weight = array_weight./(sum(array_weight(:)));
% results1 = array_weight.*results1;
if numel(find(isnan(results1)))&&numel(find(isnan(array_weight)))
    error('ZF2:Something wrong in texture part!!!!');
end
% imtool(results1,[])
% imtool(results2,[])
% mean(results1(:))
% mean(results2(:))
% imtool(Cmap1,[])
% imtool(Cmap2,[])

%当前单独结果SROCC，weighted-mean pooling是0.8608；
%  SIFT for keypoint: In order to achieve orientation invariance, the coordinates of the descriptor and the gradient orientations are rotated relative
%  to the keypoint orientation. In this work, we do not rotate the gradient
%  orientations according to the center pixel.
