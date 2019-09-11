clear;
clc;
%% Config
lambda = 1;
patch_size = 50;
overlap = 10;
%weight_patch_size = 9;
input_dir = 'G:\Database_IQA\database\LIVE\refimgs\';
save_u = 'G:\Database_IQA\database\LIVE\refimgs\u\';
save_v = 'G:\Database_IQA\database\LIVE\refimgs\v\';
inputs = dir(fullfile(input_dir, '*.bmp'));

 for k= 1:length(inputs)
     fprintf('Processing %d\n',k);
    Imgname = inputs(k).name;
    I = imread([input_dir Imgname]);
    name_length = length(inputs(k).name);
    Imgname = Imgname(1:name_length-4);
    if ndims(I) ==3
       YCbCr = rgb2ycbcr(I);
       Y = YCbCr(:, :, 1);
    end  
    % ST
    [u,v] = decompose(Y, patch_size, overlap, lambda, 1e-8);%, 0, 1-Confi);
    save([save_u Imgname '_u.mat'],'u');
    save([save_v Imgname '_v.mat'],'v');
    
end




