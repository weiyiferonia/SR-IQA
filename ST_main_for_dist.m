clear;
clc;
%% Config
lambda = 1;
patch_size = 50;
overlap = 10;

%% Path
input_dir = 'G:/Database_IQA/database/SRD/super-resolution_images/';
save_u = 'srimg_s/';
save_v = 'srimg_t/';
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
    [sr_u,sr_v] = decompose(Y, patch_size, overlap, lambda, 1e-8);%, 0, 1-Confi);

    %% Store the u and v components
    save([save_u Imgname '_u.mat'],'sr_u');
    save([save_v Imgname '_v.mat'],'sr_v');


end




