clc;
clear;
u_dir =  'G:\Database_IQA\database\LIVE\refimgs\u\';
v_dir = 'G:\Database_IQA\database\LIVE\refimgs\v\';
inputs = dir(fullfile(u_dir, '*.mat'));
alpha_H = zeros(length(inputs),1);
for k = 1:length(inputs)
    fprintf('Processing %d\n',k);
    u_name = inputs(k).name;
    v_name = ['v' u_name(2:length(u_name)-5) 'v.mat'];

    load([u_dir u_name]);
    load([v_dir v_name]);
    u_e = abs(u); 
    v_e = abs(v);
    alpha_H(k) = log(mean(u_e(:)))/log(mean(v_e(:))); %web
end
mean(alpha_H(:))
