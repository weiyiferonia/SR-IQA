function [ref,sr,scale,MOSc,DMOSc,inor,inos] = read_SRD(index, database)
%ref means the reference image
%sr means the super-resolution image
%scale means the magnification of super-resolution image
%inor means the name of reference image
%inos means the name of super-resolution image
dbn = database.name; % database name
dbpath = ['G:/Database_IQA/database/' dbn '/']; % database path 
%% 
%for super-resolution images database
%for super-resolution images database
if strcmpi(dbn,'SRD') 
    DMOSc = []; 
    pathsr = [dbpath 'super-resolution_images/']; % path for super-resolution images
    pathref = [dbpath 'reference_images/']; % path for reference images
    MOSwithnames = importdata([dbpath 'mos_with_names.txt'],'r'); % load mos with names
    ncti = MOSwithnames{index}(end-13:end); % name of the current tested image
    sr = imread([pathsr ncti]); % read the super-resolution image for testing
    ncri = [ncti(1:5) '.bmp']; % name of the current reference image
    ref = imread([pathref ncri]); % read the reference image
    scale = str2double(ncti(7)); % return the scale
    MOSc = str2double(MOSwithnames{index}(1:6));
    inor = ncri(1:end-4);%return the reference image name 
    inos = ncti(1:end-4);%return the super-resolution image name
    return;
end