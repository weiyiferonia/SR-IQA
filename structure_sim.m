function [results1,weight] = structure_sim(RefName,SRName)
% RefName: 'img08'
% SRName: 'img08_3_09'

%% Read

refpath = 'G:\Database_IQA\STD\refImg_s\';
srpath = 'G:\Database_IQA\STD\srimg_s\';
load([refpath RefName '_u.mat']);
load([srpath SRName '_u.mat'])
u1 = sr_u;
[M, N] = size(u);

%% Config
% C = 0.1;
patch_size = 5; 
padding_size = (patch_size-1)/2;
f_x = [0 0 0 0 0; 1 2 3 2 1; 0 0 0 0 0; -1 -2 -3 -2 -1; 0 0 0 0 0];
u_pad = padarray(u, [padding_size padding_size], 'replicate');
u1_pad = padarray(u1, [padding_size padding_size], 'replicate');
Gx = imfilter(u_pad, f_x, 'replicate');
Gy = imfilter(u_pad, f_x', 'replicate');
Gx1 = imfilter(u1_pad, f_x, 'replicate');
Gy1 = imfilter(u1_pad, f_x', 'replicate');
Gmax= max(max((Gx.^2+Gy.^2)).^0.5);%max(max(Gx(:)),max(Gy(:)));
results1 = zeros(M,N);
weight = zeros(M,N);
% results2 = zeros(M,N);

%% Calculate local structural similarity
for ii = 1:M
    for jj = 1:N
        u_Gx = reshape(Gx(ii:ii+patch_size-1,jj:jj+patch_size-1 ), [patch_size*patch_size,1]); 
        u1_Gx = reshape(Gx1(ii:ii+patch_size-1,jj:jj+patch_size-1 ), [patch_size*patch_size,1]);

        u_Gy = reshape(Gy(ii:ii+patch_size-1,jj:jj+patch_size-1 ), [patch_size*patch_size,1]); 
        u1_Gy = reshape(Gy1(ii:ii+patch_size-1,jj:jj+patch_size-1 ), [patch_size*patch_size,1]);
        
        tensor = [u_Gx'*u_Gx, u_Gx'*u_Gy; u_Gx'*u_Gy, u_Gy'*u_Gy]; 
        tensor1 = [u1_Gx'*u1_Gx, u1_Gx'*u1_Gy; u1_Gx'*u1_Gy, u1_Gy'*u1_Gy]; 
        
        [eigenvectors,eigenvalues] = eig(tensor);
        [eigenvectors1,eigenvalues1] = eig(tensor1);
        
        [~,index]= sort(diag(eigenvalues),'descend');
        [~,index1]= sort(diag(eigenvalues1),'descend');

        
        d = eigenvectors(:,index(2));% Maximum direction
        d1 = eigenvectors1(:,index1(2));
          
        Gm = (Gx(ii+padding_size,jj+padding_size).^2+ Gy(ii+padding_size,jj+padding_size).^2).^0.5;
        Gm1 = (Gx1(ii+padding_size,jj+padding_size).^2+ Gy1(ii+padding_size,jj+padding_size).^2).^0.5;
        cof1 = 1/(min(Gm,Gm1)/2295+1e-3);
        results1(ii,jj) = ((abs(d'*d1))+cof1)/(1+cof1);
        weight(ii,jj) = max(Gm,Gm1);
    end
end
weight = weight./(sum(weight(:)));
if numel(find(isnan(results1)))
    error('ZF2:Something wrong in structure part!!!!');
end
