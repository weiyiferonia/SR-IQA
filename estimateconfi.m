function Confi = estimateconfi(I,var_patch_size)
% 估计S-T分解的自适应权重: 结构边缘处值大，其余的值小

if ndims(I) ==3
    YCbCr = rgb2ycbcr(I);
    I = YCbCr(:, :, 1);
end
I = double(I);
[h, w] = size(I);
SE = strel('disk',3); %square
J = imerode(imdilate(I,SE),SE);
P = 255 - imerode(imdilate(255-I,SE),SE);
Confi1 =  ones(h,w);
Confi2 =  Confi1;
padded_Gxn = zeros(h+var_patch_size-1,w+var_patch_size-1);
padded_Gyn = padded_Gxn;
padded_GxJ = padded_Gxn;
padded_GyJ = padded_Gxn;
padded_GxP = padded_Gxn;
padded_GyP = padded_Gxn;
% padded_Ga  =  padded_Gm;
r =  (var_patch_size-1)/2;
state_point = r + 1;
% padded_I(state_point : h+state_point-1,state_point : w+state_point-1) = I;
%% 计算梯度
f_x = [0 0 0 0 0; 1 2 3 2 1; 0 0 0 0 0; -1 -2 -3 -2 -1; 0 0 0 0 0];
Gxn = imfilter(I, f_x, 'replicate');
Gyn = imfilter(I, f_x', 'replicate');
GxJ = imfilter(J, f_x, 'replicate');
GyJ = imfilter(J, f_x', 'replicate');
GxP = imfilter(P, f_x, 'replicate');
GyP = imfilter(P, f_x', 'replicate');
% [Gx, Gy] = gradient(I);
Gm = (Gxn.^2 + Gyn.^2).^0.5; % 梯度幅度
max_Gm = max(Gm(:));
% Ga = atan2(Gy, Gx); % 梯度角度
% Gxn = Gx;%./(Gm+1e-7);
% Gyn = Gy;%./(Gm+1e-7);
padded_Gxn(state_point : h+state_point-1,state_point : w+state_point-1) = Gxn;
padded_Gyn(state_point : h+state_point-1,state_point : w+state_point-1) = Gyn;
padded_GyJ(state_point : h+state_point-1,state_point : w+state_point-1) = GyJ;
padded_GxJ(state_point : h+state_point-1,state_point : w+state_point-1) = GxJ;
padded_GyP(state_point : h+state_point-1,state_point : w+state_point-1) = GyP;
padded_GxP(state_point : h+state_point-1,state_point : w+state_point-1) = GxP;

for ii = state_point: h+state_point-1
    for jj = state_point: w+state_point-1

        v_Gxn = reshape(padded_Gxn(ii-r:ii+r,jj-r:jj+r ), [var_patch_size*var_patch_size,1]);
        v_Gyn = reshape(padded_Gyn(ii-r:ii+r,jj-r:jj+r ), [var_patch_size*var_patch_size,1]);
        v_GxJ = reshape(padded_GxJ(ii-r:ii+r,jj-r:jj+r ), [var_patch_size*var_patch_size,1]);
        v_GyJ = reshape(padded_GyJ(ii-r:ii+r,jj-r:jj+r ), [var_patch_size*var_patch_size,1]);
        v_GxP = reshape(padded_GxP(ii-r:ii+r,jj-r:jj+r ), [var_patch_size*var_patch_size,1]);
        v_GyP = reshape(padded_GyP(ii-r:ii+r,jj-r:jj+r ), [var_patch_size*var_patch_size,1]);

        tensor = [v_Gxn'*v_Gxn, v_Gxn'*v_Gyn;v_Gxn'*v_Gyn,v_Gyn'*v_Gyn]; 
        [~,eigenvalues] = eig(tensor);
        latent = sort(diag(eigenvalues),'descend');
        aniso =  (latent(1)-latent(2))/(latent(1)+latent(2)+1e-7);
       
        tensorJ = [ v_GxJ'*v_GxJ, v_GxJ'*v_GyJ;v_GxJ'*v_GyJ,v_GyJ'*v_GyJ]; 
         [~,eigenvaluesJ] = eig(tensorJ);
        latentJ = sort(diag(eigenvaluesJ),'descend');
        anisoJ =  (latentJ (1)-latentJ (2)+1e-4)/(latentJ (1)+latentJ (2)+1e-4);
        
        tensorP = [ v_GxP'*v_GxP, v_GxP'*v_GyP;v_GxP'*v_GyP,v_GyP'*v_GyP]; 
        [~,eigenvaluesP] = eig(tensorP);
        latentP = sort(diag(eigenvaluesP),'descend');
        anisoP =  (latentP (1)-latentP (2)+1e-4)/(latentP (1)+latentP (2)+1e-4);
        
        aniso_changeJ = (min(anisoJ,aniso)+1e-4)/(max(anisoJ,aniso)+1e-4);
        aniso_changeP = (min(anisoP,aniso)+1e-4)/(max(anisoP,aniso)+1e-4);
        aniso_change = max(aniso_changeJ,aniso_changeP);
        
        Confi1(ii-r,jj-r) = aniso;
        Confi2(ii-r,jj-r) = aniso_change;

    end
end

Gm_n = (Gm+1e-3)/(max_Gm+1e-3);
Gm_n =imdilate(Gm_n ,SE);
Gm_n =imerode(Gm_n ,SE);
Confi = Confi1.*Confi2.*Gm_n;
Confi = (Confi+1e-3)/(max(Confi(:))+1e-3);

end