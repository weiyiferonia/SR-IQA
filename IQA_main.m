clc;
clear;
run('G:/vlfeat/vlfeat/toolbox/vl_setup')
cdb.name = 'SRD';
cdb.srimg = 980;
cdb.ss = [];
cdb.os = [];
os_H = [];
os_T = [];
os_S = [];
fprintf('Starting....\n');
for ii = 1:cdb.srimg
    [ref,sr,scale,MOSc,DMOSc,RefName,SRName] = read_SRD(ii, cdb);
    if isempty(MOSc)   % get the MOS/DMOS
           cdb.ss = [cdb.ss; DMOSc]; % stack the MOS
    else 
           cdb.ss = [cdb.ss; MOSc];  
    end
    [iqa,iqa_H,iqa_T,iqa_S] = SRsim(RefName,SRName);
    cdb.os = [cdb.os; iqa];  
    os_H = [os_H; iqa_H]; 
    os_T = [os_T; iqa_T];
    os_S = [os_S; iqa_S];
    fprintf('%d-th SR image out of %d ones in the database %s has been assessed... \n',ii, cdb.srimg, cdb.name);
end

[CC, RMSE, SROCC, KROCC] = resultevaluation(cdb.os,cdb.ss);   
[CC_H, RMSE_H, SROCC_H, KROCC_H] = resultevaluation(os_H,cdb.ss); 
[CC_T, RMSE_T, SROCC_T, KROCC_T] = resultevaluation(os_T,cdb.ss);
[CC_S, RMSE_S, SROCC_S, KROCC_S] = resultevaluation(os_S,cdb.ss);