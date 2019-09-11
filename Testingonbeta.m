clc;
clear;
load('os.mat');
load('cdb.mat');
SROCC_list = [];
for ii = 0: 0.1:10
    beta = ii;
    os_B = (os_H.*os_S).^beta.*os_T;
    [CC, RMSE, SROCC, KROCC] = resultevaluation(os_B,cdb.ss); 
    SROCC_list = [SROCC_list;SROCC];
end
x = 0: 0.1:10;
plot(x,SROCC_list);