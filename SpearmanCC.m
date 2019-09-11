function SROCC = SpearmanCC(x,y)
% This function is used to calucate Spearman rank-order correlation coefficient
Num = numel(x);
[temp,index1] = sort(x);
[temp,index2] = sort(y);
[temp,index1] = sort(index1);
[temp,index2] = sort(index2);
D = index1-index2;

SROCC = 1-6*(D'*D)/(Num*(Num *Num -1));