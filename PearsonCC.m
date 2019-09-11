function CC = PearsonCC(x,y)
% This function is used to calucate Pearson linear correlation coefficient
% R = corrcoef(x,y);
% CC = R(1,2);
mx = mean(x);
my = mean(y);
CC =((x'-mx)*(y-my)+1e-6)/(sqrt((x'-mx)*(x-mx)*(y'-my)*(y-my))+1e-6);