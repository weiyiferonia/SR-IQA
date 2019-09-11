function KROCC = KendallCC(x,y)
% This function is used to calucate Kendall rank-order correlation coefficient
Num = numel(x);
tau = 0;
for i = 1:Num
    tau = tau + sum(sign(x(i:Num)-x(i)).*sign(y(i:Num)-y(i)));
end
KROCC = 2*tau/(Num*(Num-1));