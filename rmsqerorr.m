function RMSE = rmsqerorr(x,y)
% This function is used to calucate RMSE

RMSE = sqrt(mean((x-y).^2));