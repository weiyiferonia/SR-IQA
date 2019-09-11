function [CC,RMSE,SROCC,KROCC,beta] = resultevaluation(OS,MOS)
% This function is used to evaluate the results obtained from the algorithm
% Input: OS！！objective sorce; MOS！！can be MOS or DMOS
% outout: CC！！Pearson linear correlation coefficient
%         RMSE！！root mean squared error
%         SROCC！！Spearman rank-order correlation coefficient
%         KROCC！！Kendall rank-order correlation coefficient
%% nonlinear fitting
beta0 = initia(OS,MOS);%[1 100 1 0 -1]; %0.5*rand(1,5); % initial value for logistic function 
% beta0 = [max(MOS) 20 mean(OS) 0.1 0.1];
logisf = inline('beta0(1).*(0.5-1./(1+exp(beta0(2).*(x-beta0(3)))))+beta0(4).*x+beta0(5)','beta0','x'); % expression of logistic function
% finalbeta = nlinfit(SS,MOS,logisf,beta); % calculate fitted coefficient
options.MaxIter = 500000;
options.MaxFunEvals = 500000;
beta = lsqcurvefit(logisf,beta0,OS,MOS,[],[],options);
wplot = 0;%true;false; %whether to plot scatter  diagrams or not
if wplot==1
    % to plot 
    figure;plot(OS,MOS,'g+');hold on;
    xx = min(OS):0.001:max(OS); 
    plot(xx,logis(xx,beta),'r');hold off;
end
%% calculate various criteria
OSnr = logis(OS,beta); % objective sorce after nonlinear regression
CC = PearsonCC(OSnr, MOS);
RMSE = rmsqerorr(OSnr, MOS);
SROCC = SpearmanCC(OSnr, MOS);
KROCC = KendallCC(OSnr, MOS);

