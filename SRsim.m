function [result,result_H, result_T, result_S] = SRsim(RefName,SRName)

alpha_H =  3.9709;
alpha_T = 1;
alpha_S =  3.9709;% Web-Fechner Law: P = k*lnS + C

[H_sim, w_h] = highf_sim(RefName,SRName);
[T_sim,w_T] =  texture_sim(RefName,SRName);
[S_sim, w_s] =  structure_sim(RefName,SRName);

%% Pooling
% result = mean(map(:));
result_T = sum(sum(w_T.*T_sim));
result_H = sum(sum(w_h.*H_sim));%1 - std(H_sim(:));
result_S = sum(sum(w_s.*S_sim));%1 - std(S_sim(:));
result = (result_H.^alpha_H) * (result_T.^alpha_T)*(result_S.^alpha_S);