clc
clear all
% 
theta = [43.83];
lambda = linspace(0.73,0.8,50)*1e-6;
[p_re, p_Im,REF]= myTMM_Partial_R_eps (lambda,theta);

%%
%E = linspace(1.5,1.7,200);
count  = 1; 
Te = [700,800,1000];
figure
for ii = 1:size(Te,2)
[dR_R(:,count)] = mydR_R_WithAngle (Te(ii),p_re, p_Im, theta,REF,lambda);
plot(lambda, dR_R(:,count)); hold on
count = count+1;
end
%%
%[~, dR_R1] = mydR_R_0angle (E, Te);
Te = [800,900,1000];
figure; 
for i = 1: size(Te,2)
[p_re(i,:), p_Im(i,:),REF(i,:)]= myTMM_Partial_R_eps (lambda,theta);
%[dR_R] = mydR_R_WithAngle (Te(i),p_re, p_Im, theta,REF,lambda);
plot(p_Im(i,:)); hold on
end