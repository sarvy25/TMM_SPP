function [dR_R] = mydR_R_WithAngle ( Te,p_re, p_Im, theta,R,lambda)  %
% 
% input

% incident angle 
% R: Reflection based on the input incident angle
% lambda : R at given lambda
% p_Im: partial derivitive of reflectivity (calculated from "function
% [p_re, p_Im, REF]= myTMM_Partial_R_eps (lambdaq,theta)" wrt imaginary of
% permitivity 
% p_re: partial derivitive of reflectivity (calculated from "function
% [p_re, p_Im, REF]= myTMM_Partial_R_eps (lambdaq,theta)" wrt real of
% permitivity 

% outputs :
% dR_R : change in the reflectivity
% lambdaq :  interpolated wavelength
% ellipsometry nand k for gold
load n_Au
load k_Au
load lambda_Au
h = 4.13*10^-15;
%h = 6.58 *10^-16; %eV

c = 3*1e8;
E = (h*c)./lambda;
%lambdaq = (h*c)./Ej;
delE = 2.37; % eV
eps0 = 8.85*10^-12;
kb = 8.6*10^-5; %eVK-1

% fermi distribution at 300 kelvin
T0= 300;
a0 = exp((E - delE)./(kb.*T0));
rho0 = 1./(1 + a0);


% calculatin change in the electronics occupancy
for ii = 1: size(Te,2)

    a(ii,:) = exp((E - delE)./(kb.*Te(1,ii)));
    rho(ii,:) = 1./(1 + a(ii,:));
    delrho(ii,:) = rho(ii,:) - rho0;  
end

% ellispometry
lambda_Au = lambda_Au*10^-9; %m
nq = interp1(lambda_Au, n_Au,lambda);
kq = interp1(lambda_Au, k_Au,lambda);

eps_r = nq.^2-kq.^2;
eps_I = 2.*nq.*kq; % T_e = 300K
eps = eps_r + 1i*eps_I;
omg = (2*pi./lambda).*c;

for ii = 1: size(Te,2)
delepsI(ii,:) = (delrho(ii,:)./rho0).*eps_I;
delchiI(ii,:) = delepsI(ii,:)./eps0;
delchir(ii,:) = kkrebook2(omg,delchiI(ii,:),0); % using Kramers -Kronig to find the real part of change in permitivity
delepsr(ii,:) = delchir(ii,:).*eps0;
end
% figure;
% plot(lambda,delepsI); hold on
% plot(lambda,delepsr)


%% with angle
%lambda0 = linspace(0.725,0.9,200)*1e-6;
%Rq = interp1(lambda, R,lambdaq,'nearest','extrap');
count = 1;
for kk = 1: size(Te,2)
dR_R(:,count) = (1./R').*(p_re.* delepsr(kk,:) + p_Im.* delepsI(kk,:));
count = count + 1;
end
