% this function gives temparature of electrons(u1) and phonons (u2) 
function [u1, u2]= TwoTempModel_Au_modified_tau(z,tspan, tau)
mahyarhandler = @(x,t,u,DuDx) mahyarpdex(x, t, u, DuDx, tau);
m = 0;
sol = pdepe(m,mahyarhandler,@pdex4ic,@pdex4bc,z ,tspan);
u1 = sol(:,:,1);
u2 = sol(:,:,2);

function [c,f,s] = mahyarpdex(x,t,u,DuDx,tau)
    % constant 
Ae = 70;%95;%70;   %%%J m-3 k-2
phai = 0.2; % J/m^2
beta = 4 * log(2);
tp = 150*1e-15;%0.49*1e-12; %0.1*1e-12; % pulse duration 100 fsec 
R = 0.14; %
N = 5.9 * 10^28 ; %m^-1
kb = 1.38*10^-23; % Jk-1
Tf = 6.4*10^4; % K
Ke0 = 315;%415; % Wm-1 K-1

   %% definitng electron and lattice heat capacity
Ce = (pi^2/2).* N .*kb.* (u(1)./Tf) ; 
CL = 2.5*1e6;%Jm-3K-1,5.7*1e6; %2.5*1e6; % Lattice heat capacity Jm-3K-1
c = [ Ce;CL]; % electron and lattece heat capacity
   %% defining electron phonon coupling constant 
%g = Ce./(5*10^(-12));

%g =    3.6*1e16;%2.1,7.5*1e16; %2.6*1e16  Wm-3K-1 
%CL = 2.5*1e6;%Jm-3K-1,5.7*1e6; %2.5*1e6; % Lattice heat capacity Jm-3K-1
%c = [ u(1)*Ae;CL]; % electron and lattece heat capacity



f = [  Ke0*u(1)./u(2); 0] .* DuDx; 
y = u(1) - u(2);

% at lambda = 730 nm 
x2 = 50*10^-9; %offset 
a2 = 0.9846;
a1 = 0.1572;
b2 = 0.01547*10^-6;
b1 = 0.02013*10^-6;
dels = 15*10^-9;

% at lambda = 735 nm 
% x2 = 50*10^-9; %offset 
% a2 = 0.9718;
% a1 = 0.1778;
% b2 = 0.01557*10^-6;
% b1 = 0.02073*10^-6;
% dels = 25*10^-9;

% at lambda = 740 nm 
% x2 = 50*10^-9; %offset 
% a2 = 0.9781;
% a1 = 0.2026;
% b2 = 0.01566*10^-6;
% b1 = 0.02106*10^-6;
% dels = 15*10^-9;


% at lambda = 745 nm 
%  x2 = 50*10^-9; %offset 
%  a2 = 0.9487;
%  a1 = 0.2211;
%  b2 = 0.01576*10^-6;
%  b1 = 0.0212*10^-6;
%  dels = 15*10^-9;

 % at lambda = 750 nm 
%  x2 = 50*10^-9; %offset 
%  a2 = 0.9203;
%  a1 = 0.24;
%  b2 = 0.01585*10^-6;
%  b1 = 0.02121*10^-6;
%  dels = 15*10^-9;



 % at lambda = 755 nm 
%  x2 = 50*10^-9; %offset 
%  a2 = 0.8929;
%  a1 = 0.259;
%  b2 = 0.01595*10^-6;
%  b1 = 0.02112*10^-6;
%  dels = 15*10^-9;

 
% at lambda = 760 nm 
% x2 = 50*10^-9; %offset 
% a2 = 0.8813;
% a1 = 0.2841;
% b2 = 0.01607*10^-6;
% b1 = 0.02095*10^-6;
% dels = 15*10^-9;



% at lambda = 765 nm 
% x2 = 50*10^-9; %offset 
% a2 = 0.8555;
% a1 = 0.3034;
% b2 = 0.01617*10^-6;
% b1 = 0.02076*10^-6;
% dels = 15*10^-9;


% at lambda = 770 nm 
% x2 = 50*10^-9; %offset 
% a2 = 0.8307;
% a1 = 0.3226;
% b2 = 0.01628*10^-6;
% b1 = 0.02053*10^-6;
% dels = 15*10^-9;

% at lambda = 775 nm 
% x2 = 50*10^-9; %offset 
% a2 = 0.8078;
% a1 = 0.341;
% b2 = 0.01639*10^-6;
% b1 = 0.02031*10^-6;
% dels = 15*10^-9;

% at lambda = 780 nm 
% x2 = 50*10^-9; %offset 
% a2 = 0.7849;
% a1 = 0.3598;
% b2 = 0.0165*10^-6;
% b1 = 0.02007*10^-6;
% dels = 15*10^-9;

% at lambda = 785 nm 
% x2 = 50*10^-9; %offset 
% a2 = 0.763;
% a1 = 0.3783;
% b2 = 0.01662*10^-6;
% b1 = 0.01983*10^-6;
% dels = 15*10^-9;

% at lambda = 790 nm 
% x2 = 50*10^-9; %offset 
% a2 = 0.7277;
% a1 = 0.39;
% b2 = 0.01674*10^-6;
% b1 = 0.01958*10^-6;
% dels = 15*10^-9;

% at lambda = 795 nm 
% x2 = 50*10^-9; %offset 
% a2 = 0.7079;
% a1 = 0.4075;
% b2 = 0.01686*10^-6;
% b1 = 0.01935*10^-6;
% dels = 15*10^-9;

Source = sqrt(beta/pi)*phai*(a2.*exp((x-x2)./b2)+a1.*exp(-x./b1)).*exp(-beta*(t/tp).^2)*((1-R)/(tp*dels));
%s = [Source - g*y; g*y]; 
s = [Source - (Ce./(tau))*y; (Ce./(tau))*y]; 


 function u0 = pdex4ic(x);

u0 = [300; 300]; 




function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
pl = [ul(1); ul(2)];  % l corresponds to the z= 0 and R corresponds to z = L
ql = [1; 1]; 
pr = [ur(1); ur(2)]; %*0 added  by sarvy 
qr = [1; 1];  



