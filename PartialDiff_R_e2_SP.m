clc;clear all;
load lambda_Au.mat % ellipsometry wavelength (nm)
load n_Au.mat % n of Au from ellipsometry measurement 
load k_Au.mat % k of Au from ellipsometry measurement 



%% First layer (glass)
% 
c = 3*10^8;
lambdaq = linspace(0.725,0.9,200)*1e-6;
k0 = 2*pi./lambdaq;
omg = c*2*pi./lambdaq;
%THETA_ext_deg= linspace(40,50,500);  %(35:.1:55)';
THETA_ext_deg = 43.83;%
theta = THETA_ext_deg*(pi/180);





%% GOLD
%%%% q = 1.6*1e-19;
%%%% omgp = 5.8*q; % bulk plasma frequency of gold
%%%% omg = sqrt(omgp.^2+c.^2*k.^2);
% Gold epsilon before interpolation

[lambda_Au, ind_l] = unique(lambda_Au);
n_Au = n_Au(ind_l);
k_Au = k_Au(ind_l);
er1 = n_Au.^2 - k_Au.^2;
eI1 = 2.*n_Au.*k_Au;
% Gold epsilon After interpolation

nq_Au = interp1(lambda_Au*1e-3',n_Au',lambdaq*1e6); % both lambda should have the same unit
kq_Au = interp1(lambda_Au*1e-3',k_Au',lambdaq*1e6);
er_Au = nq_Au.^2 - kq_Au.^2;
eI_Au = 2.*nq_Au.*kq_Au;
%p1_re = zeros(1, 10);
%p1_Im = zeros(1, 10);

%%
% prism
nprism = 1.5;
%lambda = linspace(0.73,0.8,200)*1e-6; % nm
en(1) = 1.5; %prism first layer
ek(1) = 0;

er= en(1)^2-ek(1)^2;
ei= 2*en(1)*ek(1);
e(1)=complex(er,ei);
syms er2
syms eI2
for  jj =  1:size(lambdaq,2)%250 % (jj = 40,lambda=609) (jj = 250 ,lambda=1198)%1:size(lambda,2)
    
    
    
    %edata = complex(er_Au(jj),eI_Au(jj));
   %syms e2


    d(2) = 52*1e-9;
    
%%%% air
    en(3) = 1;
    ek(3) = 0;    
    er= en(3)^2-ek(3)^2;
    ei= 2*en(3)*ek(3);
    e(3)=complex(er,ei);

    %%%%% START


        q1=sqrt(e(1)-en(1)^2*sin(theta)^2)/e(1);
        qn=sqrt(e(end)-en(1)^2*sin(theta)^2)/e(end);
        
            beta(:,jj)=(d(2)*2*pi/lambdaq(jj))*sqrt((er2+1i*eI2)-en(1)^2*sin(theta)^2);
            q=sqrt((er2+1i*eI2)-en(1)^2*sin(theta)^2)/(er2+1i*eI2);
            em(1,1)=cos(beta(:,jj));
            em(1,2)=-1i*sin(beta(:,jj))/q;
            em(2,1)=-1i*sin(beta(:,jj))*q;
            em(2,2)=cos(beta(:,jj));
  
        emtot=[1 0;
            0 1];
        
            emtot=emtot*em;
        

        rp =((emtot(1,1)+emtot(1,2)*qn)*q1-(emtot(2,1)+emtot(2,2)*qn))/...
            ((emtot(1,1)+emtot(1,2)*qn)*q1+(emtot(2,1)+emtot(2,2)*qn));
        tp=2*q1/((emtot(1,1)+emtot(1,2)*qn)*q1+(emtot(2,1)+emtot(2,2)*qn));
        ref = rp*conj(rp);
        tra=(tp*conj(tp)/cos(theta))*en(1)*real(qn);
        TRA(jj,1)=tra;
        REF(jj,1)= eval(subs(ref, {er2,eI2}, {er_Au(1,jj),eI_Au(1,jj)}));
        
        p_er = diff(ref, er2); %taking partial derivitive wrt real permitivity (symbolic)
        p_eI = diff(ref, eI2); %taking partial derivitive wrt imaginary permitivity (symbolic)
        p1_re(1,jj) = eval(subs(p_er,{er2,eI2}, {er_Au(1,jj),eI_Au(1,jj)}));
        p1_Im(1,jj) = eval(subs(p_eI,{er2,eI2}, {er_Au(1,jj),eI_Au(1,jj)}));
        fprintf('%d \n', jj);
      %p1 = eval(subs(p_eI,{er}, {edata}));

     %  p2(:,kk) = eval(subs(p_REF,{eps1, eps2}, {eps_r(:,kk), eps_I(:,kk)}));
        
    

end
figure
plot(lambdaq,p1_Im,'r'); hold on;plot(lambdaq,p1_re,'y')
figure
plot(lambdaq*10^9, REF,'k','linewidth',1.5); hold on

   