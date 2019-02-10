clc;clear all;
load lambda_Au.mat % ellipsometry wavelength (nm)
load n_Au.mat % n of Au from ellipsometry measurement 
load k_Au.mat % k of Au from ellipsometry measurement 



%% First layer (glass)
% 
c = 3*10^8;
lambda = linspace(0.7,0.9,500)*1e-6;
k0 = 2*pi./lambda;
omg = c*2*pi./lambda;
THETA_ext_deg= linspace(40,50,500);  %(35:.1:55)';
THETA = THETA_ext_deg*(pi/180);





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

nq_Au = interp1(lambda_Au*1e-3',n_Au',lambda*1e6); % both lambda should have the same unit
kq_Au = interp1(lambda_Au*1e-3',k_Au',lambda*1e6);
er_Au = nq_Au.^2 - kq_Au.^2;
eI_Au = 2.*nq_Au.*kq_Au;

%%
% prism
nprism = 1.5;
lambda = linspace(0.7,0.9,500)*1e-6; % nm
en(1) = 1.5; %prism first layer
ek(1) = 0; 

    er= en(1)^2-ek(1)^2;
    ei= 2*en(1)*ek(1);
    e(1)=complex(er,ei);
% Au layer
% en(2) = 0.38797; 
% ek(2) = 8.7971; 
for  jj =  1:size(lambda,2)%250 % (jj = 40,lambda=609) (jj = 250 ,lambda=1198)%1:size(lambda,2)


    
    e(2) = complex(er_Au(jj),eI_Au(jj));
    d(2) = 52*1e-9;
    
%%%% air
    en(3) = 1;
    ek(3) = 0;    
    er= en(3)^2-ek(3)^2;
    ei= 2*en(3)*ek(3);
    e(3)=complex(er,ei);

    %%%%% START

    for jtheta=1:length(THETA); 
          theta=THETA(jtheta);

        q1=sqrt(e(1)-en(1)^2*sin(theta)^2)/e(1);
        qn=sqrt(e(end)-en(1)^2*sin(theta)^2)/e(end);
        for j= 2:(length(e)-1)
            beta(:,jj)=(d(j)*2*pi/lambda(jj))*sqrt(e(j)-en(1)^2*sin(theta)^2);
            q=sqrt(e(j)-en(1)^2*sin(theta)^2)/e(j);
            em(j,1,1)=cos(beta(:,jj));
            em(j,1,2)=-1i*sin(beta(:,jj))/q;
            em(j,2,1)=-1i*sin(beta(:,jj))*q;
            em(j,2,2)=cos(beta(:,jj));
        end
        emtot=[1 0;
            0 1];
        for j=2:(length(e)-1)
            emtot1(:,:)=em(j,:,:);
            emtot=emtot*emtot1;
        end

        rp =((emtot(1,1)+emtot(1,2)*qn)*q1-(emtot(2,1)+emtot(2,2)*qn))/...
            ((emtot(1,1)+emtot(1,2)*qn)*q1+(emtot(2,1)+emtot(2,2)*qn));
        tp=2*q1/((emtot(1,1)+emtot(1,2)*qn)*q1+(emtot(2,1)+emtot(2,2)*qn));
        ref = rp*conj(rp);
        tra=(tp*conj(tp)/cos(theta))*en(1)*real(qn);
        TRA(jj,jtheta)=tra;
        REF(jj,jtheta)= ref;
        
  
        
    end

end
   omg = 2*pi*c./lambda;
   ABS = ones(size(REF)) - REF-TRA;
   figure;
   imagesc(THETA_ext_deg,lambda*10^9, REF); xlabel('\theta'); ylabel('\lambda (nm)');view(2); colormap hot; axis tight
   set (gca, 'fontsize',25) 
   
   figure;
   imagesc(THETA_ext_deg,lambda*10^9, ABS ); xlabel('\theta'); ylabel('\lambda (nm)');view(2); colormap hot; axis tight
   set (gca, 'fontsize',25)
   
   %%
[~,ind1]= max(ABS(101,:))
[~,ind2]= min(REF(101,:))

figure;
plot(lambda*10^9, ABS(:,ind1),'k','linewidth',1.5); hold on
plot(lambda*10^9, REF(:,ind1),'r','linewidth',1.5); hold on

xlabel('\lambda(nm)'); ylabel('A/R'); grid on
set (gca, 'fontsize',25)


   