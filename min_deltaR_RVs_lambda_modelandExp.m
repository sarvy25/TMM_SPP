clc
clear all
RR0 = [4.6,5,5.3,5.5,5.7,6,6.5,7,7.8,8.7,9.7,10.7,11.78,11.7,13,13.6,14.2,14.8,15.3,15.8,16.19]

load R730
load R740
load R745
load R750
load R755
load R760
load R765
load R770
load R775_2
load t


R730= R730- mean(R730(1:24));
R740= R740- mean(R740(1:24));
R745= R745- mean(R745(1:24));
R750= R750- mean(R750(1:24));
R755= R755- mean(R755(1:24));
R760= R760- mean(R760(1:24));
R765= R765- mean(R765(1:24));
R770= R770- mean(R770(1:24));
R775= R775_2- mean(R775_2(1:24));
p = 40;
f730= conv(R730./RR0(7),ones(p,1)/p,'same');
f740= conv(R740./RR0(9),ones(p,1)/p,'same');
f745= conv(R745./RR0(10),ones(p,1)/p,'same');
f750= conv(R750./RR0(10),ones(p,1)/p,'same');
f755= conv(R755./RR0(10),ones(p,1)/p,'same');
f760 = conv(R760./RR0(13),ones(p,1)/p,'same');
f765 = conv(R765./RR0(13),ones(p,1)/p,'same');
f770 = conv(R770./RR0(13),ones(p,1)/p,'same');
f775 = conv(R775./RR0(16),ones(p,1)/p,'same');

f1 =f730./abs(min(f745(:)));
f2 = f740./abs(min(f745(:)));
f3 = f745./abs(min(f745(:)));
f4 = f750./abs(min(f745(:)));
f5 = f755./abs(min(f745(:)));
f6 = f760./abs(min(f745(:)));
f7 = f765./abs(min(f745(:)));
f8 = f770./abs(min(f745(:)));
f9 = f775./abs(min(f745(:)));
figure; 
cc= hsv(5);
plot(t-1.1,f1,'color',cc(1,:),'linewidth',2); hold on
plot(t-1.1,f2,'color',cc(2,:),'linewidth',2); hold on
plot(t-1.1,f3 ,'color',cc(3,:),'linewidth',2); hold on
plot(t-1.1,f4,'color',cc(4,:),'linewidth',2); hold on
%plot(t-1,f5,'linewidth',1.5); hold on
plot(t-1.1,f6,'color',cc(5,:),'linewidth',2); hold on
%plot(t-1,f7,'linewidth',1.5); hold on
%plot(t-1,f8,'linewidth',1.5); hold on

ylabel('\DeltaR/R');xlabel('\Deltat(ps)'); axis tight
set(gca,'fontsize',23)
legend('\lambda = 730 nm','\lambda = 740 nm','\lambda = 745 nm','\lambda = 750 nm','\lambda = 760 nm')
xlim([-0.7 5]); ylim([-1 0])
daspect([2.3 1 1])

    %%  finding min of signal
    
    lambda_exp = [730, 740,745,750,755,760,765];
    %lambda_exp = [730, 740,745,750,755,760,765,775];

    mm = [min(f1), min(f2), min(f3),min(f4),min(f5),min(f6),min(f7)]
    figure; plot(lambda_exp, abs(mm),'p','MarkerFaceColor','m','MarkerSize',20);
    set(gca,'fontsize',25); xlabel('\lambda (nm)'); ylabel('norm {\DeltaR/R}')

    %% model 
    z = linspace(0,100*10^(-9),200); % 1 micron sample thickness
    tspan = linspace(0,5*10^(-12), 200);
    tau = 2*10^-12;
    [u1, u2]= TwoTempModel_Au_modified_tau(z,tspan,tau); % tau as the fitting paramete
    %% crop
    Te_crop = u1((1:end),(88:end));
    z_crop = z(88:end) - z(88);
    [m,n]= size(Te_crop);
    %% integration of temperature over depth 50 nm
    
for kk = 1: size(Te_crop,1)
    Te_avg(kk,1)= sum(Te_crop(kk,:))./size(Te_crop,2);
end
    theta = 43.58;
    lambda = linspace(0.55,0.9,150)*1e-6; %chnage it to 70 for seperate plot of lambdas
    [p_re, p_Im, REF]= myTMM_Partial_R_eps (lambda,theta);
    Te_avg = Te_crop(:,end);
    count = 1;
    for ii = 1:size(Te_avg,1)
    [dR_R(:,count)] = mydR_R_WithAngle (Te_avg(ii), p_re, p_Im, theta,REF,lambda);
    count = count+1;
    end
    zer = zeros(size(dR_R,1),10);
    dd = (dR_R)./abs(min(dR_R(:)));
    dd2 = [zer dd] ;
    
    %% add when you want to plot the minimums of each signal
    
      nn = 1
      for ii = 62:100
        mm(1,nn) = [min(dd(ii,:))];
        nn = nn + 1;
      end
    hold on;
    plot(lambda(62:100)*10^9, abs(mm)','k','linewidth',1.5); axis square
    xlim ([700 800])
    
    %% model on 730 740 745 750 765 nm wavelength
    
    figure; 
    plot(tspan(1:size(Te_avg,1))*10^12, dd(36,:),'color',cc(1,:),'linewidth',2); hold on %lambda=727
    plot(tspan(1:size(Te_avg,1))*10^12, dd(38,:),'color',cc(2,:),'linewidth',2); hold on %lambda = 737(118)
    plot(tspan(1:size(Te_avg,1))*10^12, dd(39,:),'color',cc(3,:),'linewidth',2); hold on% lambda = 742 (124)
    plot(tspan(1:size(Te_avg,1))*10^12, dd(42,:),'color',cc(4,:),'linewidth',2); hold on% lambda = 757(129)
    plot(tspan(1:size(Te_avg,1))*10^12, dd(43,:),'color',cc(5,:),'linewidth',2); hold on% lambda = 763 (131)
    ylabel('\DeltaR/R');xlabel('\Deltat(ps)'); axis tight
    set(gca,'fontsize',23);
    legend('\lambda = 730 nm','\lambda = 740 nm','\lambda = 745 nm','\lambda = 750 nm','\lambda = 760 nm')
    daspect([2.3 1 1])
%%
tstep =tspan(1:size(Te_avg,1));
dtstep = (tstep (2) - tstep(1))*10^12;
t2 = [(-10)*dtstep (-9)*dtstep (-8)*dtstep (-7)*dtstep (-6)*dtstep ...
    (-5)*dtstep (-4)*dtstep (-3)*dtstep (-2)*dtstep (-1)*dtstep tspan(1:size(Te_avg,1))*10^12]


    figure; 
    plot(t2, dd2(36,:),'color',cc(1,:),'linewidth',2); hold on %lambda=727
    plot(t2, dd2(38,:),'color',cc(2,:),'linewidth',2); hold on %lambda = 737(118)
    plot(t2, dd2(39,:),'color',cc(3,:),'linewidth',2); hold on% lambda = 742 (124)
    plot(t2, dd2(41,:),'color',cc(4,:),'linewidth',2); hold on% lambda = 757(129)
    plot(t2, dd2(42,:),'color',cc(5,:),'linewidth',2); hold on% lambda = 763 (131)
    ylabel('\DeltaR/R');xlabel('\Deltat(ps)'); axis tight
    set(gca,'fontsize',23);
    legend('\lambda = 730 nm','\lambda = 740 nm','\lambda = 745 nm','\lambda = 750 nm','\lambda = 760 nm')
    daspect([2.3 1 1]);ylim([-1 0])