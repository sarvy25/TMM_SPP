    clc
    clear all
    close all
    %% start by using TTM and cropping it accordingly

    z = linspace(0,100*10^(-9),200); % 1 micron sample thickness
    tspan = linspace(0,5*10^(-12), 200);
    tau = 2*10^-12;
    [u1, u2]= TwoTempModel_Au_modified_tau(z,tspan,tau); % tau as the fitting parameter
    figure ; mesh(z*10^9, tspan*10^12, u1 )
    xlabel('Depth z (nm)');ylabel('Time (ps)'); colormap hot
    set(gca,'fontsize',20); axis tight


    %% crop
    Te_crop = u1((1:end),(88:end));

    z_crop = z(88:end) - z(88);
    figure ; mesh(z_crop*10^9, tspan*10^12, Te_crop )
    xlabel('Depth z (nm)');ylabel('Time (ps)'); colormap hot
    set(gca,'fontsize',20); axis tight
    [m,n]= size(Te_crop);
    cc = jet(6);
    j = 1;
    figure
    for ii = 1:20:size(Te_crop,2)
    plot(tspan*10^12, Te_crop(:,ii),'color', cc(j,:),'linewidth',2); hold on
     j = j + 1;
    end
    xlabel('t(ps)'); ylabel('T_{e} (k)'); grid on
    set(gca,'fontsize',25);
    
    %% integration of temperature over depth 50 nm
for kk = 1: size(Te_crop,1)
    Te_avg(kk,1)= sum(Te_crop(kk,:))./size(Te_crop,2);
end
figure; plot(tspan*10^12, Te_avg,'k','linewidth',1.5)
xlabel('Delay time(ps)'); ylabel('Te(k)'); set(gca,'fontsize',25)

    %% Finding partial derivitive of reflectivity wrt real and imaginary of permitivity
   % theta = [43.4 43.6]% 43.8 44 44.2 44.4];
    theta = 43.61;
    lambda = linspace(0.55,0.9,70)*1e-6; %5 nm increment 
    figure;
    [p_re, p_Im, REF]= myTMM_Partial_R_eps (lambda,theta);
    plot(lambda,REF,'linewidth',1.5); hold on
    %%
    
    count  = 1; 
    dt = tspan(2)-tspan(1);
    Te_avg = Te_crop(:,end); %Te_crop(:,6);
    %Te_avg = [700,1000];
    for ii = 1:size(Te_avg,1)
    [dR_R(:,count)] = mydR_R_WithAngle (Te_avg(ii), p_re, p_Im, theta,REF,lambda);
    count = count+1;
    end
  
     % removing the offset mahyar added 
    %dR2 = dR_R-max(dR_R(:));
    % plotting the normalized data
    dR2 = dR_R;
    dd = (dR2)./abs(min(dR2(:)))
    figure;
    imagesc(tspan(1:size(Te_avg,1))*10^12,lambda*10^9, dd);colormap hot;
    caxis([-1 0])
    axis ([0 5 720 780]); 
    xlabel('time(ps)'); ylabel('\lambda(nm)')
    set(gca,'fontsize',25)
    %%
    
   % dR2 = dR_R-max(dR_R(:));
    dR2 = dR_R;
    dd = (dR2)./abs(min(dR2(:)))
    figure;
    [C,h] = contourf(tspan(1:size(Te_avg,1))*10^12,lambda*10^9, dd,20);colormap hot;
    h.LineWidth= 2;set ( gca, 'ydir', 'reverse' );caxis([-1 0])
    axis ([0 5 720 780]);
    xlabel('time(ps)'); ylabel('\lambda(nm)')
    set(gca,'fontsize',25)
  
    
    %%
    figure; 
    plot(tspan(1:size(Te_avg,1))*10^12, dd(36,:),'linewidth',1.5); hold on %lambda=727
    plot(tspan(1:size(Te_avg,1))*10^12, dd(38,:),'linewidth',1.5); hold on %lambda = 737(118)
    plot(tspan(1:size(Te_avg,1))*10^12, dd(39,:),'linewidth',1.5); hold on% lambda = 742 (124)
    plot(tspan(1:size(Te_avg,1))*10^12, dd(42,:),'linewidth',1.5); hold on% lambda = 757(129)
    plot(tspan(1:size(Te_avg,1))*10^12, dd(43,:),'linewidth',1.5); hold on% lambda = 763 (131)
    ylabel('\DeltaR/R');xlabel('\Deltat(ps)'); axis tight
    set(gca,'fontsize',23);
    %%  finding min of signal
    
    lambda_exp = [730, 740,745, 760, 775];
    mm = [min(dd(36,:)), min(dd(38,:)), min(dd(39,:)),min(dd(42,:)),min(dd(43,:))]
    figure; plot(lambda_exp, abs(mm),'p','MarkerFaceColor','m','MarkerSize',20);
    set(gca,'fontsize',25); xlabel('\lambda (nm)'); ylabel('norm {\DeltaR/R}')
