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

    %% Finding partial derivitive of reflectivity wrt real and imaginary of permitivity
    %theta = [43.4 43.6]% 43.8 44 44.2 44.4];
    theta = 43.6;
    lambda = linspace(0.55,0.9,200)*1e-6;
    %lambda = linspace(0.7,0.8,100)*1e-6;

    figure;
   % for kk = 1:size(theta,2)
    [p_re, p_Im, REF]= myTMM_Partial_R_eps (lambda,theta);
    plot(lambda,REF,'linewidth',1.5); hold on
    %end
    %%
    
    count  = 1; 
    dt = tspan(2)-tspan(1);
    Te_avg = Te_crop(:,end); %Te_crop(:,6);
    for ii = 1:size(Te_avg,1)
    [dR_R(:,count)] = mydR_R_WithAngle (Te_avg(ii), p_re, p_Im, theta,REF,lambda);
   % sample_temp(count,1) = Te_avg(ii)
    %axis([700*10^-9 1000*10^-9 -20*10^-3 1*10^-3])
    count = count+1;
    end
  
    %%
    
    dR2 = dR_R-max(dR_R(:));
    dd = (dR2)./abs(min(dR2(:)))
    
    figure;
    imagesc(tspan(2:size(Te_avg,1))*10^12,lambda*10^9, dd);colormap hot;
    caxis([-1 0])
    set ( gca, 'ydir', 'reverse' )
    xlabel('time(ps)'); ylabel('\lambda(nm)')
    set(gca,'fontsize',25)
    %% 
    clear count
    count = 1;
    for ii =90:1:131
        
        [m_dd(:,count),~]= min(dd(ii,:));
     count = count+1
    end
   figure; plot(lambda(90:1:131), m_dd,'-o'); hold on

    
    %%
    figure; 
    plot(tspan(1:size(Te_avg,1))*10^12, dd(109,:),'linewidth',1.5); hold on %lambda=745
    plot(tspan(1:size(Te_avg,1))*10^12, dd(112,:),'linewidth',1.5); hold on %lambda=745
    plot(tspan(1:size(Te_avg,1))*10^12, dd(118,:),'linewidth',1.5); hold on %lambda = 755
    plot(tspan(1:size(Te_avg,1))*10^12, dd(124,:),'linewidth',1.5); hold on% lambda = 766
    plot(tspan(1:size(Te_avg,1))*10^12, dd(129,:),'linewidth',1.5); hold on% lambda = 775
    plot(tspan(1:size(Te_avg,1))*10^12, dd(131,:),'linewidth',1.5); hold on% lambda = 785

    %%
    figure;
    [~,h] = contourf(tspan(2:size(Te_avg,1))*10^12,lambda*10^9, (dR2)./abs(min(dR2(:))),15);colormap hot;
    h.LineWidth = 1.7; caxis([-1 0])
    set ( gca, 'ydir', 'reverse' )
    xlabel('time(ps)'); ylabel('\lambda(nm)')
    set(gca,'fontsize',25)
    %%
    %[~, dR_R1] = mydR_R_0angle (E, Te);
    Te = [800,900,1000];
    figure; 
    for i = 1: size(Te,2)
    [p_re(i,:), p_Im(i,:),REF(i,:)]= myTMM_Partial_R_eps (lambda,theta);
    %[dR_R] = mydR_R_WithAngle (Te(i),p_re, p_Im, theta,REF,lambda);
    plot(p_Im(i,:)); hold on
    end