% ***************************************************************
% *** Matlab code for an arbritary model for finding correlation with depth and gravity profile.
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Dr. Chandra Prakash Dubey (email:p.dubey48@gmail.com)
% ***       Mr. M. Prasad (email:prasadgoud333@gmail.com)
% ***       Crustal Processes Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************
%Synthetic model for gravity field and Fourier transform
clear all
close all

%data for an arbritary synthetic basin 
x_obs=(importdata('x_obs_Model3.dat'))'; %observation points
depth= (importdata('depth_Model3.dat'))';%Depth profile
%t and c are Legendre Gaussian quadrature points for numerical integration
    [t_leg,c_leg]=lgwt(10,0,1); 
    %plotting of x_obs vs. depth
    figure(1)
    subplot(2,1,2)
    plot(x_obs,depth)
    set(gca,'Ydir','reverse')
    xlabel('Distance in meter')
    ylabel('Depth in meter')
    title('Depth profile of a complex synthetic Basin')
    box on
    %Finding Gravity field of the basin for exponential density in kg/m^3
    density=@(z)(-0.38-0.42*exp(-0.5*z*10^-3))*1000; 
    z_obs=0;   %Vertical position of observation point is in surface
    %Closed polygonic profile of the basin
    xx1=[x_obs x_obs(end) 0];
    yy1=[depth 0 0];
    %Gravity anomaly for fixed density model 
    zz1=poly_gravityrho(x_obs,z_obs,xx1,yy1,density,t_leg,c_leg);
    zz2=zz1;
    %Plotting the Gravity field anomaly 
    figure(1)
    subplot(2,1,1)
    plot(x_obs,zz2)
    xlabel('Observation points in meter')
    ylabel('Gravity anomaly in mGal')
    title('Gravity anomaly for complex synthetic Basin')
    %Fourier coefficients for Gravity anomaly
    TT=1:2*length(zz1);
    ll=(TT(end)-TT(1))/2;
    %Gravity data for Fourier series 
    grav_data=[zz1 fliplr(zz1)];
    %Depth data for Fourier series 
    depth_data=[depth fliplr(depth)];
    %1st 100 Fourier Coefficients 
    for j=1:100
      %Fourier coefficients for Gravity field 
      ss1=grav_data.*cos(j*pi*TT/ll);  
      aa_grv(j)=(1/ll)*trapz(TT,ss1);
      ss2=grav_data.*sin(j*pi*TT/ll);  
      bb_grv(j)=(1/ll)*trapz(TT,ss2);
      vv_grv(j+1)=sqrt((aa_grv(j))^2+(bb_grv(j))^2);
      
      %Fourier coefficients for Depth Profile  
      dd1=depth_data.*cos(j*pi*TT/ll);  
      aa_dep(j)=(1/ll)*trapz(TT,dd1);
      dd2=grav_data.*sin(j*pi*TT/ll);  
      bb_dep(j)=(1/ll)*trapz(TT,dd2);
      vv_dep(j+1)=sqrt((aa_grv(j))^2+(bb_grv(j))^2);
      
    end
    %Normalized amplitude of Power spectrum for Gravity anomaly
    vv_grv(1)=abs((1/ll)*trapz(TT,grav_data));
    vv_grv=vv_grv./max(vv_grv);
    %Normalized amplitude of Power spectrum for Depth Profile 
    vv_dep(1)=abs((1/ll)*trapz(TT,depth_data));
    vv_dep=vv_dep./max(vv_dep);
    %Plotting the power spectrum for Gravity anomaly
    figure(2)
    subplot(2,1,1)
    semilogy(vv_grv(1:51))
    title('Fourier Power spectrum for Gravity Anomaly')
    xlabel('Coefficients')
    ylabel('Amplitude')
    grid on
    %Plotting the power spectrum for Derpth Profile 
    subplot(2,1,2)
    semilogy(vv_dep(1:51))
    grid on
    title('Fourier Power spectrum for Depth Anomaly')
    xlabel('Coefficients')
    ylabel('Amplitude')
    
%finding correlation coefficients between power spectrum for Gravity
%anomaly and Depth profile 
cc=corrcoef(vv_grv,vv_dep);
fprintf('Correlation coefficients of power spectrum for gravity and depth profile is %f\n',cc(1,2))