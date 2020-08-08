% ***************************************************************
% *** Matlab code for synthetic model-2 with varying density distribution.
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
%Final synthetic 2nd model for gravity field and inversion
clear all
close all

%%%%    Creation of Synthetic Model and finding Gravity anomaly     %%%%
%% function for Synthetic Example of depth of the basin
z_peaks=peaks(100);
y_peaks=((z_peaks(41,:))+1.5)*200;

%t and c are Legendre Gaussian quadrature points for numerical integration
    [t_leg,c_leg]=lgwt(10,0,1); 
    %synthetic depth and observation points
    x_obs=linspace(0,5000,100);
    depth=y_peaks;
    %plotting of x_obs vs. depth
    figure(1)
    subplot(2,1,2)
    plot(x_obs,depth)
    set(gca,'Ydir','reverse')
    xlabel('Distance in meter')
    ylabel('Depth in meter')
    title('Depth profile of the synthetic basin (Model 2)')
    box on
    %Finding Gravity field of the basin for varying density
    density=@(z)(-0.38-0.42*exp(-0.5*z*10^-3))*1000; 
    z_obs=0;   %height of observation point is in surface
    %Closed polygonic profile of the basin
    xx1=[x_obs x_obs(end) 0];
    yy1=[depth 0 0];
    %Gravity anomaly for fixed density model 
    zz1=poly_gravityrho(x_obs,z_obs,xx1,yy1,density,t_leg,c_leg);
    zz2=zz1;
    %adding noise to anomaly having 0 mean and sqrt(0.5) standard deviation 
         zz1 = zz1+sqrt(0.1).*randn(size(zz1))+0;
    %for model without noise comment line 46. 
    %Plotting the depth profile 
    figure(1)
    subplot(2,1,1)
    plot(x_obs,zz1,'o-')
    hold on
    plot(x_obs,zz2)
    xlabel('Observation points in meter')
    ylabel('Gravity anomaly in mGal')
    title('Gravity anomaly for Model 2 ')
    %Fourier coefficients for Gravity anomaly
    TT=1:2*length(zz1);
    ll=(TT(end)-TT(1))/2;
    %Gravity data for Fourier series 
    grav_data=[zz1 fliplr(zz1)];
    %Depth data for Fourier series 
    depth_data=[depth fliplr(depth)];
    %1st 100 Fourier Coefficients 
    for j=1:100
      %coefficients for Gravity field 
      ss1=grav_data.*cos(j*pi*TT/ll);  
      aa_grv(j)=(1/ll)*trapz(TT,ss1);
      ss2=grav_data.*sin(j*pi*TT/ll);  
      bb_grv(j)=(1/ll)*trapz(TT,ss2);
      vv_grv(j+1)=sqrt((aa_grv(j))^2+(bb_grv(j))^2);
      
      %coefficients for Depth Profile  
      dd1=depth_data.*cos(j*pi*TT/ll);  
      aa_dep(j)=(1/ll)*trapz(TT,dd1);
      dd2=grav_data.*sin(j*pi*TT/ll);  
      bb_dep(j)=(1/ll)*trapz(TT,dd2);
      vv_dep(j+1)=sqrt((aa_grv(j))^2+(bb_grv(j))^2);
      
    end
    %normalization of Fourier Coefficients of Gravity anomaly 
    vv_grv(1)=abs((1/ll)*trapz(TT,grav_data));
    vv_grv=vv_grv./max(vv_grv);
    %normalization of Fourier Coefficients of Depth Profile 
    vv_dep(1)=abs((1/ll)*trapz(TT,depth_data));
    vv_dep=vv_dep./max(vv_dep);
    
    figure(2)
    %Plotting the Fourier Spectrum of Gravity anomaly
    subplot(2,1,1)
    semilogy(vv_grv(1:51))
    title('Fourier Power spectrum for Gravity Anomaly')
    xlabel('Coefficients')
    ylabel('Amplitude')
    %Plotting the Fourier Spectrum of Depth profile
    subplot(2,1,2)
    semilogy(vv_dep(1:51))
    title('Fourier Power spectrum for Depth Anomaly')
    xlabel('Coefficients')
    ylabel('Amplitude')
%selecting number of Fourier Coefficients    
n=find(vv_grv>=10^-2, 1, 'last' )-1;
%printing the number of parameters 
fprintf('Number of parameters for Fourier coefficients of Depth profile = %d\n',n)
data_x=x_obs; %data for horizontal locations
data=zz1;     %data for gravity anomaly  
%correlation coefficients for Fourier spectrum of gravity anomaly and depth profile
cc=corrcoef(vv_grv,vv_dep);
fprintf('Correlation coefficients of power spectrum for gravity and depth profile is %f\n',cc(1,2))

%% creating Fourier Transformation matrix for multiplication 
T=1:2*length(data_x);
l=(T(end)-T(1))/2;
for i=1:length(T)
    jj1=0;jj2=0;
    for j=1:2*n+1
        if j==1
            A_mat(i,j)=(1/2);
        elseif j>1 && j<=n+1
            jj1=jj1+1;
            A_mat(i,j)=cos(jj1*pi*T(i)/l);
        else
            jj2=jj2+1;
            A_mat(i,j)=sin(jj2*pi*T(i)/l);
        end
    end
end
     
%% Problem Definition
%%%%    Optimization of Synthetic Model using PSO     %%%%
c1=1.4; c2=1.7; 
            %Cost function with constraints for optimization
            CostFunction =@(x) myCostFunction(x,data,A_mat,t_leg,c_leg)+1000*(constrained1(x,A_mat)+constrained2(x,A_mat)); 
            nVar=(2*n+1);          %Number of Unknown Variables
            MaxIt = 300;           %Maximum number of iterations
            nPoP =  40;            %Population size or swarm size
            %Loop for 10 independent run
            for ii=1:10
                %Calling the function for optimization using PSO
                [bst_var, best_cost,iter_count,error_energy] = WIPSO(CostFunction,nVar,MaxIt,nPoP,c1,c2);
                %Saving all best cost and parameters for each independent run 
                all_best(ii)=best_cost;
                all_var(:,ii)=bst_var;
                err_enrgy(:,ii)=error_energy';
            end
            %finding best model 
            [val,id]=min(all_best);
            %Parameters for best Model 
            best_var=(squeeze(all_var(:,id)))';  
            %error energy for best model 
            best_error=(squeeze(err_enrgy(:,id)))';
            %saving error energy data for plotting 
            save error_energy_mode2_varying_wo best_error -ascii
            
%% Plotting    
%%%%    Plotting the results for Optimization model    %%%%
       x=best_var'; %Fourier coefficients for optimized model
       xx=x_obs;    % observation points
     %gravity field at 40 discrete data points  
     xx1=linspace(min(xx),max(xx),40);
     %depth profile from optimized Fourier coefficients
     yy_eval=(A_mat*x)';
     yy=yy_eval(1:length(yy_eval)/2);
     yy1=spline(xx,yy,xx1);
     %close polygonal form of depth profile 
     x1(1:length(xx1)+2)=[xx1 max(xx1) min(xx1)];
     y1(1:length(yy1)+2)=[yy1 0 0];
     %density contrast
     density=@(z)(-0.38-0.42*exp(-0.5*z*10^-3))*1000; 
     %gravity field for optimized depth profile 
     zz_eval=poly_gravityrho(xx,0,x1,y1,density,t_leg,c_leg);
    
     %true model 
     figure(3)
     subplot(2,1,2)
     hold on
     xx_p=xx;
     xx_p=[0 xx_p xx(end)];
     depth_p=[0 depth 0];
     %patch plot for true model 
     zd_p=density(depth_p);
     patch(xx_p,depth_p,zd_p)
     colormap copper
     %optimized model 
     plot(xx1,yy1,'linewidth',1.25,'color','r')
     set(gca,'Ydir','reverse')
     c = colorbar('location','southoutside');
     c.Label.String = 'Density in kg/m^3';
     xlabel('Distance in meter')
     ylabel('Depth in meter')
     legend('True Model','Best Model','location','southeast')
     title('Depth profile of synthetic basin (Model 2)')
     box on
     subplot(2,1,1)
     hold on
     %Observed 
     plot(xx,data,'-o','linewidth',1.25)
     %Inverted 
     plot(xx,zz_eval,'linewidth',1.25,'color','r')
     %Axis lebeling 
     xlabel('Observation points in meter')
     ylabel('Gravity anomaly in mGal')
     title('Gravity anomaly for Model 2 ')
     legend('Observed Data','Inverted Data','location','southeast')
     box on
     
   %% Uncertainty Estimation 
  %RMSE for gravity  
     N_g=length(data); %Number of Observation points 
     N_d=length(depth); %Number of Observation points 
     %RMSE of given model 
     RMSE_g=sqrt((sum((data-zz_eval).^2))/N_g)/(max(data(:))-min(data(:)));
     %RMSE of True model 
     RMSE_true=sqrt((sum((data-zz2).^2))/N_g)/(max(data(:))-min(data(:)));
     % RMSE error for depth profile 
     RMSE_d=sqrt((sum((depth-yy).^2))/N_d)/(max(depth(:))-min(depth(:)));
 %Printing the RMSE error for depth and gravity profile 
     fprintf('RMSE in gravity field=%f, and in depth profile=%f\n',RMSE_g,RMSE_d)
     fprintf('True RMSE in gravity field=%f\n',RMSE_true)

     
%%%%    Objective functions and Constraints   %%%%
function val=myCostFunction(x,data,A_mat,t_leg,c_leg)
     %%inputs are 
     %x= parameters for PSO algorithm
     %data= observed gravity field 
     %A_mat= Matrix with Fourier coefficients 
     %t_leg and c_leg are  Legendre Gaussian quadrature points for numerical integration
     % subroutine for t_leg and c_leg evaluation is given in lgwt.m file 
    %observation points having 100 linearly spaced data points  
     xx=linspace(0,5*10^3,100);
    %40 linearly spaced data points in closed polygonal form
     xx1=linspace(min(xx),max(xx),40);
     yy1=(A_mat*x')';
     yy=yy1(1:length(yy1)/2);
     yy1=spline(xx,yy,xx1);
     %close polygonal form of depth profile 
     x1(1:length(xx1)+2)=[xx1 xx1(end) xx1(1)];
     y1(1:length(yy1)+2)=[yy1 0 0];
     %density contrast
     density=@(z)(-0.38-0.42*exp(-0.5*z*10^-3))*1000; 
     %gravity field for given depth profile 
     zz1=poly_gravityrho(xx,0,x1,y1,density,t_leg,c_leg);
     %misfit functional for observed and inverted gravity anomaly
     val=norm(data-zz1);
end

%Constraints for upper bound 
function val=constrained1(x,A_mat)
    %%inputs are 
     %x= parameters for PSO algorithm
     %A_mat= Matrix with Fourier coefficients
     %depth profile from Fourier domain to time domain
     yy1=A_mat*x';
     yy=yy1(1:length(yy1)/2);
     %minimum of depth 
     m_yy=min(yy);
     %penalty barrier method for minimum bound 
     gg=(-m_yy+0);
     val=(max(0,gg))^2;
     
end
%Constraints for lower bound 
function val=constrained2(x,A_mat)
     %%inputs are 
     %x= parameters for DE algorithm
     %A_mat= Matrix with Fourier coefficients
     %depth profile from Fourier domain to time domain
     yy1=A_mat*x';
     yy=yy1(1:length(yy1)/2);
     %maximum of depth 
     m_yy=max(yy);
     %penalty barrier method for minimum bound 
     gg=(m_yy-6000);
     val=(max(0,gg))^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  