% Spacecraft Guidance and Navigation (2023/2024)
% Assignment # 2, Exercise 1
% Author: Samuele Vincenzi
tic
clear
clc
close all
format long g

% CODE INIZIALIZATION
% Load spice kernels
cspice_furnsh('assignment02.tm');

% PROBLEM INITIAL DATA
% Separation epoch
t_sep_utc = '2010-08-12T05:27:39.114';
t_sep_et =  cspice_str2et(t_sep_utc);

% Satellite ID
Sat1_ID = '36599';
Sat2_ID = '36827';
sat_name = {'Mango','Tango'};

% Earth gravitational constant
GM = cspice_bodvrd('Earth','GM',1); % [km^3/s^2]

% Mean values of the state estimate at t_sep
r0_mean_sat1 = [4622.232026629, 5399.3369588058, -0.0212138165769957]'; % [km]
v0_mean_sat1 = [0.812221125483763, -0.721512914578826, 7.42665302729053]'; % [km/s]
x0_mean_sat1 = [r0_mean_sat1;v0_mean_sat1];
r0_mean_sat2 = [4621.69343340281, 5399.26386352847, -3.09039248714313]'; % [km]
v0_mean_sat2 = [0.813960847513811, -0.719449862738607, 7.42706066911294]'; % [km/s]
x0_mean_sat2 = [r0_mean_sat2;v0_mean_sat2];

% Rotation matrix from ECI to perifocal frame
Rot_eci2pf = rotation_eci2pf(r0_mean_sat1,v0_mean_sat1,GM);

% Covariance matrix of the state estimate at t_sep
P0 = [+5.6e-7 +3.5e-7 -7.1e-8 0        0        0;
      +3.5e-7 +9.7e-7 +7.6e-8 0        0        0;
      -7.1e-8 +7.6e-8 +8.1e-8 0        0        0;
      0       0       0       +2.8e-11 0        0;
      0       0       0       0        +2.7e-11 0
      0       0       0       0        0        +9.6e-12]; % [km2, km2/s, km2/s2]

%% POINT 1
% INITIAL DATA FOR POINT 1
N = 10; % Grid points
a1 = 1/(-(norm(v0_mean_sat1))^2/GM + 2/norm(r0_mean_sat1)); % [km] semi-major axis
T1 = 2*pi*sqrt(a1^3/GM); % [s] % Orbital period of s/c 1
t_g = linspace(t_sep_et,(t_sep_et + N*T1),N+1); % [s] Time grid
len = length(t_g);

% Struct definition and inizialization
% mean --- struct for mean value
% cov --- struct for covariance 
  % [].sat1 ---- mean value for satellite 1 propagation
  % [].sat2 ---- mean value for satellite 2 propagation
    % [].[].lincov --- uncertainty propagation with lincov
    % [].[].UT --- uncertainty propagation with unscented transform
    % [].[].MC --- uncertainty propagation with Montecarlo simulation   
mean.sat1.lincov = zeros(6,N+1); % [6x11] each column is one revolution
mean.sat1.UT = zeros(6,N+1);
mean.sat1.MC = zeros(6,N+1);
mean.sat2.lincov = zeros(6,N+1);
mean.sat2.UT = zeros(6,N+1);
mean.sat2.MC = zeros(6,N+1); 
cov.sat1.lincov = cell(1,N+1);  % {11} each element is one revolution
cov.sat1.UT = cell(1,N+1);
cov.sat1.MC = cell(1,N+1);
cov.sat2.lincov = cell(1,N+1);
cov.sat2.UT = cell(1,N+1);
cov.sat2.MC = cell(1,N+1);

% LINCOV APPROACH
tic
% Inizialization
% The first elements are the initial mean and covariance
mean.sat1.lincov(:,1) = x0_mean_sat1;
mean.sat2.lincov(:,1) = x0_mean_sat2;
cov.sat1.lincov{1} = P0;
cov.sat2.lincov{1} = P0;

for i = 2 : len
    t_in = t_g(i-1); % [s] Initial time for i-th propagation
    t_fin = t_g(i); % [s] Final time for i-th propagation
    
    % Uncertainty propagation for satellite 1
    x_star = mean.sat1.lincov(:,i-1); % Initial mean value
    [xf_sat1,PHIf_sat1,~,~,~]  = propagate([t_in t_fin],x_star,'Earth'); % 2BP propagation with STM
    mean.sat1.lincov(:,i) = xf_sat1; % Update of the mean
    cov.sat1.lincov{i} = PHIf_sat1*cov.sat1.lincov{i-1}*PHIf_sat1'; % Update of the covariance
    
    % Uncertainty propagation for satellite 2
    x_star = mean.sat2.lincov(:,i-1); % Initial mean value
    [xf_sat2,PHIf_sat2,~,~,~]  = propagate([t_in t_fin],x_star,'Earth'); % 2BP propagation with STM
    mean.sat2.lincov(:,i) = xf_sat2; % Update of the mean
    cov.sat2.lincov{i} = PHIf_sat2*cov.sat2.lincov{i-1}*PHIf_sat2'; % Update of the covariance
end
t_lincov = toc;

% UNSCENTED TRANSFORM APPROACH
tic
% Inizialization
alpha = 0.1;
beta = 2;
k = 0;
par = [alpha,beta,k];
 
% Uncertainty propagation for satellite 1
[y_mean,Py] = UT(par,x0_mean_sat1,P0,t_g); % UT solution
mean.sat1.UT = y_mean; % Update of the mean
cov.sat1.UT = Py; % Update of the covariance

% Uncertainty propagation for satellite 2
[y_mean,Py] = UT(par,x0_mean_sat2,P0,t_g); % UT solution
mean.sat2.UT = y_mean; % Update of the mean
cov.sat2.UT = Py; % Update of the covariance

t_UT = toc;
%% POINT 2
% CRITICAL CONDITION EVALUATION
% LINCOV APPROACH
N_rev = (t_g-t_g(1))/T1;
m = 0;
delta_r = zeros(1,len);
sigma3_max = zeros(1,len);
for i = 1 : len
    P_sum = cov.sat1.lincov{i}(1:3,1:3) + cov.sat2.lincov{i}(1:3,1:3); % delta_r
    delta_r(i) = norm(mean.sat1.lincov(1:3,i) - mean.sat2.lincov(1:3,i));
    sigma3_max(i) = 3*sqrt(max(eig(P_sum))); % 3sigma
    res = delta_r(i) - sigma3_max(i);
    if res < 0 && m==0
        res_Lincov = res;
        m = m+1;
        Nc_lincov = (t_g(i) - t_g(1))/T1;
    end
end

% PLOT OF THE CRITICAL CONDITION
figure
subplot(1,2,1)
hold on
grid on
plot(N_rev,delta_r,'Marker','o','LineStyle','-','LineWidth',1.5)
plot(N_rev,sigma3_max,'Marker','o','LineStyle','-','LineWidth',1.5)
title('Critical condition analysis : Lincov','Interpreter','latex')
legend('$\Delta r$','$3\sqrt{max(eig(P_{sum}))}$','Interpreter','latex')
xlabel('$N_{rev} [-]$','Interpreter','latex');
ylabel('$\Delta r,3\sigma_{max}[km]$','Interpreter','latex')
set(gca,'FontSize',22)

% UT APPROACH
m = 0;
delta_r = zeros(1,len);
sigma3_max = zeros(1,len);
for i = 1 : len
    P_sum = cov.sat1.UT{i}(1:3,1:3) + cov.sat2.UT{i}(1:3,1:3);
    delta_r(i) = norm(mean.sat1.UT(1:3,i) - mean.sat2.UT(1:3,i)); % delta_r
    sigma3_max(i) = 3*sqrt(max(eig(P_sum))); % 3sigma
    res = delta_r(i) - sigma3_max(i);
    if res < 0 && m==0
        res_UT = res;
        m = m+1;
        Nc_UT = (t_g(i) - t_g(1))/T1;
    end
end

% PLOT OF THE CRITICAL CONDITION
subplot(1,2,2)
hold on
grid on
plot(N_rev,delta_r,'Marker','o','LineStyle','-','LineWidth',1.5)
plot(N_rev,sigma3_max,'Marker','o','LineStyle','-','LineWidth',1.5)
title('Critical condition analysis : UT','Interpreter','latex')
legend('$\Delta r$','$3\sqrt{max(eig(P_{sum}))}$','Interpreter','latex')
xlabel('$N_{rev} [-]$','Interpreter','latex');
ylabel('$\Delta r,3\sigma_{max}[km]$','Interpreter','latex')
set(gca,'FontSize',22)

%% POINT 3
% MONTECARLO SIMULATION APPROACH
tic
% Inizialization
N_MC = 100;

% SAMPLE POINTS COMPUTATION
% Uncertainty propagation for satellite 1
[y_mean,Py,Y1] = MC(N_MC,x0_mean_sat1,P0,t_g); % MC solution
mean.sat1.MC = y_mean; % Update of the mean
cov.sat1.MC = Py; % Update of the covariance

% Uncertainty propagation for satellite 2
[y_mean,Py,Y2] = MC(N_MC,x0_mean_sat2,P0,t_g); % MC solution
mean.sat2.MC = y_mean; % Update of the mean
cov.sat2.MC = Py; % Update of the covariance

t_MC = toc;

% CRITICAL CONDITION EVALUATION WITH MC
res = 1;
i = 0;
while res>=0 && i<11
    i = i + 1;
    P_sum = cov.sat1.MC{i}(1:3,1:3) + cov.sat2.MC{i}(1:3,1:3);
    delta_r = norm(mean.sat1.MC(1:3,i) - mean.sat2.MC(1:3,i));
    res = delta_r - 3*sqrt(max(eig(P_sum)));
end
Nc_MC = (t_g(i) - t_g(1))/T1;

% COMPUTATION OF 3*sqrt(max(eig(P_{r,i})) AND 3*sqrt(max(eig(P_{v,i}))
sigma3.lincov = zeros(4,len);
sigma3.UT = zeros(4,len);
sigma3.MC = zeros(4,len);
for i = 1 : len
    % Lincov
    sigma3.lincov(1,i) = 3*sqrt(max(eig(cov.sat1.lincov{i}(1:3,1:3)))); % P_{r,1}
    sigma3.lincov(2,i) = 3*sqrt(max(eig(cov.sat1.lincov{i}(4:6,4:6)))); % P_{v,1}
    sigma3.lincov(3,i) = 3*sqrt(max(eig(cov.sat2.lincov{i}(1:3,1:3)))); % P_{r,2}
    sigma3.lincov(4,i) = 3*sqrt(max(eig(cov.sat2.lincov{i}(4:6,4:6)))); % P_{v,2}
    % Unscented transform
    sigma3.UT(1,i) = 3*sqrt(max(eig(cov.sat1.UT{i}(1:3,1:3)))); % P_{r,1}
    sigma3.UT(2,i) = 3*sqrt(max(eig(cov.sat1.UT{i}(4:6,4:6)))); % P_{v,1}
    sigma3.UT(3,i) = 3*sqrt(max(eig(cov.sat2.UT{i}(1:3,1:3)))); % P_{r,2}
    sigma3.UT(4,i) = 3*sqrt(max(eig(cov.sat2.UT{i}(4:6,4:6)))); % P_{v,2}
    % Montecarlo simulation
    sigma3.MC(1,i) = 3*sqrt(max(eig(cov.sat1.MC{i}(1:3,1:3)))); % P_{r,1}
    sigma3.MC(2,i) = 3*sqrt(max(eig(cov.sat1.MC{i}(4:6,4:6)))); % P_{v,1}
    sigma3.MC(3,i) = 3*sqrt(max(eig(cov.sat2.MC{i}(1:3,1:3)))); % P_{r,2}
    sigma3.MC(4,i) = 3*sqrt(max(eig(cov.sat2.MC{i}(4:6,4:6)))); % P_{v,2}
end

% PLOT OF 3*sqrt(max(eig(P_{r,i})) AND 3*sqrt(max(eig(P_{v,i}))
figure
title_name = {'Position-Mango','Velocity-Mango','Position-Tango','Velocity-Tango'};
y_label = {'[km]','[km/s]','[km]','[km/s]'};
for i = 1 : 4
    subplot(2,2,i)
    plot(N_rev,sigma3.lincov(i,:),'*','MarkerSize',13,'LineWidth',1.5,'MarkerEdgeColor',"#4DBEEE");
    hold on
    plot(N_rev,sigma3.UT(i,:),'o','MarkerSize',13,'LineWidth',1.5,'MarkerEdgeColor',"#EDB120");
    plot(N_rev,sigma3.MC(i,:),'+','MarkerSize',13,'LineWidth',1.5,'MarkerEdgeColor',"#7E2F8E");
    title(num2str(title_name{i}),'Interpreter','latex')
    grid on
    legend('LinCov','UT','MC','Interpreter','latex','Location','best');
    xlabel('$N_{rev} [-]$','Interpreter','latex');
    ylabel(['$3\sqrt{max(eig(P))}\quad$',num2str(y_label{i})],'Interpreter','latex')
    set(gca,'FontSize',18)
end

% PLOT OF MEAN VALUES TIME EVOLUTION
figure
subplot(1,2,1)
plot(N_rev,vecnorm(mean.sat1.lincov(1:3,:)),'*','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor',"#4DBEEE");
hold on
plot(N_rev,vecnorm(mean.sat1.UT(1:3,:)),'o','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor',"#EDB120");
plot(N_rev,vecnorm(mean.sat1.MC(1:3,:)),'+','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor',"#7E2F8E");
plot(N_rev,vecnorm(mean.sat2.lincov(1:3,:)),'*','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor',"#EDB120");
plot(N_rev,vecnorm(mean.sat2.UT(1:3,:)),'o','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor',"#77AC30");
plot(N_rev,vecnorm(mean.sat2.MC(1:3,:)),'+','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor',"#A2142F");
title('Time evolution of the mean value (position)','Interpreter','latex')
grid on
legend('LinCov-Mango','UT-Mango','MC-Mango','LinCov-Tango','UT-Tango','MC-Tango','Interpreter','latex','Location','best');
xlabel('$N_{rev} [-]$','Interpreter','latex');
ylabel('$||\hat{r}||\quad [km]$','Interpreter','latex')
ylim([7107 7107.8]);
set(gca,'FontSize',18)
subplot(1,2,2)
plot(N_rev,vecnorm(mean.sat1.lincov(4:6,:)),'*','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor',"#4DBEEE");
hold on
plot(N_rev,vecnorm(mean.sat1.UT(4:6,:)),'o','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor',"#EDB120");
plot(N_rev,vecnorm(mean.sat1.MC(4:6,:)),'+','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor',"#7E2F8E");
plot(N_rev,vecnorm(mean.sat2.lincov(4:6,:)),'*','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor',"#EDB120");
plot(N_rev,vecnorm(mean.sat2.UT(4:6,:)),'o','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor',"#77AC30");
plot(N_rev,vecnorm(mean.sat2.MC(4:6,:)),'+','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor',"#A2142F");
title('Time evolution of the mean value (velocity)','Interpreter','latex')
grid on
legend('LinCov-Mango','UT-Mango','MC-Mango','LinCov-Tango','UT-Tango','MC-Tango','Interpreter','latex','Location','best');
xlabel('$N_{rev} [-]$','Interpreter','latex');
ylabel('$||\hat{v}||\quad [km/s]$','Interpreter','latex')
%ylim([7107 7107.8]);
set(gca,'FontSize',18)

% MC SAMPLES + LINCOV,UT,MC MEAN AND COVARIANCE ROTATION ON ORBITAL PLANE
% Satellite 1
index = 11; % Change the index to plot for different revolutions
a = 4; % a = 1 for position, a = 4 for velocity
b = 6; % b = 3 for position, b = 6 for velocity
Y_2D{1} = Rot_eci2pf*Y1{index}(a:b,:); 
mean_lincov_2D{1} = Rot_eci2pf*mean.sat1.lincov(a:b,index);
mean_UT_2D{1} = Rot_eci2pf*mean.sat1.UT(a:b,index);
mean_MC_2D{1} = Rot_eci2pf*mean.sat1.MC(a:b,index);
ellipsoid_lincov_2D{1} = Rot_eci2pf*ellipsoid_computation(cov.sat1.lincov{index}(a:b,a:b),mean.sat1.lincov(a:b,index));
ellipsoid_lincov_2D{1} = ellipsoid_lincov_2D{1}(1:2,:)';
ellipsoid_UT_2D{1} = Rot_eci2pf*ellipsoid_computation(cov.sat1.UT{index}(a:b,a:b),mean.sat1.UT(a:b,index));
ellipsoid_UT_2D{1} = ellipsoid_UT_2D{1}(1:2,:)';
ellipsoid_MC_2D{1} = Rot_eci2pf*ellipsoid_computation(cov.sat1.MC{index}(a:b,a:b),mean.sat1.MC(a:b,index));
ellipsoid_MC_2D{1} = ellipsoid_MC_2D{1}(1:2,:)';
% Satellite 2
Y_2D{2} = Rot_eci2pf*Y2{index}(a:b,:); 
mean_lincov_2D{2} = Rot_eci2pf*mean.sat2.lincov(a:b,index);
mean_UT_2D{2} = Rot_eci2pf*mean.sat2.UT(a:b,index);
mean_MC_2D{2} = Rot_eci2pf*mean.sat2.MC(a:b,index);
ellipsoid_lincov_2D{2} = Rot_eci2pf*ellipsoid_computation(cov.sat2.lincov{index}(a:b,a:b),mean.sat2.lincov(a:b,index));
ellipsoid_lincov_2D{2} = ellipsoid_lincov_2D{2}(1:2,:)';
ellipsoid_UT_2D{2} = Rot_eci2pf*ellipsoid_computation(cov.sat2.UT{index}(a:b,a:b),mean.sat2.UT(a:b,index));
ellipsoid_UT_2D{2} = ellipsoid_UT_2D{2}(1:2,:)';
ellipsoid_MC_2D{2} = Rot_eci2pf*ellipsoid_computation(cov.sat2.MC{index}(a:b,a:b),mean.sat2.MC(a:b,index));
ellipsoid_MC_2D{2} = ellipsoid_MC_2D{2}(1:2,:)';

% PLOT OF MC SAMPLES +  LINCOV,UT,MC MEAN AND COVARIANCE
% Satellite 1 and 2
for j = 1 : 2
    figure
    for k = 1 : 3
        subplot(1,3,k)
        hold on
        switch k
            case 1
            plot(ellipsoid_lincov_2D{j}(:,1),ellipsoid_lincov_2D{j}(:,2),'-','Color',"#77AC30",'LineWidth',2);
            plot(mean_lincov_2D{j}(1),mean_lincov_2D{j}(2),'o','LineWidth',2,'MarkerSize',13,'Color','k','MarkerFaceColor',	"#7E2F8E");
            case 2
            plot(ellipsoid_UT_2D{j}(:,1),ellipsoid_UT_2D{j}(:,2),'-','Color',"#D95319",'LineWidth',2);
            plot(mean_UT_2D{j}(1),mean_UT_2D{j}(2),'o','LineWidth',2,'MarkerSize',13,'Color','k','MarkerFaceColor','g');  
            case 3
            plot(ellipsoid_MC_2D{j}(:,1),ellipsoid_MC_2D{j}(:,2),'-','Color',"#EDB120",'LineWidth',2);
            plot(mean_MC_2D{j}(1),mean_MC_2D{j}(2),'o','LineWidth',2,'MarkerSize',13,'Color','k','MarkerFaceColor','r');
        end
        for i = 1:length(Y_2D{j}(1,:))
            plot(Y_2D{j}(1,i),Y_2D{j}(2,i),'o','LineWidth',1.3,'Color',[0.3,0.3,1])
        end
        switch k
            case 1
            legend('3$\sigma$ Ellipse-Lincov','Mean-Lincov','MC Samples','Interpreter','latex','Location','best')
            case 2
            legend('3$\sigma$ Ellipse-UT','Mean-UT','MC Samples','Interpreter','latex','Location','best')
            case 3
            legend('3$\sigma$ Ellipse-MC','Mean-MC','MC Samples','Interpreter','latex','Location','best')
        end
        title(['Uncertainty propagation for ',num2str(sat_name{j})],'Interpreter','latex')
        if a == 1
            xlabel('$X_{perifocal}$[km]','Interpreter','latex')
            ylabel('$Y_{perifocal}$[km]','Interpreter','latex')
        else
            xlabel('$V_{x_{perifocal}}$[km/s]','Interpreter','latex')
            ylabel('$V_{y_{perifocal}}$[km/s]','Interpreter','latex')
        end
        axis equal
        grid on
        set(gca,'FontSize',18)
    end
end
Total_CPU_time = toc;
fprintf('Overall CPU time = %f [s]\n',Total_CPU_time)

% Clear kernels
cspice_kclear();
%% FUNCTIONS

function [xf,PHIf,tf,xx, tt]  = propagate(t_grid,x0,attractor)
% OBJECTIVE
% This function propagates the state from the initial condition x0 and the
% STM for a 3D 2BP

% INPUTS
% t_grid ---- [1xn] Time grid for propagation
% x0 ---- [6x1] Initial state
% attractor --- 'Name' Name of the attractor

% OUTPUTS 
% xf ---- [6x1] Final state
% PHIf ---- [6x6] Final STM
% tf --- [1x1] Final time
% xx ---- Result of the ODE at each step
% tt ---- All the time steps

    % Initialize propagation data
    if isfloat(attractor)
        GM = attractor;
    else
        GM = cspice_bodvrd(attractor,'GM',1);
    end

    % Initialize State Transition Matrix at t0
    Phi0 = eye(6);

    % Initial conditions the conditions for the STM
    x0Phi0 = [x0;Phi0(:)];
    
    % Perform integration
    options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12);
    [tt, xx] = ode78(@(t,x) keplerian_STM_rhs(t,x,GM),t_grid, x0Phi0, options_STM);

    % Extract state vector and State Transition Matrix
    xf = xx(end,1:6)';
    PHIf = reshape(xx(end,7:end),6,6);
    tf = tt(end);
    xx = xx';
end


function [dxdt] = keplerian_rhs(~, x, GM)
%KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian) propagator
%   Evaluates the right-hand-side of a newtonian 2-body propagator.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 24/10/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   t   : [ 1, 1] epoch (unused)
%   x   : [6, 1] cartesian state vector wrt Solar-System-Barycentre and
%                 State Transition Matrix elements
%   GM  : [ 1, 1] gravitational constant of the body
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
%

    % Initialize right-hand-side
    dxdt = zeros(6,1);
    
    % Extract positions
    rr = x(1:3);
    
    % Compute square distance and distance
    dist2 = dot(rr, rr);
    dist = sqrt(dist2);
    
    % Position detivative is object's velocity
    dxdt(1:3) = x(4:6);   
    % Compute the gravitational acceleration using Newton's law
    dxdt(4:6) = - GM * rr /(dist*dist2);
end


function [dxdt] = keplerian_STM_rhs(~,x,GM)
%KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian)
%               propagator with STM
%   Evaluates the right-hand-side of a newtonian 2-body propagator with STM.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 11/10/2023
%   Contact: alessandro.morselli@polimi.it
%   CoPy{j}right: (c) 2023 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2023/2024.
%
%
% Inputs:
%   t   : [ 1, 1] epoch (unused)
%   x   : [42, 1] cartesian state vector wrt Solar-System-Barycentre and
%                 State Transition Matrix elements
%   GM  : [ 1, 1] gravitational constant of the body
%
% Outputs:
%   dxdt   : [42,1] RHS, newtonian gravitational acceleration only
%

    % Initialize right-hand-side
    dxdt = zeros(42,1);
    
    % Extract positions
    rr = x(1:3);
    Phi = reshape(x(7:end),6,6);
    
    % Compute square distance and distance
    dist2 = dot(rr, rr);
    dist = sqrt(dist2);
    
    % Compute the gravitational acceleration using Newton's law
    aa_grav =  - GM * rr /(dist*dist2);
    
    % Compute the derivative of the flow
    dfdv = 3*GM/(dist^5)*(rr*rr') - GM/dist^3*eye(3);
    
    % Assemble the matrix A(t)=dfdx
    dfdx = [zeros(3), eye(3);
            dfdv, zeros(3)];
    % Compute the derivative of the state transition matrix
    Phidot = dfdx*Phi;
    
    dxdt(1:3) = x(4:6);   % Position detivative is object's velocity
    dxdt(4:6) = aa_grav;  % Sum up acceleration to right-hand-side
    dxdt(7:end) = Phidot(:);
end


function [y_mean,Py] = UT(par,x_mean,Px,t_grid)
% OBJECTIVE
% Propagation of mean value and covariance with unscented transform

% INPUTS
% par --- [1x3] [alpha,beta,k] parameters for UT computation
% x_mean --- [6x1] mean value of the previous step
% Px --- [6x6] Covariance of the previous step
% t_grid ---- [1xn] Time grid for propagation

% OUTPUTS
% y_mean --- [6x1] mean value of the current step
% Py --- [6x6] Covariance of the current step

    % Initialization
    alpha = par(1);
    beta = par(2);
    k = par(3);
    n = length(x_mean);
    c = alpha^2*(n+k);
    X = zeros(n,2*n+1);
    delta_X = zeros(n,2*n);
    Y = cell(1,length(t_grid));
    
    % SIGMA POINTS COMPUTATION
    delta_X(:,1:n) = sqrtm(c*Px);
    delta_X(:,n+1:2*n) = -sqrtm(c*Px);
    X(:,1) = x_mean;
    for i = 2 : (2*n+1)
        X(:,i) = x_mean + delta_X(:,i-1);
    end
    
    % WEIGHTS COMPUTATION
    W_m0 = 1 - n/c;
    W_c0 = (2-alpha^2+beta) - n/c;
    W_mi = 1/(2*c);
    W_ci = W_mi;
    
    % SIGMA POINTS PROPAGATION
    for i = 1 : (2*n+1)
        % Perform integration
        GM = cspice_bodvrd('Earth','GM',1);
        options = odeset('reltol', 1e-12, 'abstol', 1e-12);
        [~, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM),t_grid,X(:,i),options);
        YY = xx';
        for j = 1 : length(t_grid)
            Y{j}(:,i) = YY(1:6,j); 
        end
    end
    
    % WEIGHTED SAMPLE MEAN AND COVARIANCE COMPUTATION
    y_mean = zeros(6,length(t_grid));
    Py = cell(1,length(t_grid));
    y_mean(:,1) = x_mean;
    Py{1} = Px;
    for j = 2 : length(t_grid)
        y_mean(:,j) = W_m0*Y{j}(:,1) + W_mi*sum(Y{j}(:,2:end),2);
        Py{j} = W_c0*(Y{j}(:,1)-y_mean(:,j))*(Y{j}(:,1)-y_mean(:,j))';
        for i = 2 : (2*n+1)
            Py{j} = Py{j} + W_ci*(Y{j}(:,i)-y_mean(:,j))*(Y{j}(:,i)-y_mean(:,j))';
        end
    end
end


function [y_mean,Py,Y] = MC(N,x_mean,Px,t_grid)
% OBJECTIVE
% Propagation of mean value and covariance with Montecarlo simulation

% INPUTS
% N --- [1x1] Number of MC samples
% x_mean --- [6x1] mean value of the previous step
% Px --- [6x6] Covariance of the previous step
% t_grid ---- [1xn] Time grid for propagation

% OUTPUTS
% y_mean --- [6x1] mean value of the current step
% Py --- [6x6] Covariance of the current step
% Y --- [6xN] Matrix with propagated samples 


    % Inizialization
    Y = cell(1,length(t_grid));

    % GENERATE SAMPLES
    X = mvnrnd(x_mean,Px,N);
    X = X';
    
    % SAMPLES PROPAGATION
    for i = 1 : N
        % Perform integration
        GM = cspice_bodvrd('Earth','GM',1);
        options = odeset('reltol', 1e-12, 'abstol', 1e-12);
        [~, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM),t_grid,X(:,i),options);
        YY = xx';
        for j = 1 : length(t_grid)
            Y{j}(:,i) = YY(1:6,j); 
        end
    end
    
    % SAMPLE MEAN AND COVARIANCE COMPUTATION
    y_mean = zeros(6,length(t_grid));
    Py = cell(1,length(t_grid));
    y_mean(:,1) = x_mean;
    Py{1} = Px;
    for j = 2 : length(t_grid)
        y_mean(:,j) = sum(Y{j},2)/N;
        Py{j} = zeros(6);
        for i = 1 : N
            Py{j} = Py{j} + (Y{j}(:,i)-y_mean(:,j))*(Y{j}(:,i)-y_mean(:,j))';
        end
        Py{j} = Py{j}/(N-1);
    end
end


function R_eci2pf = rotation_eci2pf(r,v,GM)
% OBJECTIVE
% Computation of rotation matrix from ECI to Perifocal frame

% INPUTS 
% r --- [3x1] Orbital position vector
% v --- [3x1] Orbital velocity vector
% GM --- [1x1] Gravitational constant of the attractor

% OUTPUT
% R_eci2pf [3x3] Rotation matrix from ECI to Perifocal frame

    % KEPLERIAN ELEMENTS COMPUTATION
    K       = [0 0 1]';
    r_norm  = norm(r);                                 
    % Angular momentum
    h       = cross(r,v);                         
    h_norm  = norm(h);
    % Eccentricity
    e_vect  = (1/GM) * (cross(v,h) - GM*(r./r_norm));
    e       = norm(e_vect);
    % Inclination
    i       = acos(h(3)/h_norm);
    n       = cross(K,h)./norm(cross(K,h));
    % RAAN
    if n(2) > 0
        OM = acos(n(1));
    elseif n(2)<0
        OM = 2*pi - acos(n(1));
    end
    % Argument of pericentre
    if e_vect(3)>0
        om = acos(dot(n,e_vect)/e);
    elseif e_vect(3)<0
        om = 2*pi - acos(dot(n,e_vect)/e);
    end
    
    % ROTATION MATRIX FROM ECI TO PERIFOCAL FRAME
    R_OM = [cos(OM) sin(OM) 0;-sin(OM) cos(OM) 0;0 0 1];
    R_i = [1 0 0;0 cos(i) sin(i);0 -sin(i) cos(i)];
    R_om = [cos(om) sin(om) 0;-sin(om) cos(om) 0;0 0 1];
    R_eci2pf = R_om*R_i*R_OM;
end


function ellip = ellipsoid_computation(Pr,center)
% OBJECTIVE
% Computation the 3D ellispoid of a covariance matrix

% INPUT
% Pr --- [3x3] Covariance matrix
% center --- [3x1] Center of the ellipsoid (mean state)

% OUTPUT
% ellips [3xnum_points] 3D ellispoid of the covariance matrix

    % Find eigenvectors and eigenvalues of the covariance matrix
    [eigenvecs, eigenvals] = eig(Pr); 

    % Calculate scale factors for ellipsoid along principal axes
    scale_fact = 3*sqrt(diag(eigenvals)); 

    % Define number of points for ellipsoid
    num_points = 1000; 

    % Create sphere
    [x, y, z] = sphere(num_points); 

    % Scale and rotate the sphere to form ellipsoid
    ellip = eigenvecs * diag(scale_fact) * [x(:), y(:), z(:)]' + center;
end



