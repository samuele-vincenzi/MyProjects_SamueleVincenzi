% Spacecraft Guidance and Navigation (2023/2024)
% Assignment # 2, Exercise 2
% Author: Samuele Vincenzi

clear
clc
close all
format long g

% CODE INIZIALIZATION
addpath('sgp4')
cspice_furnsh('assignment02.tm');
typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

% PROBLEM INITIAL DATA
% Separation epoch
t_sep_utc = '2010-08-12T05:27:39.114';
t_sep_et =  cspice_str2et(t_sep_utc);
% Time window definition
t0_UTC = '2010-08-12T05:30:00.000';
t0_et = cspice_str2et(t0_UTC);
tf_UTC = '2010-08-12T11:00:00.000';
tf_et = cspice_str2et(tf_UTC);
t_g = t0_et:60:tf_et;

% Mean values of the state estimate at t0
% Mango
r0_mean_MANGO = [4622.232026629, 5399.3369588058, -0.0212138165769957]'; % [km]
v0_mean_MANGO = [0.812221125483763, -0.721512914578826, 7.42665302729053]'; % [km/s]
x0_mean_MANGO = [r0_mean_MANGO;v0_mean_MANGO];
% Tango
r0_mean_TANGO = [4621.69343340281, 5399.26386352847, -3.09039248714313]'; % [km]
v0_mean_TANGO = [0.813960847513811, -0.719449862738607, 7.42706066911294]'; % [km/s]
x0_mean_TANGO = [r0_mean_TANGO;v0_mean_TANGO];

% Define a struct for Mango data
MANGO(1).Satellite_name = 'MANGO';
MANGO(2).Satellite_name = 'MANGO';
MANGO(1).ID = '36599';
MANGO(2).ID = '36599';
MANGO(1).Station = 'KOUROU';
MANGO(2).Station = 'SVALBARD';

% Define a struct for Tango data
TANGO(1).Satellite_name = 'TANGO';
TANGO(2).Satellite_name = 'TANGO';
TANGO(1).ID = '36827';
TANGO(2).ID = '36827';
TANGO(1).Station = 'KOUROU';
TANGO(2).Station = 'SVALBARD';

% Define a struct for navigation solution
Nav_sol(1).TITLE = 'SOLUTION';
Nav_sol(2).TITLE = 'AVERAGE WEIGHTED RESIDUALS: Az[-] -- El[-] -- Range [-]';
Nav_sol(3).TITLE = 'ESTIMATE ERROR';
Nav_sol(4).TITLE = '3-SIGMA-POSITON [km], 3-SIGMA-VELOCITY [km/s]';
Nav_sol(1).SATELLITE = 'MANGO';
Nav_sol(2).SATELLITE = 'MANGO';
Nav_sol(3).SATELLITE = 'MANGO';
Nav_sol(4).SATELLITE = 'MANGO';

% Ground stations data
station(1).Name = 'KOUROU';
station(2).Name = 'SVALBARD';
station(1).LAT = 5.25144; % [deg]
station(2).LAT = 8.229772; % [deg]
station(1).LON = -52.80466; % [deg]
station(2).LON = 15.407786; % [deg]
station(1).ALT = -14.67; % [m]
station(2).ALT = 458; % [m]
station(1).noise_Az_El = 100*1e-3*cspice_rpd; % [rad]
station(2).noise_Az_El = 125*1e-3*cspice_rpd; % [rad]
station(1).noise_range = 0.01; % [km]
station(2).noise_range = 0.01; % [km]
station(1).min_Elev = 10; % [deg]
station(2).min_Elev = 5; % [deg]

% Earth gravitational constant
GM = cspice_bodvrd('Earth','GM',1); % [km^3/s^2]

%% POINT 1
% SOLUTION FOR MANGO
% MEAN STATE PROPAGATION WITH 2BP
[~,~,Mean_state,~]  = propagate(t_sep_et,x0_mean_MANGO,t_g,'Earth');
Mean_state = Mean_state(2:end,:)';

% COMPUTATION OF AZIMUTH, ELEVATION AND RANGE FROM t0 TO tf
Az = zeros(2,length(t_g));
El = zeros(2,length(t_g));
Range = zeros(2,length(t_g));
for i = 1 : length(t_g)
    [Az(1,i),El(1,i),Range(1,i)] = measurements_mean_value(station(1).Name,Mean_state(1:3,i),t_g(i));
    [Az(2,i),El(2,i),Range(2,i)] = measurements_mean_value(station(2).Name,Mean_state(1:3,i),t_g(i));
end

% PLOT OF ELEVATION TIME EVOLUTION WITH CONSTRAINT
figure
for i = 1 : 2
    subplot(1,2,i)
    hold on
    plot(t_g/cspice_spd(),El(i,:)*cspice_dpr(),'DisplayName',MANGO(1).Satellite_name,'LineWidth',1.5)
    yline(station(i).min_Elev,'Color','r','LineWidth',1.5);
    yl = ylim;
    xlim([t_g(1)/cspice_spd() t_g(end)/cspice_spd()])
    xl = xlim;
    xBox = [xl(1),xl(1),xl(2),xl(2),xl(1)];
    yBox = [station(i).min_Elev,yl(1),yl(1),station(i).min_Elev,station(i).min_Elev];
    patch(xBox,yBox,'black','FaceColor','k','FaceAlpha',0.2,'EdgeColor','none')
    title(['@',station(i).Name])
    xlabel('Epoch [MJD2000]','Interpreter','latex')
    ylabel('Elevation [deg]','Interpreter','latex')
    legend(MANGO(1).Satellite_name,'Minimum elevation','Non-visible','Interpreter','latex','Location','best')
    set(gca,'FontSize',20)
    grid on
end

% VISIBILITY WINDOWS COMPUTATION
i_visibility = El > [station(1).min_Elev;station(2).min_Elev]*cspice_rpd;
pos = localize_position(i_visibility); % Retrieve the indexes correspondent to windows
for i = 1 : 2 % For both ground stations
    for j = 1 : length(pos{i})
        MANGO(i).Vis_windows{j} = t_g(pos{i}{j}); % Division and selection of windows
    end
end

% VISIBILITY WINDOWS INITIAL AND FINAL DATES (UTC)
for i = 1 : 2 % For both ground stations
    fprintf('\nVISIBILITY WINDOWS FROM: %s\n',station(i).Name)
    for j = 1 : length(MANGO(i).Vis_windows)
        fprintf('%d) Initial date: %s UTC\n',j,cspice_et2utc(MANGO(i).Vis_windows{j}(1),'C',4))
        fprintf('   Final date: %s UTC\n',cspice_et2utc(MANGO(i).Vis_windows{j}(end),'C',4))
    end
end


% PLOT OF AZIMUTH TIME EVOLUTION WITH VISIBILITY WINDOWS
figure
for i = 1 : 2
    subplot(1,2,i)
    hold on
    plot(t_g/cspice_spd(),Az(i,:)*cspice_dpr(),'DisplayName',MANGO(1).Satellite_name,'LineWidth',1.5)
    yl = ylim;
    xlim([t_g(1)/cspice_spd() t_g(end)/cspice_spd()])
    xl = xlim;
    yBox = [yl(1),yl(1),yl(2),yl(2),yl(1)];
    xBox1 = [MANGO(i).Vis_windows{1}(1)/cspice_spd(),xl(1),xl(1),MANGO(i).Vis_windows{1}(1)/cspice_spd(),MANGO(i).Vis_windows{1}(1)/cspice_spd()];
    xBox2 = [MANGO(i).Vis_windows{1}(end)/cspice_spd(),MANGO(i).Vis_windows{2}(1)/cspice_spd(),MANGO(i).Vis_windows{2}(1)/cspice_spd(),MANGO(i).Vis_windows{1}(end)/cspice_spd(),MANGO(i).Vis_windows{1}(end)/cspice_spd()];
    patch(xBox1,yBox,'black','FaceColor','k','FaceAlpha',0.2,'EdgeColor','none')
    patch(xBox2,yBox,'black','FaceColor','k','FaceAlpha',0.2,'EdgeColor','none')
    if i == 1
        xBox3 = [MANGO(i).Vis_windows{2}(end)/cspice_spd(),xl(2),xl(2),MANGO(i).Vis_windows{2}(end)/cspice_spd(),MANGO(i).Vis_windows{2}(end)/cspice_spd()];
        patch(xBox3,yBox,'black','FaceColor','k','FaceAlpha',0.2,'EdgeColor','none')
    end
    if i == 2
        xBox3 = [MANGO(i).Vis_windows{2}(end)/cspice_spd(),MANGO(i).Vis_windows{3}(1)/cspice_spd(),MANGO(i).Vis_windows{3}(1)/cspice_spd(),MANGO(i).Vis_windows{2}(end)/cspice_spd(),MANGO(i).Vis_windows{2}(end)/cspice_spd()];
        xBox4 = [MANGO(i).Vis_windows{3}(end)/cspice_spd(),MANGO(i).Vis_windows{4}(1)/cspice_spd(),MANGO(i).Vis_windows{4}(1)/cspice_spd(),MANGO(i).Vis_windows{3}(end)/cspice_spd(),MANGO(i).Vis_windows{3}(end)/cspice_spd()];
        xBox5 = [MANGO(i).Vis_windows{4}(end)/cspice_spd(),xl(2),xl(2),MANGO(i).Vis_windows{4}(end)/cspice_spd(),MANGO(i).Vis_windows{4}(end)/cspice_spd()];
        patch(xBox3,yBox,'black','FaceColor','k','FaceAlpha',0.2,'EdgeColor','none')
        patch(xBox4,yBox,'black','FaceColor','k','FaceAlpha',0.2,'EdgeColor','none')
        patch(xBox5,yBox,'black','FaceColor','k','FaceAlpha',0.2,'EdgeColor','none')
    end
    title(['@',station(i).Name])
    xlabel('Epoch [MJD2000]','Interpreter','latex')
    ylabel('Azimuth [deg]','Interpreter','latex')
    set(gca,'FontSize',20)
    legend(MANGO(1).Satellite_name,'Outside visibility window','Interpreter','latex','Location','best','FontSize',15)
    grid on
end

% PLOT OF AZIMUTH AND ELEVATION PROFILE FOR EACH STATION AND VISIBILITY WINDOW
figure
for i = 1 : 2
    subplot(1,2,i)
    hold on
    len = length(MANGO(i).Vis_windows);
    for j = 1 : len
        in = pos{i}{j}(1);
        fin = pos{i}{j}(end);
        plot(Az(i,in:fin)*cspice_dpr(),El(i,in:fin)*cspice_dpr(),'*','MarkerSize',10,'LineWidth',1)
    end
    if i == 1
    legend('First visibility window','Second visibility window','Interpreter','latex')
    else
        legend('First visibility window','Second visibility window','Third visibility window','Fourth visibility window','Interpreter','latex')
    end
    axis([-180,180,0,70]);
    title(station(i).Name,'Interpreter','latex')
    xlabel('Azimuth [deg]','Interpreter','latex')
    ylabel('Elevation [deg]','Interpreter','latex')
    set(gca,'FontSize',19)
    grid on
end

%% POINT 2
% DATA COMPUTATION FROM TLE
% Load TLE file and create "satrec" structure
TLE_file = 'tle\36599.3le';
satrec = read_3LE(MANGO(1).ID,TLE_file, whichconst);
% Compute the epoch
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

% MEASUREMENTS SIMULATION
pos_new = cell(1,2);
Meas_Az_plot = cell(1,2); % Ausiliary Azimuth variable
Meas_El_plot = cell(1,2); % Ausiliary Elevation variable
R = cell(1,2);
for i = 1 : 2 % Measurements computation for each station
    len = length(MANGO(i).Vis_windows);
    % Measurements Covariance matrix
    R{i} = diag([station(i).noise_Az_El^2 ; ...
                        station(i).noise_Az_El^2; ...
                        station(i).noise_range^2]);
    % Measurements initialization
    MANGO(i).Meas_Az = [];
    MANGO(i).Meas_El = [];
    MANGO(i).Meas_Range = [];
    pos_new{i} = [];
    for j = 1 : len % Measurements computation for each visibility window
        [Meas_Az,Meas_El,Meas_Range,i_visibility] = TLE_2_measurements(satrec,sat_epoch_et,station(i),MANGO(i).Vis_windows{j},R{i});
        pos_new{i} = [pos_new{i} pos{i}{j}(i_visibility)];
        Meas_Az_plot{i}{j} = Meas_Az*cspice_dpr; % Ausiliary Azimuth variable
        Meas_El_plot{i}{j} = Meas_El*cspice_dpr; % Ausiliary Elevation variable
        MANGO(i).Meas_Az = [MANGO(i).Meas_Az Meas_Az*cspice_dpr]; % Azimuth
        MANGO(i).Meas_El = [MANGO(i).Meas_El Meas_El*cspice_dpr]; % Elevation
        MANGO(i).Meas_Range = [MANGO(i).Meas_Range Meas_Range]; % Range
    end
    MANGO(i).Meas_Time_et = t_g(pos_new{i}); % Measurements time (ET)
end

% PLOT OF SIMULATED MEASUREMENTS FOR EACH STATION AND VISIBILITY WINDOW
pos_Name = {'First','Second','Third','Fourth'};
for i = 1 : 2
    figure
    len = length(MANGO(i).Vis_windows);
    for j = 1 : len
        subplot(len/2,2,j)
         hold on
        in = pos{i}{j}(1);
        fin = pos{i}{j}(end);
        plot(Az(i,in:fin)*cspice_dpr(),El(i,in:fin)*cspice_dpr(),'.','MarkerSize',35,'Color',"#00FFFF")
        plot(Meas_Az_plot{i}{j},Meas_El_plot{i}{j},'*','MarkerSize',11,'LineWidth',1.2)
        yline(station(i).min_Elev,'Color',"#77AC30",'LineWidth',1);
        axis([-180,180,0,40])
        xl = xlim;
        yl = ylim;
        xBox = [xl(1),xl(1),xl(2),xl(2),xl(1)];
        yBox = [station(i).min_Elev,yl(1),yl(1),station(i).min_Elev,station(i).min_Elev];
        patch(xBox,yBox,'black','FaceColor','k','FaceAlpha',0.1,'EdgeColor','none')
        title(['@',station(i).Name,' - ',pos_Name{j},' visibility window'])
        xlabel('Azimuth [deg]','Interpreter','latex')
        ylabel('Elevation [deg]','Interpreter','latex')
        legend('Predicted Az-El profile','Simulated measurements','',['$El<$',num2str(station(i).min_Elev),'$^{\circ}$'],'Interpreter','latex')
        set(gca,'FontSize',19)
        grid on
    end
end

%% POINT 3
% BATCH FILTER : NAVIGATION SOLUTION
% Computation of the state at t0 from TLE propagation
[r0_star,v0_star] = TLE_2_ECI(satrec,sat_epoch_et,t0_et);
x_star = [r0_star;v0_star];

% COMPUTATION OF COST FUNCTIONs FOR EACH CASE
fun_kep = cell(1,2);
fun_j2 = cell(1,2);
meas_real = cell(1,2);
tspan = cell(1,2);
for i = 1 : 2 % For both stations
    W_m = inv(sqrtm(R{i}));
    tspan{i} = [t0_et MANGO(i).Meas_Time_et];
    % Define a matrix with estimated real measurements
    meas_real{i} = [MANGO(i).Meas_Az*cspice_rpd; 
                 MANGO(i).Meas_El*cspice_rpd ;
                 MANGO(i).Meas_Range]';
    % Cost function computation
    fun_kep{i} = @(x) cost_function_lsq(x,tspan{i},W_m,meas_real{i},station(i).Name,0); % Without J2
    fun_j2{i} = @(x) cost_function_lsq(x,tspan{i},W_m,meas_real{i},station(i).Name,1); % With J2
end

% SUB-POINT a) NAVIGATION SOLUTION FOR KOUROU AND KEP. DYNAMICS
fun_a = @(x) fun_kep{1}(x);
% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
[x_lsq_a,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun_a,x_star,[],[],options); % Without J2
Nav_sol(1).KOUROU = x_lsq_a'; % Navigation solution for Kourou
I = eye(3);
Nav_sol(2).KOUROU = mean(abs(residual),1); % Average value of the residuals
Nav_sol(3).KOUROU = abs((x_lsq_a - x_star)'); % Navigation solution error
Jac = full(jacobian);
P_ls = resnorm/(length(residual)-length(x_star)).*inv(Jac.'*Jac);
Nav_sol(4).KOUROU = [3*sqrt(max(eig(P_ls(1:3,1:3)))),3*sqrt(max(eig(P_ls(4:6,4:6))))]; % 3sigma


% SUB-POINT b) NAVIGATION SOLUTION FOR KOUROU+SVALBARD AND KEP. DYNAMICS
fun_b = @(x) [fun_kep{1}(x); fun_kep{2}(x)];
[x_lsq_b,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun_b,x_star,[],[],options); % Without J2
Nav_sol(1).KOUROU_SVALBARD = x_lsq_b'; % Navigation solution for Kourou+Svalbard
Nav_sol(2).KOUROU_SVALBARD = mean(abs(residual),1); % Average value of the residuals
Nav_sol(3).KOUROU_SVALBARD = abs((x_lsq_b - x_star)'); % Navigation solution error
Jac = full(jacobian);
P_ls = resnorm/(length(residual)-length(x_star)).*inv(Jac.'*Jac);
Nav_sol(4).KOUROU_SVALBARD = [3*sqrt(max(eig(P_ls(1:3,1:3)))),3*sqrt(max(eig(P_ls(4:6,4:6))))]; % 3sigma

% SUB-POINT c) NAVIGATION SOLUTION FOR KOUROU+SVALBARD AND KEP. DYNAMICS+J2
fun_c = @(x) [fun_j2{1}(x); fun_j2{2}(x)];
[x_lsq_c,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun_c,x_star,[],[],options); % With J2
Nav_sol(1).KOUROU_SVALBARD_J2 = x_lsq_c'; % Navigation solution for Kourou+Svalbard+J2
Nav_sol(2).KOUROU_SVALBARD_J2 = mean(abs(residual),1); % Average value of the residuals
Nav_sol(3).KOUROU_SVALBARD_J2 = abs((x_lsq_c - x_star)'); % Navigation solution error
Jac = full(jacobian);
P_ls = resnorm/(length(residual)-length(x_star)).*inv(Jac.'*Jac);
Nav_sol(4).KOUROU_SVALBARD_J2 = [3*sqrt(max(eig(P_ls(1:3,1:3)))),3*sqrt(max(eig(P_ls(4:6,4:6))))]; % 3sigma

%% POINT 5
% SOLUTION FOR TANGO
% MEAN STATE PROPAGATION WITH 2BP
[~,~,Mean_state,~]  = propagate(t_sep_et,x0_mean_TANGO,t_g,'Earth');
Mean_state = Mean_state(2:end,:)';

% COMPUTATION OF AZIMUTH, ELEVATION AND RANGE FROM t0 TO tf
Az = zeros(2,length(t_g));
El = zeros(2,length(t_g));
Range = zeros(2,length(t_g));
for i = 1 : length(t_g)
    [Az(1,i),El(1,i),Range(1,i)] = measurements_mean_value(station(1).Name,Mean_state(1:3,i),t_g(i));
    [Az(2,i),El(2,i),Range(2,i)] = measurements_mean_value(station(2).Name,Mean_state(1:3,i),t_g(i));
end

% VISIBILITY WINDOWS COMPUTATION
i_visibility = El > [station(1).min_Elev;station(2).min_Elev]*cspice_rpd;
pos = localize_position(i_visibility);
for i = 1 : 2 % For both ground stations
    for j = 1 : length(pos{i})
        TANGO(i).Vis_windows{j} = t_g(pos{i}{j});
    end
end

% DATA COMPUTATION FROM TLE
% Load TLE file and create "satrec" structure
TLE_file = 'tle\36827.3le';
satrec = read_3LE(TANGO(1).ID,TLE_file, whichconst);
% Compute the epoch
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

% MEASUREMENTS SIMULATION
pos_new = cell(1,2);
Meas_Az_plot = cell(1,2); % Ausiliary Azimuth variable
Meas_El_plot = cell(1,2); % Ausiliary Elevation variable
R = cell(1,2);
for i = 1 : 2 % Measurements computation for each station
    m = 1;
    len = length(TANGO(i).Vis_windows);
    % Measurements Covariance matrix
    R{i} = diag([station(i).noise_Az_El^2 ; ...
                 station(i).noise_Az_El^2; ...
                 station(i).noise_range^2]);
    % Measurements initialization
    TANGO(i).Meas_Az = [];
    TANGO(i).Meas_El = [];
    TANGO(i).Meas_Range = [];
    pos_new{i} = [];
    for j = 1 : len % Measurements computation for each visibility window
        [Meas_Az,Meas_El,Meas_Range,i_visibility] = TLE_2_measurements(satrec,sat_epoch_et,station(i),TANGO(i).Vis_windows{j},R{i});
        pos_new{i} = [pos_new{i} pos{i}{j}(i_visibility)];
        Meas_Az_plot{i}{j} = Meas_Az*cspice_dpr; % Ausiliary Azimuth variable
        Meas_El_plot{i}{j} = Meas_El*cspice_dpr; % Ausiliary Elevation variable
        TANGO(i).Meas_Az = [TANGO(i).Meas_Az Meas_Az*cspice_dpr]; % Azimuth
        TANGO(i).Meas_El = [TANGO(i).Meas_El Meas_El*cspice_dpr]; % Elevation
        TANGO(i).Meas_Range = [TANGO(i).Meas_Range Meas_Range]; % Range
    end
    TANGO(i).Meas_Time_et = t_g(pos_new{i}); % Measurements time (ET)
end

% BATCH FILTER : NAVIGATION SOLUTION
% Computation of the state at t0 from TLE propagation
[r0_star,v0_star] = TLE_2_ECI(satrec,sat_epoch_et,t0_et);
x_star = [r0_star;v0_star];

% Computation of the solution
fun_kep = cell(1,2);
fun_j2 = cell(1,2);
meas_real = cell(1,2);
tspan = cell(1,2);
for i = 1 : 2 % For both stations
    W_m = inv(sqrtm(R{i}));
    tspan{i} = [t0_et TANGO(i).Meas_Time_et];
    % Define a matrix with estimated real measurements
    meas_real{i} = [TANGO(i).Meas_Az*cspice_rpd; 
                 TANGO(i).Meas_El*cspice_rpd ;
                 TANGO(i).Meas_Range]';
    % Cost function computation
    fun_j2{i} = @(x) cost_function_lsq(x,tspan{i},W_m,meas_real{i},station(i).Name,1); % With J2
end

% NAVIGATION SOLUTION FOR KOUROU+SVALBARD AND KEP. DYNAMICS+J2
Nav_sol(1).SATELLITE_ = 'TANGO';
Nav_sol(2).SATELLITE_ = 'TANGO';
Nav_sol(3).SATELLITE_ = 'TANGO';
Nav_sol(4).SATELLITE_ = 'TANGO';
fun_c = @(x) [fun_j2{1}(x); fun_j2{2}(x)];
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
[x_lsq_c,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun_c,x_star,[],[],options); % With J2
Nav_sol(1).KOUROU_SVALBARD_J2_ = x_lsq_c'; % Navigation solution for Kourou+Svalbard+J2
Nav_sol(2).KOUROU_SVALBARD_J2_ = mean(abs(residual),1); % Average value of the residuals
Nav_sol(3).KOUROU_SVALBARD_J2_ = abs((x_lsq_c - x_star)'); % Navigation solution error
Jac = full(jacobian);
P_ls = resnorm/(length(residual)-length(x_star)).*inv(Jac.'*Jac);
Nav_sol(4).KOUROU_SVALBARD_J2_ = [3*sqrt(max(eig(P_ls(1:3,1:3)))),3*sqrt(max(eig(P_ls(4:6,4:6))))]; % 3sigma

% CLEAR KERNELS
cspice_kclear()
%% FUNCTIONS
function [xf,tf,xx, tt]  = propagate(t0,x0,tf,attractor)
% DESCRIPTION
% This function propagates the state from the initial condition x0 for a 3D 2BP

% INPUTS
% t0 ---- [1x1] Initial propagation time
% x0 ---- [6x1] Initial state
% tf ---- [1x1] Final propagation time
% attractor ---- 'Name', Name of the attractor

% OUTPUTS 
% xf ---- [6x1] Final state
% tf --- [1x1] Final time
% xx ---- Result of the ODE at each step
% tt ---- All the time steps

    % Initialize propagation data
    GM = cspice_bodvrd(attractor,'GM',1);
    
    % Perform integration
    options = odeset('reltol', 1e-12, 'abstol', 1e-12);
    [tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM), [t0 tf], x0, options);

    % Extract state vector and State Transition Matrix
    xf = xx(end,1:6)';
    tf = tt(end);
end


function [dxdt] = keplerian_rhs(t, x, GM,type)
% DESCRIPTION
% If type=0 the function evaluates the right-hand-side of a 2-body 
% (keplerian) propagator.
% If type=1 the function evaluates the right-hand-side of a 2-body 
% (keplerian) propagator with J2 perturbation.

% INPUTS:
% t --- [1x1] Epoch (used only for type=1)
% x --- [6x1] Cartesian state vector wrt attractor Barycentre
% GM --- [1x1] Gravitational constant of the attractor
% type --- [1x1] type=0 -> keplerian propagator 
%                type=1 -> keplerian propagator + J2

% OUTPUTS:
% dxdt --- [6x1] RHS

    if nargin == 3
            type = 0;
    end
    
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
    
    switch type
        case 0
            a_J2_eci = [0;0;0];
        case 1
            % Compute J2 contribution to acceleration
            Re_vec = cspice_bodvrd('Earth','RADII',3);
            Re = Re_vec(1); % Earth equatorial radius
            J2 = 0.0010826269;
            rotm = cspice_pxform('J2000','ITRF93',t);
            pos_ecef = rotm*rr;
            r = norm(pos_ecef);
            a_J2_ecef = 3/2*GM*J2 * pos_ecef/r^3 * (Re/r)^2 .* (5*(pos_ecef(3)/r)^2 - [1;1;3]);
            a_J2_eci = rotm'*a_J2_ecef;
    end
    
    % Update the acceleration
    dxdt(4:6) = dxdt(4:6) + a_J2_eci;
end


function [meas_Az,meas_El,meas_Range,i_visibility] = TLE_2_measurements(satrec,TLE_epoch_et,station,Vis_window,R)
% DESCRIPTION
% Given the TLE evaluation, the function computes the estimated real
% measurements

% INPUTS
% satrec --- Struct obtained from TLE extraction
% TLE_epoch_et --- [1x1] Date (ET) of the TLE data
% station --- Struct with ground station name and minimum visible elevation
% Vis_window --- [1xnpoints] Visibility window
% R --- [3x3] Measurements covariance matrix

% OUTPUTS 
% meas_Az --- [1xnpoints] Azimuth for each point of the visibility window
% meas_El --- [1xnpoints] Elevation for each point of the visibility window
% meas_Range --- [1xnpoints] Range for each point of the visibility window
% i_visibility --- [1xnpoints] Analysis of visibility constraint

    % Orbit propagation with SGP4
    tspan = Vis_window;
    npoints = length(tspan);
    r_eci = zeros(3,npoints);
    v_eci = zeros(3,npoints);
    Az_mean = zeros(1,npoints);
    El_mean = zeros(1,npoints);
    Range_mean = zeros(1,npoints);
    Az = zeros(1,npoints);
    El = zeros(1,npoints);
    Range = zeros(1,npoints);
    for i = 1:npoints
        % ECI position and velocity
        [r_eci(:,i),v_eci(:,i)] = TLE_2_ECI(satrec,TLE_epoch_et,tspan(i));
        
        % Measurements simulation
        [Az_mean(i),El_mean(i),Range_mean(i)] = measurements_mean_value(station.Name,r_eci(:,i),tspan(i));

        % Simulation of measurement uncertainty
        y_hat = [Az_mean(i);El_mean(i);Range_mean(i)];
        Measurements = mvnrnd(y_hat,R);
        Az(i) = Measurements(1);
        El(i) = Measurements(2);
        Range(i) = Measurements(3);
    end
    % Visibility window check and measurements selection
    i_visibility = El*cspice_dpr>station.min_Elev;
    meas_Az = Az(i_visibility);
    meas_El = El(i_visibility);
    meas_Range = Range(i_visibility);
end


function [Az,El,Range] = measurements_mean_value(station_name,rr_sat_ECI,time_et)
% DESCRIPTION
% Given the position of the satellite the function computes Azimuth, 
% Elevation and Range.

% INPUTS
% station_name --- 'Name', Name of thr ground station
% rr_sat_ECI --- [3x1] Satellite position in ECI coordinates
% time_et --- [1x1] Date defined in ET

% OUTPUTS
% Az --- [1x1] Azimuth
% El --- [1x1] Elevation
% Range --- [1x1] Range

    % TOPOCENTRIC CARTESIAN COORDINATES OF THE SATELLITE
    % Define station name
    topoFrame = [station_name,'_TOPO'];
    % Transformation from ECI to topocentric frame
    ROT_ECI2TOPO = cspice_pxform('J2000',topoFrame,time_et);
    rr_sat_TOPO = ROT_ECI2TOPO*rr_sat_ECI;
    
    % TOPOCENTRIC CARTESIAN COORDINATES OF THE STATION
    rr_stat_ECI = cspice_spkpos(station_name,time_et,'J2000','NONE','EARTH');
    rr_stat_TOPO = ROT_ECI2TOPO*rr_stat_ECI;
    
    % COMPUTATION OF AZIMUTH,ELEVATION,RANGE
    rr_stat_sat_TOPO = rr_sat_TOPO - rr_stat_TOPO;
    [Range,Az,El] = cspice_reclat(rr_stat_sat_TOPO);
end


function [r_eci,v_eci] = TLE_2_ECI(satrec,TLE_epoch_et,time_et)
% DESCRIPTION
% Given the TLE evaluation, the function computes the ECI state of the
% target sateellite.

% INPUTS
% satrec --- Struct obtained from TLE extraction
% TLE_epoch_et --- [1x1] Date (ET) of the TLE data
% time_et --- [1x1] Date (ET) of the state computation

% OUTPUTS
% r_eci --- [3x1] Satellite position in a ECI reference frame
% v_eci --- [3x1] Satellite velocity in a ECI reference frame

    % Set nutation corrections parameters (from Celestrac EOP data)
    arcsec2rad = pi / (180*3600); % Constant for arcseconds to radians conversions
    ddpsi = -0.073296*arcsec2rad; %  [rad] (from Celestrac EOP data)
    ddeps = -0.009373*arcsec2rad; %  [rad] (from Celestrac EOP data)
    
    % SGP4 propagation
    tsince = (time_et - TLE_epoch_et)/60.0; % minutes from TLE epoch
    [~,r_teme,v_teme] = sgp4(satrec,tsince);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt = cspice_unitim(time_et,'ET','TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [r_eci,v_eci] = teme2eci(r_teme,v_teme,[0.0;0.0;0.0],ttt,ddpsi,ddeps);
end


function residual = cost_function_lsq(x,tspan,W_m,meas_real,station_name,type)
% DESCRIPTION
% Cost function for a weighted least-square solution of the navigation 
% problem.

% INPUTS
% x --- [3x1] Unknown of the navigation problem
% W_m --- [3x3] Matrix of the weights
% meas_real --- [nx3] Estimated real measurements
% station_name --- 'Name', Name of the station
% type --- [1x1] type = 0 for 2BP dynamical model
%                type = 1 for 2BP+J2 dynamical model

% OUTPUTS
% residual --- [nx3] Difference between predicted and real measurements

    residual = zeros(size(meas_real)); % Initialize output variable
    
    % Propagate x to the epochs of the measurements
    GM = cspice_bodvrd('Earth','GM',1);
    fun = @(t,x) keplerian_rhs(t,x,GM,type);
    options = odeset('Reltol',1.e-13,'Abstol',1.e-20);
    [~,x_prop] = ode78(fun,tspan,x,options);
    tspan = tspan(2:end);
    r_prop = x_prop(2:end,1:3);
    len = length(r_prop(:,1));
    
    % Compute predicted measurements 
    meas_pred = zeros(len-1,3);
    for i = 1 : len
        [Az,El,Range] = measurements_mean_value(station_name,r_prop(i,:)',tspan(i));
        meas_pred(i,:) = [Az El Range];
    end
    
    % Compute the residual of the measurements and append it to the output
    for k = 1 : length(tspan)
        diff_meas = [angdiff(meas_pred(k,1),meas_real(k,1)); 
                    angdiff(meas_pred(k,2),meas_real(k,2));
                    meas_pred(k,3) - meas_real(k,3)];
        diff_meas_weighted = W_m*diff_meas;
        residual(k,:) = diff_meas_weighted';
    end
end


function pos = localize_position(i_visibility)
% DESCRIPTION
% Ausiliary function for index detection of visibility window

% INPUT 
% i_visibility --- [1xn] Vector with 1 in positions that respect the
% constraint, 0 elsewhere

% OUTPUT
% pos --- cell(1,2) Variable that contains the visibility windows indexes
% wrt the initial time grid
    len = length(i_visibility(:,1));
    pos = cell(1,len);
    for i = 1 : len
        k = 0;
        j = 1;
        while j<=length(i_visibility(i,:))
            m = 1;
            if i_visibility(i,j) == 1
                k = k + 1;
            end
            while j<=length(i_visibility(i,:)) && i_visibility(i,j) == 1
                pos{i}{k}(m) = j;
                m = m + 1;
                j = j + 1;
            end
            while j<=length(i_visibility(i,:)) && i_visibility(i,j) == 0 
                j = j + 1;
            end
        end
    end
end

