% Spacecraft Guidance and Navigation (2023/2024)
% Assignment # 2, Exercise 3
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
tf_UTC = '2010-08-12T06:30:00.000';
tf_et = cspice_str2et(tf_UTC);
t_step = 5; % [s]
t_g = t0_et:t_step:tf_et;

% Mean values of the state estimate at t_sep
% Mango
r0_mean_MANGO = [4622.232026629, 5399.3369588058, -0.0212138165769957]'; % [km]
v0_mean_MANGO = [0.812221125483763, -0.721512914578826, 7.42665302729053]'; % [km/s]
x0_mean_MANGO = [r0_mean_MANGO;v0_mean_MANGO];
% Tango
r0_mean_TANGO = [4621.69343340281, 5399.26386352847, -3.09039248714313]'; % [km]
v0_mean_TANGO = [0.813960847513811, -0.719449862738607, 7.42706066911294]'; % [km/s]
x0_mean_TANGO = [r0_mean_TANGO;v0_mean_TANGO];

% Covariance matrix of the state estimate at t_sep
P0 = [+5.6e-7 +3.5e-7 -7.1e-8 0        0        0;
      +3.5e-7 +9.7e-7 +7.6e-8 0        0        0;
      -7.1e-8 +7.6e-8 +8.1e-8 0        0        0;
      0       0       0       +2.8e-11 0        0;
      0       0       0       0        +2.7e-11 0
      0       0       0       0        0        +9.6e-12]; % [km2, km2/s, km2/s2]

% Ground stations data -> Svalbard
station.Name = 'SVALBARD';
station.LAT = 8.229772; % [deg]
station.LON = 15.407786; % [deg]
station.ALT = 458; % [m]
station.noise_Az_El = 125*1e-3*cspice_rpd; % [rad]
station.noise_range = 0.01; % [km]
station.noise_Az_El_new = 1*cspice_rpd; % [rad]
station.noise_range_new = 1e-5; % [km]
station.min_Elev = 5; % [deg]

% Define a struct for Mango data
MANGO.Satellite_name = 'MANGO';
MANGO.ID = '36599';

% Define a struct for Tango data
TANGO.Satellite_name = 'TANGO';
TANGO.ID = '36827';

% Earth gravitational constant
GM = cspice_bodvrd('Earth','GM',1); % [km^3/s^2]

%% POINT 1: ESTIMATE MANGO ABSOLUTE STATE
tic
% ------- SUB-POINT a) PREDICTED VISIBILITY WINDOW ------- %
% MEAN STATE PROPAGATION WITH 2BP
[~,~,Mean_state,~]  = propagate(t_sep_et,x0_mean_MANGO,t_g,'Earth');
Mean_state = Mean_state(2:end,:)';

% COMPUTATION OF AZIMUTH, ELEVATION AND RANGE FROM t0 TO tf
Az = zeros(1,length(t_g));
El = zeros(1,length(t_g));
Range = zeros(1,length(t_g));
for i = 1 : length(t_g)
    [Az(i),El(i),Range(i)] = measurements_mean_value(station.Name,Mean_state(1:3,i),t_g(i));
end

% VISIBILITY WINDOW COMPUTATION
i_visibility = El > station.min_Elev*cspice_rpd;
pos = localize_position(i_visibility); % Retrieve the indexes correspondent to the window
vis_window_pred = t_g(pos{1}{1}); % Predicted visibility window
Pred_Az = Az(pos{1}{1}); % Predicted Azimuth
Pred_El = El(pos{1}{1}); % Predicted Elevation
% Print visibility window inital and final date (UTC)
fprintf('\nVISIBILITY WINDOW FROM: %s\n',station.Name)
fprintf('Initial date: %s UTC\n',cspice_et2utc(vis_window_pred(1),'C',4))
fprintf('Final date: %s UTC\n',cspice_et2utc(vis_window_pred(end),'C',4))

% ------- SUB-POINT b) ACTUAL MEASUREMENTS SIMULATION ------- %
% DATA COMPUTATION FROM TLE
% Load TLE file and create "satrec" structure
TLE_file = 'tle\36599.3le';
satrec = read_3LE(MANGO.ID,TLE_file, whichconst);
% Compute the epoch
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

% MEASUREMENTS SIMULATION
% Measurements Covariance matrix
R = diag([station.noise_Az_El^2 ; ...
          station.noise_Az_El^2; ...
          station.noise_range^2]);
% Measurements computation (actual measurements)
[Meas_Az,Meas_El,Meas_Range,i_visibility] = TLE_2_measurements(satrec,sat_epoch_et,station,vis_window_pred,R);
pos_new = pos{1}{1}(i_visibility); % Discard the indexes that does not respect the constraint
vis_window_act = t_g(pos_new); % Actual visibility window

% PLOT OF PREDICTED AND ACTUAL MEASUREMENTS
figure
hold on
plot(Pred_Az*cspice_dpr(),Pred_El*cspice_dpr(),'o','MarkerSize',9,'LineWidth',1)
plot(Meas_Az*cspice_dpr(),Meas_El*cspice_dpr(),'*','MarkerSize',9,'LineWidth',1)
axis([10,145,2,27])
title('Satellite: Mango, Ground station: Svalabard','Interpreter','latex')
legend('Predicted measurements','Actual measurements','Interpreter','latex','Location','best')
xlabel('Azimuth [deg]','Interpreter','latex')
ylabel('Elevation [deg]','Interpreter','latex')
set(gca,'FontSize',22)
grid on

% ------- SUB-POINT c) MANGO STATE ESTIMATE WITH UKF ------- %
% Parameters for UT and UKF
alpha = 0.1;
beta = 2;
k = 0;
par = [alpha,beta,k];

% Variables inizialization
len = length(vis_window_act);
x_mean_abs_Mango = zeros(6,len);
P_abs_Mango = cell(1,len);
sigma3_pos = zeros(1,len);
sigma3_vel = zeros(1,len);
err_pos = zeros(1,len);
err_vel = zeros(1,len);

% UNCERTAINTY PROPAGATION WITH UT (from t_sep to 5 seconds before the first
% measurement)
t_in = t_sep_et;
t_fin = vis_window_act(1) - t_step;
[x_mean_plus,P_plus] = UT(par,x0_mean_MANGO,P0*1e4,t_in,t_fin);

% STATE ESTIMATE WITH UKF (starting from 5 seconds before the first measurement)
t_in = t_fin; % Time at step 0
for i = 1 : len
    t_fin = vis_window_act(i); % Time at step k 
    y = [Meas_Az(i) Meas_El(i) Meas_Range(i)]'; % Actual measuremnts
    [x_mean_abs_Mango(:,i),P_abs_Mango{i}] = UKF(x_mean_plus,P_plus,y,R,t_in,t_fin,station.Name,par,0);

    % Update intial conditions for the following step
    x_mean_plus = x_mean_abs_Mango(:,i);
    P_plus = P_abs_Mango{i};
    t_in = t_fin; % Time at step k-1

    % 3sigma of the estimated covariavce
    sigma3_pos(i) = 3*sqrt(trace(P_plus(1:3,1:3)));
    sigma3_vel(i) = 3*sqrt(trace(P_plus(4:6,4:6)));

    % Error estimate
    [r_real,v_real] = TLE_2_ECI(satrec,sat_epoch_et,vis_window_act(i)); 
    err_pos(i) = norm(r_real - x_mean_abs_Mango(1:3,i));
    err_vel(i) = norm(v_real - x_mean_abs_Mango(4:6,i));
end

% PLOT OF 3*SIGMA AND ERROR ESTIMATE
figure
subplot(1,2,1)
plot(vis_window_act/cspice_spd(),sigma3_pos,'LineWidth',2)
hold on
plot(vis_window_act/cspice_spd(),err_pos,'LineWidth',2)
xlim([vis_window_act(1)/cspice_spd() vis_window_act(end)/cspice_spd()])
grid on
title('UKF solution for position','Interpreter','latex')
xlabel('Epoch [MJD2000]','Interpreter','latex')
ylabel('$3\sigma$,Error [km]','Interpreter','latex')
legend('$3\sigma$','Error estimate','Interpreter','latex')
set(gca,'FontSize',18)
subplot(1,2,2)
plot(vis_window_act/cspice_spd(),sigma3_vel,'LineWidth',2)
hold on
plot(vis_window_act/cspice_spd(),err_vel,'LineWidth',2)
xlim([vis_window_act(1)/cspice_spd() vis_window_act(end)/cspice_spd()])
grid on
title('UKF solution for velocity','Interpreter','latex')
xlabel('Epoch [MJD2000]','Interpreter','latex')
ylabel('$3\sigma$,Error [km/s]','Interpreter','latex')
legend('$3\sigma$','Error estimate','Interpreter','latex','Location','best')
set(gca,'FontSize',18)

%% POINT 2: ESTIMATE TANGO RELATIVE STATE
% ------- SUB_POINT a) RELATIVE STATE AT t0 ------- %
% Mango absolute state in ECI (at t0)
[r_eci,v_eci] = TLE_2_ECI(satrec,sat_epoch_et,t0_et);
x_t0_mango_ECI = [r_eci; v_eci];

% Load TANGO TLE file and create "satrec" structure
TLE_file = 'tle\36827.3le';
satrec = read_3LE(TANGO.ID,TLE_file, whichconst);
% Compute the epoch
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

% Tango absolute state in ECI (at t0)
[r_eci,v_eci] = TLE_2_ECI(satrec,sat_epoch_et,t0_et);
x_t0_tango_ECI = [r_eci; v_eci];

% Computation of the rotation matrix
rotmat_eci2lvlh = ECI_2_LVLH(x_t0_mango_ECI);

% Tango relative state in ECI frame (at t_0)
x_t0_rel_ECI = x_t0_tango_ECI - x_t0_mango_ECI;

% Tango relative state in LVLH frame (at t_0)
x_t0_rel_LVLH = rotmat_eci2lvlh*x_t0_rel_ECI;

% ------- SUB-POINT b) MEASUREMENTS SIMULATION WITH CW EQUATIONS ------- %
% CW equations integration
R = norm(x_t0_mango_ECI(1:3));
n_ang = sqrt(GM/(R^3));
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[~, xx] = ode78(@(t,x) CW_rhs(t,x,n_ang),t_g,x_t0_rel_LVLH, options);
Mean_state_rel = xx';
% Measurements Covariance matrix
R_new = diag([station.noise_Az_El_new^2 ; ...
              station.noise_Az_El_new^2; ...
              station.noise_range_new^2]);
% Measurements simulation
Meas_Az_rel = zeros(1,length(t_g)); % Azimuth
Meas_El_rel = zeros(1,length(t_g)); % Elevation
Meas_Range_rel = zeros(1,length(t_g)); % Range
for i = 1 : length(t_g)
    % Expected measurements
    [Range,Az,El] = cspice_reclat(Mean_state_rel(1:3,i));
    y_hat = [Az,El,Range]';
    % Actual measurements
    Measurements = mvnrnd(y_hat,R_new);
    Meas_Az_rel(i) = Measurements(1);
    Meas_El_rel(i) = Measurements(2);
    Meas_Range_rel(i) = Measurements(3);
end

% ------- SUB-POINT c) TANGO RELATIVE STATE ESTIMATE WITH UKF ------- %
% Time grid definition
index = pos_new(end);
t_grid_in = t_g(index+1); % Starting date of the time grid
t_grid_fin = t_grid_in + 20*60; % Final date of the time grid
t_grid = t_grid_in:5:t_grid_fin;
len = length(t_grid);

% Variables inizialization
x_mean_rel_Tango = zeros(6,len);
P_rel_Tango = cell(1,len);
sigma3_pos = zeros(1,len);
sigma3_vel = zeros(1,len);

% Starting the propagation from the last date of the actual visibility window
t_in = t_g(index); % Time at step 0
x_mean_plus = Mean_state_rel(:,index);
P0_new = diag([0.01,0.01,0.1,0.0001,0.0001,0.001]); % Covariance matrix at t_in
P_plus = P0_new;
for i = 1 : len
    t_fin = t_grid(i); % Time at step k
    y = [Meas_Az_rel(index+i) Meas_El_rel(index+i) Meas_Range_rel(index+i)]'; % Actual measurements
    [x_mean_rel_Tango(:,i),P_rel_Tango{i}] = UKF(x_mean_plus,P_plus,y,R_new,t_in,t_fin,n_ang,par,1);

    % Update intial conditions for the following step
    x_mean_plus = x_mean_rel_Tango(:,i);
    P_plus = P_rel_Tango{i};
    t_in = t_fin; % Time at step k-1

    % 3sigma of the estimated covariavce
    sigma3_pos(i) = 3*sqrt(trace(P_plus(1:3,1:3)));
    sigma3_vel(i) = 3*sqrt(trace(P_plus(4:6,4:6)));
end
% Error estimate
r_real = Mean_state_rel(1:3,index+1:index+len);
v_real = Mean_state_rel(4:6,index+1:index+len);
err_pos = vecnorm(r_real - x_mean_rel_Tango(1:3,:));
err_vel = vecnorm(v_real - x_mean_rel_Tango(4:6,:));

% PLOT OF 3*SIGMA AND ERROR ESTIMATE
figure
subplot(1,2,1)
plot(t_grid/cspice_spd(),sigma3_pos,'LineWidth',2)
hold on
plot(t_grid/cspice_spd(),err_pos,'LineWidth',2)
xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
grid on
title('UKF solution for position','Interpreter','latex')
xlabel('Epoch [MJD2000]','Interpreter','latex')
ylabel('$3\sigma$,Error [km]','Interpreter','latex')
legend('$3\sigma$','Error estimate','Interpreter','latex')
set(gca,'FontSize',20)
subplot(1,2,2)
plot(t_grid/cspice_spd(),sigma3_vel,'LineWidth',2)
hold on
plot(t_grid/cspice_spd(),err_vel,'LineWidth',2)
xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
grid on
title('UKF solution for velocity','Interpreter','latex')
xlabel('Epoch [MJD2000]','Interpreter','latex')
ylabel('$3\sigma$,Error [km/s]','Interpreter','latex')
legend('$3\sigma$','Error estimate','Interpreter','latex')
set(gca,'FontSize',20)

%% POINT 3: RECONSTRUCT TANGO ABSOLUTE COVARIANCE

% Variables inizialization
sigma3_pos = zeros(1,len);
sigma3_vel = zeros(1,len);

% Initial conditions for UT propagation
P_abs_Mango = P_abs_Mango{end};
x_mean_abs_Mango = x_mean_abs_Mango(:,end);
t_in = vis_window_act(end);

for i = 1 : len
    t_fin = t_grid(i); % Update the final time

    % ------- SUB-POINT a) PROPAGATE MANGO ESTIMATED COVARIANCE ------- %
    [x_mean_abs_Mango,P_abs_Mango] = UT(par,x_mean_abs_Mango,P_abs_Mango,t_in,t_fin);

    % ------- SUB-POINT b) ROTATE TANGO ESTIMATED COVARIANCE ------- %
    [~,rotmat_lvlh2eci] = ECI_2_LVLH(x_mean_abs_Mango);
    P_Tango_rel2abs =  rotmat_lvlh2eci*P_rel_Tango{i}* rotmat_lvlh2eci';

    % ------- SUB-POINT C) SUM THE TWO CONTRIBUTIONS ------- %
    P_abs_Tango = P_abs_Mango + P_Tango_rel2abs;
    % 3sigma computation
    sigma3_pos(i) = 3*sqrt(trace(P_abs_Tango(1:3,1:3)));
    sigma3_vel(i) = 3*sqrt(trace(P_abs_Tango(4:6,4:6)));
    
    t_in = t_fin; % Update the initial time
end

% PLOT OF 3*SIGMA AND ESTIMATE ERROR
figure
subplot(1,2,1)
plot(t_grid/cspice_spd(),sigma3_pos,'LineWidth',2)
xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
grid on
title('Solution for position','Interpreter','latex')
xlabel('Epoch [MJD2000]','Interpreter','latex')
ylabel('$3\sigma$[km]','Interpreter','latex')
legend('$3\sigma$','Interpreter','latex')
set(gca,'FontSize',20)
subplot(1,2,2)
plot(t_grid/cspice_spd(),sigma3_vel,'LineWidth',2)
xlim([t_grid(1)/cspice_spd() t_grid(end)/cspice_spd()])
grid on
title('Solution for velocity','Interpreter','latex')
xlabel('Epoch [MJD2000]','Interpreter','latex')
ylabel('$3\sigma$[km/s]','Interpreter','latex')
legend('$3\sigma$','Interpreter','latex')
set(gca,'FontSize',20)

Total_CPU_time = toc;
fprintf('\nOverall CPU time = %f [s]\n',Total_CPU_time)
% CLEAR KERNELS
cspice_kclear()
%% FUNCTIONS
function [xf,tf,xx, tt]  = propagate(t0,x0,tf,attractor,type)
% DESCRIPTION
% This function propagates the state from the initial condition x0 for a 3D 2BP

% INPUTS
% t0 ---- [1x1] Initial propagation time
% x0 ---- [6x1] Initial state
% tf ---- [1x1] Final propagation time
% attractor ---- 'Name', Name of the attractor
% type --- [1x1] type=0 -> keplerian propagator 
%                type=1 -> keplerian propagator + J2

% OUTPUTS 
% xf ---- [6x1] Final state
% tf ---- [1x1] Final date
% xx ---- Result of the ODE at each step
% tt ---- All the time steps

    if nargin == 4
            type = 0;
    end

    % Initialize propagation data
    GM = cspice_bodvrd(attractor,'GM',1);
    
    % Perform integration
    options = odeset('reltol', 1e-12, 'abstol', 1e-12);
    [tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM,type), [t0 tf], x0, options);

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


function dxdt = CW_rhs(~,x,n)
% DESCRIPTION
% The function evaluates the right-hand-side of a Clohessy Wiltshire propagator.

% INPUTS:
% t --- [1x1] Epoch (used only for type=1)
% x --- [6x1] Cartesian state vector wrt attractor Barycentre
% n --- [1x1] Mean motion of the lvlh frame

% OUTPUTS:
% dxdt --- [6x1] RHS

    % Initialize right-hand-side
    dxdt = zeros(6,1);
    
    % Extract position and velocity
    rr = x(1:3);
    vv = x(4:6);
    
    % Position detivative is object's velocity
    dxdt(1:3) = x(4:6);   
    % Compute the acceleration using CW formula
    dxdt(4) = 3*n^2*rr(1) + 2*n*vv(2);
    dxdt(5) = -2*n*vv(1);
    dxdt(6) = -n^2*rr(3);
end


function [meas_Az,meas_El,meas_Range,i_visibility] = TLE_2_measurements(satrec,TLE_epoch_et,station,Vis_window,R)
% DESCRIPTION
% Given the TLE evaluation, the function computes the actual
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


function [Az,El,Range] = measurements_mean_value(station_name,rr_sat_ECI,time_et)
% DESCRIPTION
% Given the position of the satellite the function computes Azimuth, 
% Elevation and Range.

% INPUTS
% station_name --- 'Name', Name of the ground station
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


function [y_mean,Py] = UT(par,x_mean,Px,t_in,t_fin)
% OBJECTIVE
% Propagation of mean value and covariance with unscented transform

% INPUTS
% par --- [1x3] [alpha,beta,k] parameters for UT computation
% x_mean --- [6x1] mean value of the previous step
% Px --- [6x6] Covariance of the previous step
% t_in --- [1x1] Previous time 
% t_fin --- [1x1] Current time 

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
    Y = zeros(n,2*n+1);
    
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
        Y(:,i) = propagate(t_in,X(:,i),t_fin,'Earth',1);
    end
    
    % WEIGHTED SAMPLE MEAN AND COVARIANCE COMPUTATION
    y_mean = W_m0*Y(:,1) + W_mi*sum(Y(:,2:end),2);
    Py = W_c0*(Y(:,1)-y_mean)*(Y(:,1)-y_mean)';
    for i = 2 : (2*n+1)
        Py = Py + W_ci*(Y(:,i)-y_mean)*(Y(:,i)-y_mean)';
    end
end


function [x_mean_plus,P_plus] = UKF(x0_mean_plus,P0_plus,y,R,t_in,t_fin,spec_var,par,type)
% DESCRIPTION
% The function performes a state estimate with UKF

% INPUTS
% x0_mean_plus --- [6x1] mean value of the previous step
% P0_plus --- [6x6] Covariance of the previous step
% y --- [3x1] Actual measurements (Azimuth,Elevation and Range) of the current step
% R --- [3x3] Measurements covariance matrix
% t_in --- [1x1] Previous time 
% t_fin --- [1x1] Current time 
% spec_var = station_name --- 'Name', Name of the ground station -> Propagation with absolute dynamics
%          = n --- [1x1] Mean motion of Mango
% par --- [1x3] [alpha,beta,k] parameters for UT computation -> Propagation with relative dynamics
% type --- [1x1] type=0 -> keplerian propagator 
%                type=1 -> keplerian propagator + J2

% OUTPUTS
% x_mean_plus --- [6x1] mean value of the current step
% P_plus --- [6x6] Covariance of the current step

    % Initialization
    alpha = par(1);
    beta = par(2);
    k = par(3);
    n = length(x0_mean_plus);
    c = alpha^2*(n+k);
    XX = zeros(n,2*n+1);
    delta_X = zeros(n,2*n);
    X = zeros(n,2*n+1);
    Y = zeros(3,2*n+1);
    
    
    % SIGMA POINTS COMPUTATION
    delta_X(:,1:n) = sqrtm(c*P0_plus);
    delta_X(:,n+1:2*n) = -sqrtm(c*P0_plus);
    XX(:,1) = x0_mean_plus;
    for i = 2 : (2*n+1)
        XX(:,i) = x0_mean_plus + delta_X(:,i-1);
    end
    
    % WEIGHTS COMPUTATION
    W_m0 = 1 - n/c;
    W_c0 = (2-alpha^2+beta) - n/c;
    W_mi = 1/(2*c);
    W_ci = W_mi;
    
    % SIGMA POINTS PROPAGATION 
    for i = 1 : (2*n+1)
        switch type
            case 0
            % Compute state samples
            X(:,i) = propagate(t_in,XX(:,i),t_fin,'Earth',1);
            % Compute measurements samples
            [Az,El,Range] = measurements_mean_value(spec_var,X(1:3,i),t_fin);
            Y(:,i) = [Az;El;Range];
            case 1
            % Compute state samples
            options = odeset('reltol', 1e-12, 'abstol', 1e-12);
            [~, xx] = ode78(@(t,x) CW_rhs(t,x,spec_var),[t_in t_fin],XX(:,i), options);
            X(:,i) = xx(end,:)';
            % Compute measurements samples
            [Range,Az,El] = cspice_reclat(X(1:3,i));
            Y(:,i) = [Az;El;Range];
        end
    end
    
    % WEIGHTED SAMPLE MEAN AND COVARIANCE COMPUTATION
    % State Covariance
    x_mean_minus = W_m0*X(:,1) + W_mi*sum(X(:,2:end),2);
    P_minus = W_c0*(X(:,1)-x_mean_minus)*(X(:,1)-x_mean_minus)';
    for i = 2 : (2*n+1)
        P_minus = P_minus + W_ci*(X(:,i)-x_mean_minus)*(X(:,i)-x_mean_minus)';
    end
    % Measurements Covariance
    y_mean_minus = W_m0*Y(:,1) + W_mi*sum(Y(:,2:end),2);
    diff_1 = [angdiff(y_mean_minus(1),Y(1,1)); 
             angdiff(y_mean_minus(2),Y(2,1));
                Y(3,1)-y_mean_minus(3)];
    Pee = (W_c0*diff_1)*diff_1';
    diff_i = cell(1,2*n);
    for i = 2 : (2*n+1)
        diff_i{i-1} = [angdiff(y_mean_minus(1),Y(1,i)); 
                    angdiff(y_mean_minus(2),Y(2,i));
                    Y(3,i)-y_mean_minus(3)];
        Pee = Pee + (W_ci*diff_i{i-1})*(diff_i{i-1})';
    end
    Pee = Pee + R;
    % Cross-covariance
    Pxy = W_c0*(X(:,1)-x_mean_minus)*diff_1';
    for i = 2 : (2*n+1)
        Pxy = Pxy + W_ci*(X(:,i)-x_mean_minus)*(diff_i{i-1})';
    end
    
    % UPDATE OF MEAN VALUE AND COVARIANCE
    K = Pxy/Pee;
    diff = [angdiff(y_mean_minus(1),y(1)); 
            angdiff(y_mean_minus(2),y(2));
            y(3)-y_mean_minus(3)];
    x_mean_plus = x_mean_minus + K*diff;
    P_plus = P_minus - K*Pee*K';
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


function [rotmat_eci2lvlh,rotmat_lvlh2eci] = ECI_2_LVLH(xx_ref_eci)
% DESCRIPTION
% The function computes the rotation matrixes from ECI to LVLH and from
% LVLH to ECI

% INPUT
% xx_ref_eci --- [6x1] State of the body associated to lvlh frame

% OUTPUTS
% rotmat_eci2lvlh --- [6x6] rotation matrix from ECI to LVLH
% rotmat_lvlh2eci --- [6x6] rotation matrix from LVLH to ECI


% Position and velocity in eci
rr_ref_eci = xx_ref_eci(1:3);
r_ref_eci = norm(rr_ref_eci);
vv_ref_eci = xx_ref_eci(4:6);

    % Position rotation matrix
    i = rr_ref_eci/r_ref_eci;
    k = cross(rr_ref_eci,vv_ref_eci)/norm(cross(rr_ref_eci,vv_ref_eci));
    j = cross(k,i);
    R = [i';j';k'];
    
    % Derivative of position rotation matrix
    di = 1/r_ref_eci*(vv_ref_eci - (dot(i,vv_ref_eci))*i);
    dj = cross(k,di);
    dk = zeros(3,1);
    R_dot = [di';dj';dk'];
    
    % Output rotation matrix from eci to lvlh
    rotmat_eci2lvlh = [R      zeros(3);
                        R_dot    R];
    
    rotmat_lvlh2eci = [R'      zeros(3);
                        R_dot'    R'];
end




