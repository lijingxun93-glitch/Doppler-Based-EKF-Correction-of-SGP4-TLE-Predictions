clear;
measured_data_1 = readtable('C:\Users\22743\Desktop\Master project\satellite data\ISS 20260622\epoch 20260516 SGP4_Az_El_Doppler_Visible.csv');
measured_data_2 = readtable('C:\Users\22743\Desktop\Master project\satellite data\ISS 20260622\epoch 20260621 SGP4_Az_El_Doppler_Visible.csv');
% AWGN
sigma_D = 10;
rng(1);
measured_data_2.Doppler_Hz = measured_data_2.Doppler_Hz + sigma_D * randn(size(measured_data_2.Doppler_Hz));
% Use data_2 as reference, interpolate data_2 onto data_1 time
t1 = measured_data_1.Time_UTC;
t2 = measured_data_2.Time_UTC;

idx2_vis = measured_data_2.Elevation_deg > 0;

t2_vis = t2(idx2_vis);

doppler2_vis = measured_data_2.Doppler_Hz(idx2_vis);
az2_vis      = measured_data_2.Azimuth_deg(idx2_vis);
el2_vis      = measured_data_2.Elevation_deg(idx2_vis);

% Create NaN arrays first, so times outside visible interval will stay empty
doppler_ref = NaN(size(t1));
az_ref      = NaN(size(t1));
el_ref      = NaN(size(t1));

% Only interpolate t1 inside the visible time range of measured_data_2
idx_interp = t1 >= min(t2_vis) & t1 <= max(t2_vis);

doppler_ref(idx_interp) = interp1(t2_vis, doppler2_vis, t1(idx_interp), 'linear');
az_ref(idx_interp)      = interp1(t2_vis, az2_vis,      t1(idx_interp), 'linear');
el_ref(idx_interp)      = interp1(t2_vis, el2_vis,      t1(idx_interp), 'linear');

doppler_error = measured_data_1.Doppler_Hz - doppler_ref;
az_error      = wrapTo180_local(measured_data_1.Azimuth_deg - az_ref);
el_error      = measured_data_1.Elevation_deg - el_ref;

% Mean squared error over valid samples only
idx_doppler_error = isfinite(doppler_error);
idx_az_error      = isfinite(az_error);
idx_el_error      = isfinite(el_error);

doppler_rmse = sqrt(mean(doppler_error(idx_doppler_error).^2));
az_rmse      = sqrt(mean(az_error(idx_az_error).^2));
el_rmse      = sqrt(mean(el_error(idx_el_error).^2));

fig1 = figure(1);

subplot(3,2,1);
plot(measured_data_1.Time_UTC, measured_data_1.Doppler_Hz);
hold on; grid on;
plot(measured_data_2.Time_UTC, measured_data_2.Doppler_Hz);
xlabel('Time (UTC)');
ylabel('Doppler shift (Hz)');
title('Doppler shifts');
legend('TLE epoch at 2026-05-16', 'TLE epoch at 2026-06-21');
hold off;

subplot(3,2,2);
plot(t1, doppler_error);
grid on;
xlabel('Time (UTC)');
ylabel('Doppler residual (Hz)');
title(sprintf('Doppler residual (RMSE = %.3f Hz)', doppler_rmse));

subplot(3,2,3);
plot(measured_data_1.Time_UTC, measured_data_1.Azimuth_deg);
hold on; grid on;
plot(measured_data_2.Time_UTC, measured_data_2.Azimuth_deg);
xlabel('Time (UTC)');
ylabel('Azimuth (deg)');
title('Az angles');
legend('TLE epoch at 2026-05-16', 'TLE epoch at 2026-06-21');
hold off;

subplot(3,2,4);
plot(t1, az_error);
grid on;
xlabel('Time (UTC)');
ylabel('Azimuth error (deg)');
title(sprintf('Azimuth Error (RMSE = %.6f deg)', az_rmse));

subplot(3,2,5);
plot(measured_data_1.Time_UTC, measured_data_1.Elevation_deg);
hold on; grid on;
plot(measured_data_2.Time_UTC, measured_data_2.Elevation_deg);
xlabel('Time (UTC)');
ylabel('Elevation (deg)');
title('El angles');
legend('TLE epoch at 2026-05-16', 'TLE epoch at 2026-06-21');
hold off;

subplot(3,2,6);
plot(t1, el_error);
grid on;
xlabel('Time (UTC)');
ylabel('Elevation error (deg)');
title(sprintf('Elevation Error (RMSE = %.6f deg)', el_rmse));

format_datetime_axes_english(fig1);

function a = wrapTo180_local(a)
% Wrap angle in degrees to [-180, 180].
    a = mod(a + 180, 360) - 180;
end
%%
measured_data = measured_data_2;

t_real = measured_data.Time_UTC;
idx_real = find(measured_data.Elevation_deg > 0);



%% ------------------------------------Model initialization-----------------------
% 状态量删除doppler
format long g

global const
SAT_Const

ge = 398600.8; % Earth gravitational constant [km3/s2]
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);

% TLE file name
fname = 'ISS.txt';

% Open the TLE file and read TLE elements
fid = fopen(fname, 'r');

% read first line
tline = fgetl(fid);
Cnum = tline(3:7);
SC   = tline(8);
ID   = tline(10:17);
year = str2num(tline(19:20));
doy  = str2num(tline(21:32));
epoch = str2num(tline(19:32));
TD1   = str2num(tline(34:43));
TD2   = str2num(tline(45:50));
ExTD2 = str2num(tline(51:52));
BStar = str2num(tline(54:59));
ExBStar = str2num(tline(60:61));
BStar = BStar*1e-5*10^ExBStar;
Etype = tline(63);
Enum  = str2num(tline(65:end));

% read second line
tline = fgetl(fid);
i = str2num(tline(9:16));
raan = str2num(tline(18:25));
e = str2num(strcat('0.',tline(27:33)));
omega = str2num(tline(35:42));
M = str2num(tline(44:51));
no = str2num(tline(53:63));
a = ( ge/(no*2*pi/86400)^2 )^(1/3);
rNo = str2num(tline(65:end));

fclose(fid);

satdata.epoch             = epoch;
satdata.norad_number      = Cnum;
satdata.bulletin_number   = ID;
satdata.classification    = SC;
satdata.revolution_number = rNo;
satdata.ephemeris_type    = Etype;
satdata.xmo               = M * (pi/180);
satdata.xnodeo            = raan * (pi/180);
satdata.omegao            = omega * (pi/180);
satdata.xincl             = i * (pi/180);
satdata.eo                = e;
satdata.xno               = no * TWOPI / MINUTES_PER_DAY;
satdata.xndt2o            = TD1 * TWOPI / MINUTES_PER_DAY_SQUARED;
satdata.xndd6o            = TD2 * 10^ExTD2 * TWOPI / MINUTES_PER_DAY_CUBED;
satdata.bstar             = BStar;

satdata_err = satdata;

TimeShift = measured_data_2.Time_UTC(1)-measured_data_1.Time_UTC(1);
if abs(minutes(TimeShift)) > 1
    satdata_err.epoch = satdata_err.epoch + minutes(TimeShift)/60/24;
    doy = doy + minutes(TimeShift)/60/24;
end


% read Earth orientation parameters
fid = fopen('EOP-All.txt','r');
while ~feof(fid)
    tline = fgetl(fid);
    ktmp = strfind(tline,'NUM_OBSERVED_POINTS');
    if (ktmp == 1)
        numrecsobs = str2num(tline(21:end));
        tline = fgetl(fid);
        for ii=1:numrecsobs
            eopdata(:,ii) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        for ii=1:4
            tline = fgetl(fid);
        end
        numrecspred = str2num(tline(22:end));
        tline = fgetl(fid);
        for ii=numrecsobs+1:numrecsobs+numrecspred
            eopdata(:,ii) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        break
    end
end
fclose(fid);

if (year < 57)
    year = year + 2000;
else
    year = year + 1900;
end
[mon,day,hr,minute,sec] = days2mdh(year,doy);
MJD_Epoch = Mjday(year,mon,day,hr,minute,sec);

% ------ Ground station
gs.lat_deg  = 35.6833; %Tehran
gs.lon_deg  = 51.35;
gs.alt_m    = 1191;

% ------ choose model start time
startTime = datetime(2026,6,22,5,43,0);
t_epoch   = datetime(year,mon,day,hr,minute,sec);
t_start   = startTime - t_epoch;   % duration

duration_min = floor(seconds(t_real(end) - startTime)/60);
dt = 1;   % second
t = (0:duration_min*60)/60;   % in minutes for sgp4 tsince

fc = 100e6;
c = 299792.458;   % km/s

N = length(t);
dopplerShift = zeros(N,1);
azimuth      = zeros(N,1);
elevation    = zeros(N,1);

doppler_error = nan(N,1);        % measured - predicted, Hz
azimuth_error = nan(N,1);        % measured - predicted, °
elevation_error = nan(N,1);      % measured - predicted, °
doppler_meas_match = nan(N,1);   % matched measured Doppler, Hz
azimuth_meas_match = nan(N,1);   % matched measured Azimuth, Hz
elevation_meas_match = nan(N,1); % matched measured Elevation, Hz
update_flag = false(N,1);        % whether EKF update is used

% ---------- EKF ----------
% State:
% x = [dM, dRAAN, dARGP, dINC, de, dn, dBstar]^T
% Units:
% dM, dRAAN, dARGP, dINC : rad
% de                      : dimensionless
% dn                      : rad/min
% dBstar                  : same unit as satdata.bstar

x = zeros(7,1);

% Random-walk model for TLE parameter correction
% x_k = x_{k-1} + w
F = eye(7);

% Initial uncertainty P

P = diag([ ...
    deg2rad(5), ...
    deg2rad(2), ...
    deg2rad(2), ...
    deg2rad(0.5), ...
    1e-3, ...
    satdata.xno*5e-5, ...
    max(abs(satdata.bstar),2e-8) ...
    ].^2);
% Process noise Q
% 先设小一点，避免 Doppler-only 过度改轨道
Q = diag([ ...
    deg2rad(1e-5), ...       % dM
    deg2rad(1e-7), ...       % dRAAN
    deg2rad(1e-7), ...       % dARGP
    deg2rad(1e-8), ...       % dINC
    1e-8, ...                % de
    satdata.xno*1e-9, ...    % dn
    max(abs(satdata.bstar),1e-8)*1e-4 ...
].^2);

% Doppler measurement noise
R = 20^2 ;   % Hz^2, start with 10 Hz std

j = 1;
for tsince = t

    % current model time
    t_model = startTime + seconds((j-1)*dt);

    % ---------- time update ----------
    x = F * x;
    P = F * P * F' + Q;

    % ---------- UTC / EOP ----------
    MJD_UTC = MJD_Epoch + (minutes(t_start) + tsince)/1440;

    [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
    [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    MJD_UT1 = MJD_UTC + UT1_UTC/86400;
    MJD_TT  = MJD_UTC + TT_UTC/86400;
    T = (MJD_TT-const.MJD_J2000)/36525;

    % ---------- nominal orbit from SGP4 ----------
    % [rteme, vteme] = sgp4(minutes(t_start) + tsince, satdata);
    % [recef, vecef] = teme2ecef(rteme, vteme, T, MJD_UT1+2400000.5, LOD, x_pole, y_pole, 2);

    % ---------- SGP4 with EKF-corrected TLE parameters ----------
    tsince_total = minutes(t_start) + tsince;
    
    satdata_est = apply_tle_state_correction(satdata_err, x);
    
    [rteme_est, vteme_est] = sgp4(tsince_total, satdata_est);
    [recef_est, vecef_est] = teme2ecef(rteme_est, vteme_est, T, MJD_UT1+2400000.5, LOD, x_pole, y_pole, 2);
    
    % ---------- ground station ----------
    % [rsite_ecef, vsite_ecef] = groundstation_ecef(gs.lat_deg, gs.lon_deg, gs.alt_m/1000);
    [rsite_ecef, ~] = groundstation_ecef(gs.lat_deg, gs.lon_deg, gs.alt_m/1000);
    vsite_ecef = [0; 0; 0];   % fixed ground station in ECEF frame
    % ---------- corrected state from corrected TLE ----------
    r_est = recef_est;
    v_est = vecef_est;

    rho_ecef    = r_est - rsite_ecef;
    rhodot_ecef = v_est - vsite_ecef;

    [enu, enu_dot] = ecef2enu_vector(rho_ecef, rhodot_ecef, gs.lat_deg, gs.lon_deg);

    E = enu(1);
    Nn = enu(2);
    U = enu(3);

    range_km = norm(enu);

    el = asin(U / range_km);
    az = atan2(E, Nn);
    if az < 0
        az = az + 2*pi;
    end

    elevation(j) = rad2deg(el);
    azimuth(j)   = rad2deg(az);

    rho_unit = rho_ecef / norm(rho_ecef);
    range_rate = dot(rhodot_ecef, rho_unit);   % km/s
    dopplerShift(j) = - (range_rate / c) * fc;

    % ---------- find nearest measurement ----------
    [dt_meas,k_meas] = min(abs(t_real - t_model));

    % ---------- record error before / after EKF output ----------
    doppler_meas_match(j) = measured_data.Doppler_Hz(k_meas);
    azimuth_meas_match(j) = measured_data.Azimuth_deg(k_meas);
    elevation_meas_match(j) = measured_data.Elevation_deg(k_meas);

    if dt_meas<seconds(0.5) && measured_data.Elevation_deg(k_meas) > 0
        doppler_error(j) = doppler_meas_match(j) - dopplerShift(j);
        azimuth_error(j) = azimuth_meas_match(j) - azimuth(j);
        elevation_error(j) = elevation_meas_match(j) - elevation(j);
    end

    % ---------- update only when satellite is above horizon in measured data ----------
    if dt_meas<seconds(0.5) && measured_data.Elevation_deg(k_meas) > 0

        z = measured_data.Doppler_Hz(k_meas);

        z_hat = dopplerShift(j);

        H = numerical_H_doppler_tle( ...
                x, satdata_err, tsince_total, ...
                T, MJD_UT1, LOD, x_pole, y_pole, ...
                rsite_ecef, vsite_ecef, fc, c);
        innov = z - z_hat;

        S = H * P * H' + R;
        K = (P * H') / S;

        x = x + K * innov;
        % P = (eye(7) - K * H) * P;
        I7 = eye(7);
        P = (I7 - K*H) * P * (I7 - K*H)' + K * R * K';
        % recompute after EKF update using corrected TLE
        satdata_est = apply_tle_state_correction(satdata_err, x);
        
        [rteme_est, vteme_est] = sgp4(tsince_total, satdata_est);
        [recef_est, vecef_est] = teme2ecef(rteme_est, vteme_est, T, MJD_UT1+2400000.5, LOD, x_pole, y_pole, 2);
        
        r_est = recef_est;
        v_est = vecef_est;

        rho_ecef    = r_est - rsite_ecef;
        rhodot_ecef = v_est - vsite_ecef;

        [enu, enu_dot] = ecef2enu_vector(rho_ecef, rhodot_ecef, gs.lat_deg, gs.lon_deg);

        E = enu(1);
        Nn = enu(2);
        U = enu(3);

        range_km = norm(enu);

        el = asin(U / range_km);
        az = atan2(E, Nn);
        if az < 0
            az = az + 2*pi;
        end

        elevation(j) = rad2deg(el);
        azimuth(j)   = rad2deg(az);

        rho_unit = rho_ecef / norm(rho_ecef);
        range_rate = dot(rhodot_ecef, rho_unit);
        dopplerShift(j) = - (range_rate / c) * fc;

        doppler_error(j) = doppler_meas_match(j) - dopplerShift(j);
        azimuth_error(j) = azimuth_meas_match(j) - azimuth(j);
        elevation_error(j) = elevation_meas_match(j) - elevation(j);

        update_flag(j) = true;
    end

    j = j + 1;
end
fig3 = figure(2);

% Mean squared error over valid above-horizon matched samples only
idx_doppler_error = isfinite(doppler_error);
idx_azimuth_error = isfinite(azimuth_error);
idx_elevation_error = isfinite(elevation_error);

doppler_rmse_ekf = sqrt(mean(doppler_error(idx_doppler_error).^2));
azimuth_rmse_ekf = sqrt(mean(azimuth_error(idx_azimuth_error).^2));
elevation_rmse_ekf = sqrt(mean(elevation_error(idx_elevation_error).^2));

idx_EKF = find(elevation>0); 
t_EKF = minutes(t)+startTime;
subplot(3,2,1);
plot(t_EKF(idx_EKF),dopplerShift(idx_EKF));
hold on;

subplot(3,2,2);
idx_valid = ~isnan(doppler_error);
plot(minutes(t(idx_valid))+startTime, doppler_error(idx_valid), 'LineWidth', 1.2);
grid on;
xlabel('Time (UTC)');
ylabel('Doppler residual (Hz)');
title(sprintf('Doppler residual (RMSE = %.3f Hz)', doppler_rmse_ekf));
yline(0,'--');

subplot(3,2,3);
plot(t_EKF(idx_EKF),azimuth(idx_EKF));
hold on;

subplot(3,2,4);
idx_valid = ~isnan(azimuth_error);
plot(minutes(t(idx_valid))+startTime, azimuth_error(idx_valid), 'LineWidth', 1.2);
grid on;
xlabel('Time (UTC)');
ylabel('Azimuth error (°)');
title(sprintf('Azimuth Error (RMSE = %.6f deg)', azimuth_rmse_ekf));
yline(0,'--');

subplot(3,2,5);
plot(t_EKF(idx_EKF),elevation(idx_EKF));
hold on;

subplot(3,2,6);
idx_valid = ~isnan(elevation_error);
plot(minutes(t(idx_valid))+startTime, elevation_error(idx_valid), 'LineWidth', 1.2);
grid on;
xlabel('Time (UTC)');
ylabel('Elevation error (°)');
title(sprintf('Elevation Error (RMSE = %.6f deg)', elevation_rmse_ekf));
yline(0,'--');

subplot(3,2,1);
plot(t_real(idx_real), measured_data.Doppler_Hz(idx_real)); 
legend('TLE epoch at 2026-05-16', 'TLE epoch at 2026-06-21', 'Location','best'); 
title('Doppler shifts');
grid on; 
ylabel('Doppler (Hz)');
hold off;

subplot(3,2,3);
plot(t_real(idx_real), measured_data.Azimuth_deg(idx_real));  

legend('TLE epoch at 2026-05-16', 'TLE epoch at 2026-06-21', 'Location','best'); 
title('Azimuth');
grid on; 
ylabel('Angles (deg)');
title('Az angles');
hold off;

subplot(3,2,5);
plot(t_real(idx_real), measured_data.Elevation_deg(idx_real));
grid on; 
ylabel('Angles (deg)');
title('El angles');
legend('TLE epoch at 2026-05-16', 'TLE epoch at 2026-06-21', 'Location','best'); 
hold off;

format_datetime_axes_english(fig3);
function satdata_new = apply_tle_state_correction(satdata_base, x)
% Apply EKF correction to TLE/SGP4 mean elements.
%
% x = [dM, dRAAN, dARGP, dINC, de, dn, dBstar]^T

    satdata_new = satdata_base;

    satdata_new.xmo    = wrapTo2Pi_local(satdata_base.xmo    + x(1));
    satdata_new.xnodeo = wrapTo2Pi_local(satdata_base.xnodeo + x(2));
    satdata_new.omegao = wrapTo2Pi_local(satdata_base.omegao + x(3));
    satdata_new.xincl  = satdata_base.xincl + x(4);

    % Keep inclination physically valid
    satdata_new.xincl = max(min(satdata_new.xincl, pi), 0);

    % Keep eccentricity valid for SGP4
    satdata_new.eo = satdata_base.eo + x(5);
    satdata_new.eo = max(min(satdata_new.eo, 0.99), 1e-8);

    % Mean motion must be positive
    satdata_new.xno = satdata_base.xno + x(6);
    satdata_new.xno = max(satdata_new.xno, 1e-10);

    satdata_new.bstar = satdata_base.bstar + x(7);
end


function z_hat = meas_model_doppler_tle( ...
        x, satdata_base, tsince_total, ...
        T, MJD_UT1, LOD, x_pole, y_pole, ...
        rsite_ecef, vsite_ecef, fc, c)
% Doppler measurement model using corrected TLE parameters.

    satdata_est = apply_tle_state_correction(satdata_base, x);

    [rteme, vteme] = sgp4(tsince_total, satdata_est);
    [recef, vecef] = teme2ecef(rteme, vteme, T, MJD_UT1+2400000.5, LOD, x_pole, y_pole, 2);

    rho_ecef    = recef - rsite_ecef;
    rhodot_ecef = vecef - vsite_ecef;

    rho_unit = rho_ecef / norm(rho_ecef);
    range_rate = dot(rhodot_ecef, rho_unit);   % km/s

    z_hat = - (range_rate / c) * fc;           % Hz
end


function H = numerical_H_doppler_tle( ...
        x, satdata_base, tsince_total, ...
        T, MJD_UT1, LOD, x_pole, y_pole, ...
        rsite_ecef, vsite_ecef, fc, c)
% Numerical Jacobian:
% H = d(Doppler) / d([dM,dRAAN,dARGP,dINC,de,dn,dBstar])

    h0 = meas_model_doppler_tle( ...
            x, satdata_base, tsince_total, ...
            T, MJD_UT1, LOD, x_pole, y_pole, ...
            rsite_ecef, vsite_ecef, fc, c);

    n_state = length(x);
    H = zeros(1,n_state);

    % Finite difference steps
    dx_step = [ ...
        deg2rad(1e-5); ...                         % dM
        deg2rad(1e-5); ...                         % dRAAN
        deg2rad(1e-5); ...                         % dARGP
        deg2rad(1e-5); ...                         % dINC
        1e-7; ...                                  % de
        max(abs(satdata_base.xno)*1e-8, 1e-12); ...% dn
        max(abs(satdata_base.bstar)*1e-4, 1e-10) ...% dBstar
    ];

    for ii = 1:n_state
        xp = x;
        xp(ii) = xp(ii) + dx_step(ii);

        hp = meas_model_doppler_tle( ...
                xp, satdata_base, tsince_total, ...
                T, MJD_UT1, LOD, x_pole, y_pole, ...
                rsite_ecef, vsite_ecef, fc, c);

        H(ii) = (hp - h0) / dx_step(ii);
    end
end


function angle_out = wrapTo2Pi_local(angle_in)
% Avoid Mapping Toolbox dependency.

    angle_out = mod(angle_in, 2*pi);

    if angle_out < 0
        angle_out = angle_out + 2*pi;
    end
end
function format_datetime_axes_english(fig)
% Show time at each tick and one English date per subplot.

    axes_list = findall(fig, 'Type', 'axes');

    for ii = 1:numel(axes_list)
        ax = axes_list(ii);
        
        if isa(ax.XAxis, 'matlab.graphics.axis.decorator.DatetimeRuler')
            % Set the DatetimeRuler property directly so every major tick
            % displays only the time, even if MATLAB previously selected an
            % automatic date-time format.
            ax.XAxis.TickLabelFormat = 'HH:mm';

            % Display the date once, using an English month abbreviation.
            date_label = string( ...
                ax.XLim(1), 'dd-MMM-yyyy', 'en_US');
            ax.XAxis.SecondaryLabel.String = char(date_label);
            ax.XAxis.SecondaryLabel.Interpreter = 'none';
            ax.XAxis.SecondaryLabel.Visible = 'on';
        end
    end
end
