clear; %clc; close all;
format long g

%% ===================== Measured data =====================
measured_data = readtable('C:\Users\22743\Desktop\Master project\satellite data\20251101\SGP4_Az_El_Doppler_All_Tabriz.csv');

t_real = measured_data.Time_UTC;
idx_real = find(measured_data.Elevation_deg > 0);

figure(12);

subplot(3,2,1);
plot(t_real(idx_real), measured_data.Doppler_Hz(idx_real)); 
grid on; hold on;
ylabel('Doppler (Hz)');
title('Doppler shift');

subplot(3,2,3);
plot(t_real(idx_real), measured_data.Azimuth_deg(idx_real));  
grid on; hold on;
ylabel('Angles (deg)');
title('Az angles');

subplot(3,2,5);
plot(t_real(idx_real), measured_data.Elevation_deg(idx_real));
grid on; hold on;
ylabel('Angles (deg)');
title('El angles');

%% ===================== Model initialization =====================
global const
SAT_Const

ge = 398600.8;                 % km^3/s^2
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = MINUTES_PER_DAY^2;
MINUTES_PER_DAY_CUBED = MINUTES_PER_DAY^3;

%% ===================== Read TLE =====================
fname = 'COSMOS_2428.txt';
fid = fopen(fname, 'r');

if fid < 0
    error("Cannot open TLE file: %s", fname);
end

% ----- Line 1 -----
tline = fgetl(fid);

Cnum   = tline(3:7);
SC     = tline(8);
ID     = tline(10:17);
year   = str2num(tline(19:20));
doy    = str2num(tline(21:32));
epoch  = str2num(tline(19:32));
TD1    = str2num(tline(34:43));
TD2    = str2num(tline(45:50));
ExTD2  = str2num(tline(51:52));
BStar  = str2num(tline(54:59));
ExBStar = str2num(tline(60:61));
BStar  = BStar * 1e-5 * 10^ExBStar;
Etype  = tline(63);
Enum   = str2num(tline(65:end)); %#ok<NASGU>

% ----- Line 2 -----
tline = fgetl(fid);

i      = str2num(tline(9:16));
raan   = str2num(tline(18:25));
e      = str2num(strcat('0.', tline(27:33)));
omega  = str2num(tline(35:42));
M      = str2num(tline(44:51));
no     = str2num(tline(53:63));
a      = (ge / (no*2*pi/86400)^2)^(1/3); %#ok<NASGU>
rNo    = str2num(tline(65:end));

fclose(fid);

satdata.epoch             = epoch;
satdata.norad_number      = Cnum;
satdata.bulletin_number   = ID;
satdata.classification    = SC;
satdata.revolution_number = rNo;
satdata.ephemeris_type    = Etype;
satdata.xmo               = M     * pi/180;
satdata.xnodeo            = raan  * pi/180;
satdata.omegao            = omega * pi/180;
satdata.xincl             = i     * pi/180;
satdata.eo                = e;
satdata.xno               = no * TWOPI / MINUTES_PER_DAY;
satdata.xndt2o            = TD1 * TWOPI / MINUTES_PER_DAY_SQUARED;
satdata.xndd6o            = TD2 * 10^ExTD2 * TWOPI / MINUTES_PER_DAY_CUBED;
satdata.bstar             = BStar;

%% ===================== Add initial TLE error =====================
satdata_err = satdata;

satdata_err.xmo    = satdata.xmo    + deg2rad(5);
satdata_err.xnodeo = satdata.xnodeo + deg2rad(2);
satdata_err.omegao = satdata.omegao + deg2rad(2);
satdata_err.xincl  = satdata.xincl  + deg2rad(0.5);
satdata_err.eo     = satdata.eo     + 1e-3;
satdata_err.xno    = satdata.xno    * (1 + 5e-5);
satdata_err.bstar  = satdata.bstar  * (1 + 1.0);

% satdata_err.xmo    = 1.01677646233434;
% satdata_err.xnodeo = 5.68788316992161;
% satdata_err.omegao = 5.26689928236506;
% satdata_err.xincl  = 1.7196153161247;
% satdata_err.eo     = 0.0009236;
% satdata_err.xno    = 0.0622665667143146;
% satdata_err.bstar  = 0.00010966;

%% ===================== Read EOP =====================
fid = fopen('EOP-All.txt','r');

if fid < 0
    error("Cannot open EOP-All.txt");
end

eopdata = [];

while ~feof(fid)
    tline = fgetl(fid);
    ktmp = strfind(tline, 'NUM_OBSERVED_POINTS');

    if ktmp == 1
        numrecsobs = str2num(tline(21:end));
        tline = fgetl(fid);

        for ii = 1:numrecsobs
            eopdata(:,ii) = fscanf(fid, '%i %d %d %i %f %f %f %f %f %f %f %f %i', [13 1]);
        end

        for ii = 1:4
            tline = fgetl(fid);
        end

        numrecspred = str2num(tline(22:end));
        tline = fgetl(fid);

        for ii = numrecsobs+1 : numrecsobs+numrecspred
            eopdata(:,ii) = fscanf(fid, '%i %d %d %i %f %f %f %f %f %f %f %f %i', [13 1]);
        end

        break
    end
end

fclose(fid);

%% ===================== Epoch =====================
if year < 57
    year = year + 2000;
else
    year = year + 1900;
end

[mon, day, hr, minute, sec] = days2mdh(year, doy);
MJD_Epoch = Mjday(year, mon, day, hr, minute, sec);

%% ===================== Ground station =====================
gs.lat_deg = 38.0833; %Tabriz
gs.lon_deg = 46.2833;
gs.alt_m   = 1361;

% Fixed ground-station position in ECEF, km.
% No ground-station velocity is needed because Doppler is calculated in ECEF.
rsite_ecef = groundstation_ecef(gs.lat_deg, gs.lon_deg, gs.alt_m/1000);

%% ===================== Time setting =====================
startTime = datetime(2025,10,31,13,19,11);
t_epoch   = datetime(year, mon, day, hr, minute, sec);
t_start   = startTime - t_epoch;

duration_min = floor(seconds(t_real(end) - startTime) / 60);
dt = 1;                         % second
t = (0:duration_min*60) / 60;   % minutes for SGP4

fc = 100e6;
c  = 299792.458;                % km/s

N = length(t);

dopplerShift = zeros(N,1);
azimuth      = zeros(N,1);
elevation    = zeros(N,1);

doppler_error   = nan(N,1);
azimuth_error   = nan(N,1);
elevation_error = nan(N,1);

doppler_meas_match   = nan(N,1);
azimuth_meas_match   = nan(N,1);
elevation_meas_match = nan(N,1);

update_flag = false(N,1);

%% ===================== EKF setting =====================
% State:
% x = [dM, dRAAN, dARGP, dINC, de, dn, dBstar]^T
%
% Units:
% dM, dRAAN, dARGP, dINC : rad
% de                      : dimensionless
% dn                      : rad/min
% dBstar                  : same unit as satdata.bstar

x = zeros(7,1);
F = eye(7);

P = diag([ ...
    deg2rad(5), ...
    deg2rad(2), ...
    deg2rad(2), ...
    deg2rad(0.5), ...
    1e-3, ...
    satdata.xno*5e-5, ...
    max(abs(satdata.bstar),2e-8) ...
].^2);

Q = diag([ ...
    deg2rad(1e-5), ...
    deg2rad(1e-7), ...
    deg2rad(1e-7), ...
    deg2rad(1e-8), ...
    1e-8, ...
    satdata.xno*1e-9, ...
    max(abs(satdata.bstar),1e-8)*1e-4 ...
].^2);

R = 20^2;      % Hz^2

%% ===================== Main loop =====================
for j = 1:N

    tsince = t(j);
    t_model = startTime + seconds((j-1)*dt);

    % ---------- EKF time update ----------
    x = F * x;
    P = F * P * F' + Q;

    % ---------- UTC / EOP ----------
    tsince_total = minutes(t_start) + tsince;
    MJD_UTC = MJD_Epoch + tsince_total / 1440;

    [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC] = ...
        IERS(eopdata, MJD_UTC, 'l'); %#ok<ASGLU>

    [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = ...
        timediff(UT1_UTC, TAI_UTC); %#ok<ASGLU>

    MJD_UT1 = MJD_UTC + UT1_UTC / 86400;
    MJD_TT  = MJD_UTC + TT_UTC  / 86400;
    T = (MJD_TT - const.MJD_J2000) / 36525;

    % ---------- Prediction before EKF update ----------
    [dopplerShift(j), azimuth(j), elevation(j)] = predict_obs_tle( ...
        x, satdata_err, tsince_total, ...
        T, MJD_UT1, LOD, x_pole, y_pole, ...
        rsite_ecef, gs.lat_deg, gs.lon_deg, fc, c);

    % ---------- Match nearest measurement ----------
    [~, k_meas] = min(abs(t_real - t_model));

    doppler_meas_match(j)   = measured_data.Doppler_Hz(k_meas);
    azimuth_meas_match(j)   = measured_data.Azimuth_deg(k_meas);
    elevation_meas_match(j) = measured_data.Elevation_deg(k_meas);

    if measured_data.Elevation_deg(k_meas) > 0
        doppler_error(j)   = doppler_meas_match(j)   - dopplerShift(j);
        azimuth_error(j)   = azimuth_meas_match(j)   - azimuth(j);
        elevation_error(j) = elevation_meas_match(j) - elevation(j);
    end

    % ---------- EKF measurement update ----------
    if measured_data.Elevation_deg(k_meas) > 0

        z = measured_data.Doppler_Hz(k_meas);
        z_hat = dopplerShift(j);
        innov = z - z_hat;

        H = numerical_H_doppler_tle( ...
            x, satdata_err, tsince_total, ...
            T, MJD_UT1, LOD, x_pole, y_pole, ...
            rsite_ecef, gs.lat_deg, gs.lon_deg, fc, c);

        S = H * P * H' + R;
        K = (P * H') / S;

        x = x + K * innov;

        I7 = eye(7);
        P = (I7 - K*H) * P * (I7 - K*H)' + K * R * K';

        % ---------- Recompute output after EKF update ----------
        [dopplerShift(j), azimuth(j), elevation(j)] = predict_obs_tle( ...
            x, satdata_err, tsince_total, ...
            T, MJD_UT1, LOD, x_pole, y_pole, ...
            rsite_ecef, gs.lat_deg, gs.lon_deg, fc, c);

        doppler_error(j)   = doppler_meas_match(j)   - dopplerShift(j);
        azimuth_error(j)   = azimuth_meas_match(j)   - azimuth(j);
        elevation_error(j) = elevation_meas_match(j) - elevation(j);

        update_flag(j) = true;
    end
end

%% ===================== Plot =====================
idx_EKF = find(elevation > 0);
t_EKF = minutes(t) + startTime;

subplot(3,2,1);
plot(t_EKF(idx_EKF), dopplerShift(idx_EKF));
legend('Doppler real', 'Doppler SGP4 & EKF', 'Location', 'best');
hold off;

subplot(3,2,2);
idx_valid = ~isnan(doppler_error);
plot(t_EKF(idx_valid), doppler_error(idx_valid), 'LineWidth', 1.2);
grid on;
xlabel('Time (UTC)');
ylabel('Doppler error (Hz)');
title('Doppler Error above horizon');
yline(0,'--');

subplot(3,2,3);
plot(t_EKF(idx_EKF), azimuth(idx_EKF));
legend('Az real', 'Az SGP4 & EKF', 'Location', 'best');
hold off;

subplot(3,2,4);
idx_valid = ~isnan(azimuth_error);
plot(t_EKF(idx_valid), azimuth_error(idx_valid), 'LineWidth', 1.2);
grid on;
xlabel('Time (UTC)');
ylabel('Azimuth error (deg)');
title('Azimuth Error above horizon');
yline(0,'--');

subplot(3,2,5);
plot(t_EKF(idx_EKF), elevation(idx_EKF));
legend('El real', 'El SGP4 & EKF', 'Location', 'best');
hold off;

subplot(3,2,6);
idx_valid = ~isnan(elevation_error);
plot(t_EKF(idx_valid), elevation_error(idx_valid), 'LineWidth', 1.2);
grid on;
xlabel('Time (UTC)');
ylabel('Elevation error (deg)');
title('Elevation Error above horizon');
yline(0,'--');

%% ===================== Local functions =====================

function [doppler_hz, az_deg, el_deg] = predict_obs_tle( ...
        x, satdata_base, tsince_total, ...
        T, MJD_UT1, LOD, x_pole, y_pole, ...
        rsite_ecef, lat_deg, lon_deg, fc, c)
% Predict Doppler, azimuth, and elevation from corrected TLE.
%
% Important:
%   recef, vecef are assumed to be ECEF-frame position and velocity.
%   For a fixed ground station in the ECEF frame:
%       rho_dot_ecef = vecef
%   No vsite_ecef is required.

    satdata_est = apply_tle_state_correction(satdata_base, x);

    [rteme, vteme] = sgp4(tsince_total, satdata_est);

    [recef, vecef] = teme2ecef( ...
        rteme, vteme, T, MJD_UT1 + 2400000.5, ...
        LOD, x_pole, y_pole, 2);

    rho_ecef    = recef - rsite_ecef;
    rhodot_ecef = vecef;

    [az_deg, el_deg] = azel_from_ecef(rho_ecef, lat_deg, lon_deg);

    rho_unit = rho_ecef / norm(rho_ecef);
    range_rate = dot(rhodot_ecef, rho_unit);     % km/s

    doppler_hz = - (range_rate / c) * fc;        % Hz
end


function H = numerical_H_doppler_tle( ...
        x, satdata_base, tsince_total, ...
        T, MJD_UT1, LOD, x_pole, y_pole, ...
        rsite_ecef, lat_deg, lon_deg, fc, c)
% Numerical Jacobian:
% H = d(Doppler) / d([dM,dRAAN,dARGP,dINC,de,dn,dBstar])

    h0 = meas_model_doppler_tle( ...
        x, satdata_base, tsince_total, ...
        T, MJD_UT1, LOD, x_pole, y_pole, ...
        rsite_ecef, lat_deg, lon_deg, fc, c);

    n_state = length(x);
    H = zeros(1,n_state);

    dx_step = [ ...
        deg2rad(1e-5); ...
        deg2rad(1e-5); ...
        deg2rad(1e-5); ...
        deg2rad(1e-5); ...
        1e-7; ...
        max(abs(satdata_base.xno)*1e-8, 1e-12); ...
        max(abs(satdata_base.bstar)*1e-4, 1e-10) ...
    ];

    for ii = 1:n_state
        xp = x;
        xp(ii) = xp(ii) + dx_step(ii);

        hp = meas_model_doppler_tle( ...
            xp, satdata_base, tsince_total, ...
            T, MJD_UT1, LOD, x_pole, y_pole, ...
            rsite_ecef, lat_deg, lon_deg, fc, c);

        H(ii) = (hp - h0) / dx_step(ii);
    end
end


function z_hat = meas_model_doppler_tle( ...
        x, satdata_base, tsince_total, ...
        T, MJD_UT1, LOD, x_pole, y_pole, ...
        rsite_ecef, lat_deg, lon_deg, fc, c)
% Doppler-only measurement model for EKF.

    [z_hat, ~, ~] = predict_obs_tle( ...
        x, satdata_base, tsince_total, ...
        T, MJD_UT1, LOD, x_pole, y_pole, ...
        rsite_ecef, lat_deg, lon_deg, fc, c);
end


function satdata_new = apply_tle_state_correction(satdata_base, x)
% Apply EKF correction to TLE/SGP4 mean elements.
%
% x = [dM, dRAAN, dARGP, dINC, de, dn, dBstar]^T

    satdata_new = satdata_base;

    satdata_new.xmo    = wrapTo2Pi_local(satdata_base.xmo    + x(1));
    satdata_new.xnodeo = wrapTo2Pi_local(satdata_base.xnodeo + x(2));
    satdata_new.omegao = wrapTo2Pi_local(satdata_base.omegao + x(3));
    satdata_new.xincl  = satdata_base.xincl + x(4);

    satdata_new.xincl = max(min(satdata_new.xincl, pi), 0);

    satdata_new.eo = satdata_base.eo + x(5);
    satdata_new.eo = max(min(satdata_new.eo, 0.99), 1e-8);

    satdata_new.xno = satdata_base.xno + x(6);
    satdata_new.xno = max(satdata_new.xno, 1e-10);

    satdata_new.bstar = satdata_base.bstar + x(7);
end


function rsite_ecef = groundstation_ecef(lat_deg, lon_deg, alt_km)
% Ground station ECEF position using WGS-84 ellipsoid.
% Output unit: km.
%
% In ECEF coordinates, the ground station is fixed.
% Therefore no ground-station velocity is returned.

    a = 6378.137;                       % km
    f = 1 / 298.257223563;
    e2 = f * (2 - f);

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    N = a / sqrt(1 - e2 * sin(lat)^2);

    x = (N + alt_km) * cos(lat) * cos(lon);
    y = (N + alt_km) * cos(lat) * sin(lon);
    z = (N * (1 - e2) + alt_km) * sin(lat);

    rsite_ecef = [x; y; z];
end


function [az_deg, el_deg] = azel_from_ecef(rho_ecef, lat_deg, lon_deg)
% Convert ECEF relative vector to azimuth/elevation.

    enu = ecef2enu_vector(rho_ecef, lat_deg, lon_deg);

    E = enu(1);
    N = enu(2);
    U = enu(3);

    range_km = norm(enu);

    el = asin(U / range_km);
    az = atan2(E, N);

    if az < 0
        az = az + 2*pi;
    end

    az_deg = rad2deg(az);
    el_deg = rad2deg(el);
end


function enu = ecef2enu_vector(vec_ecef, lat_deg, lon_deg)
% Rotate ECEF vector to local ENU vector.

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    R = [ ...
        -sin(lon),              cos(lon),             0; ...
        -sin(lat)*cos(lon),    -sin(lat)*sin(lon),    cos(lat); ...
         cos(lat)*cos(lon),     cos(lat)*sin(lon),    sin(lat) ...
    ];

    enu = R * vec_ecef;
end


function angle_out = wrapTo2Pi_local(angle_in)
% Avoid Mapping Toolbox dependency.

    angle_out = mod(angle_in, 2*pi);

    if angle_out < 0
        angle_out = angle_out + 2*pi;
    end
end