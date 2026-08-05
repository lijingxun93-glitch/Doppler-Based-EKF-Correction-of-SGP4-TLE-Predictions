%% ========================================================================
%  Two-ground-station Doppler EKF for SGP4/TLE correction
%
%  Rewritten for CSV files with columns:
%    Time_UTC, Azimuth_deg, Elevation_deg, Range_m, Doppler_Hz
%
%  State:
%    x = [dM, dRAAN, dARGP, dINC, de, dn, dBstar]^T
%
%  Measurement:
%    z = available Doppler_Hz from Tehran / Tabriz
%
%  Required external files/functions:
%    HST_20260523.tle
%    EOP-All.txt
%    SAT_Const.m, sgp4.m, teme2ecef.m, IERS.m, timediff.m,
%    days2mdh.m, Mjday.m
%% ========================================================================

clear; clc;
format long g

global const
SAT_Const

%% ----------------------------- User settings ------------------------------
% Put this script in the same folder as the two CSV files, or replace these 
% with full paths.
measFile(1) = "C:\Users\22743\Desktop\Master project\satellite data\HST 20260704\2 station Tehran SGP4_Az_El_Doppler_Visible.csv";
measFile(2) = "C:\Users\22743\Desktop\Master project\satellite data\HST 20260704\2 station Tabriz SGP4_Az_El_Doppler_Visible.csv";
stationName = ["Tehran", "Tabriz"];

% Carrier frequency used for Doppler calculation.
fc = 100e6;                    % Hz
c  = 299792.458;               % km/s

% EKF measurement noise. Unit: Hz^2.
% Increase these values if the correction is too aggressive.
sigma_D = 10;  %AWGN % Hz
R_full = diag([(sigma_D+10)^2, (sigma_D+10)^2]);

%% ----------------------------- Ground stations ----------------------------
% Station 1: Tehran
gs(1).lat_deg = 35.6833;
gs(1).lon_deg = 51.35;
gs(1).alt_m   = 1191;

% Station 2: Tabriz
gs(2).lat_deg = 38.0833;
gs(2).lon_deg = 46.2833;
gs(2).alt_m   = 1361;

nStation = numel(gs);

rsite_ecef = zeros(3,nStation);
vsite_ecef = zeros(3,nStation);

% for s = 1:nStation
%     [rsite_ecef(:,s), vsite_ecef(:,s)] = groundstation_ecef_local( ...
%         gs(s).lat_deg, gs(s).lon_deg, gs(s).alt_m/1000);
% end
for s = 1:nStation
    [rsite_ecef(:,s), ~] = groundstation_ecef_local( ...
        gs(s).lat_deg, gs(s).lon_deg, gs(s).alt_m/1000);
end

% Fixed ground stations have zero velocity in the ECEF frame.
% Do not subtract omega x r again if teme2ecef already returns ECEF velocity.
vsite_ecef(:,:) = 0;
%% ----------------------------- CSV data -----------------------------------
% This version uses absolute UTC time directly from Time_UTC.
% It does NOT use measStartTime or t_video_sec_.

measured_data = cell(nStation,1);
t_real        = cell(nStation,1);
idx_visible   = cell(nStation,1);
idx_doppler   = cell(nStation,1);

for s = 1:nStation
    opts = detectImportOptions(measFile(s));
    opts.VariableNamingRule = 'preserve';
    Tcsv = readtable(measFile(s), opts);

    measured_data{s}.Time_UTC      = parse_utc_time_local(Tcsv.("Time_UTC"));
    measured_data{s}.Azimuth_deg   = Tcsv.("Azimuth_deg");
    measured_data{s}.Elevation_deg = Tcsv.("Elevation_deg");
    measured_data{s}.Range_m       = Tcsv.("Range_m");
    %AWGN
    rng(s);       % for reproducibility
    measured_data{s}.Doppler_Hz    = Tcsv.("Doppler_Hz") + sigma_D * randn(size(Tcsv.("Doppler_Hz")));

    t_real{s} = measured_data{s}.Time_UTC;

    idx_visible{s} = find(measured_data{s}.Elevation_deg > 0);
    idx_doppler{s} = find(measured_data{s}.Elevation_deg > 0 & ...
                          ~isnan(measured_data{s}.Doppler_Hz));

    fprintf('%s CSV: %s to %s, visible samples = %d, Doppler samples = %d\n', ...
        stationName(s), string(t_real{s}(1)), string(t_real{s}(end)), ...
        numel(idx_visible{s}), numel(idx_doppler{s}));
end

% Common time interval of both station files.
startTime = max([t_real{1}(1),   t_real{2}(1)]);
endTime   = min([t_real{1}(end), t_real{2}(end)]);

if endTime <= startTime
    error('The two CSV files have no overlapping UTC time interval.');
end

%% ----------------------------- TLE initialization -------------------------
ge = 398600.8;                 % Earth gravitational constant [km^3/s^2]
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = MINUTES_PER_DAY^2;
MINUTES_PER_DAY_CUBED = MINUTES_PER_DAY^3;

% TLE file used by MatlabPredict.m
tleFile = "HST_20260523.tle";

% Read a standard .tle file. The file may contain an optional satellite-name
% line before the two element lines.
tleLines = strip(readlines(tleFile));
tleLines(tleLines == "") = [];

line1Index = find(startsWith(tleLines, "1 "), 1, "first");
line2Index = find(startsWith(tleLines, "2 "), 1, "first");

if isempty(line1Index) || isempty(line2Index) || line2Index <= line1Index
    error("Invalid TLE file '%s': TLE line 1 or line 2 is missing.", tleFile);
end

% Parse TLE line 1
tline = char(tleLines(line1Index));
Cnum = tline(3:7);
SC   = tline(8);
ID   = tline(10:17);
year = str2double(tline(19:20));
doy  = str2double(tline(21:32));
epoch = str2double(tline(19:32));
TD1   = str2double(tline(34:43));
TD2   = str2double(tline(45:50));
ExTD2 = str2double(tline(51:52));
BStar = str2double(tline(54:59));
ExBStar = str2double(tline(60:61));
BStar = BStar*1e-5*10^ExBStar;
Etype = tline(63);
Enum  = str2double(tline(65:end)); %#ok<NASGU>

% Read TLE line 2
tline = char(tleLines(line2Index));
i     = str2double(tline(9:16));
raan  = str2double(tline(18:25));
e     = str2double(strcat('0.',tline(27:33)));
omega = str2double(tline(35:42));
M     = str2double(tline(44:51));
no    = str2double(tline(53:63));
a     = ( ge/(no*2*pi/86400)^2 )^(1/3); %#ok<NASGU>
rNo   = str2double(tline(65:end));

satdata.epoch             = epoch;
satdata.norad_number      = Cnum;
satdata.bulletin_number   = ID;
satdata.classification    = SC;
satdata.revolution_number = rNo;
satdata.ephemeris_type    = Etype;
satdata.xmo               = M * pi/180;
satdata.xnodeo            = raan * pi/180;
satdata.omegao            = omega * pi/180;
satdata.xincl             = i * pi/180;
satdata.eo                = e;
satdata.xno               = no * TWOPI / MINUTES_PER_DAY;
satdata.xndt2o            = TD1 * TWOPI / MINUTES_PER_DAY_SQUARED;
satdata.xndd6o            = TD2 * 10^ExTD2 * TWOPI / MINUTES_PER_DAY_CUBED;
satdata.bstar             = BStar;

%% ----------------------------- Inject initial TLE error --------------------
% This block intentionally creates an inaccurate initial TLE estimate, so
% the EKF has something to correct. If you want to start from the original
% TLE, set satdata_err = satdata and comment out the perturbations.

satdata_err = satdata;

%% ----------------------------- EOP data -----------------------------------
fid = fopen('EOP-All.txt','r');
if fid < 0
    error('Cannot open EOP-All.txt');
end

eopdata = [];
while ~feof(fid)
    tline = fgetl(fid);
    ktmp = strfind(tline,'NUM_OBSERVED_POINTS');
    if ktmp == 1
        numrecsobs = str2double(tline(21:end));
        tline = fgetl(fid); %#ok<NASGU>
        for ii = 1:numrecsobs
            eopdata(:,ii) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]); %#ok<SAGROW>
        end
        for ii = 1:4
            tline = fgetl(fid); %#ok<NASGU>
        end
        numrecspred = str2double(tline(22:end));
        tline = fgetl(fid); %#ok<NASGU>
        for ii = numrecsobs+1:numrecsobs+numrecspred
            eopdata(:,ii) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]); %#ok<SAGROW>
        end
        break
    end
end
fclose(fid);

%% ----------------------------- Epoch time ---------------------------------
if year < 57
    year_full = year + 2000;
else
    year_full = year + 1900;
end

TimeShift = -minutes(1.116666666666);
if abs(minutes(TimeShift)) > 1
    satdata_err.epoch = satdata_err.epoch + minutes(TimeShift)/60/24;
    doy = doy + minutes(TimeShift)/60/24;
end
[mon,day,hr,minute,sec] = days2mdh(year_full,doy);
MJD_Epoch = Mjday(year_full,mon,day,hr,minute,sec);
t_epoch   = datetime(year_full,mon,day,hr,minute,sec,'TimeZone','UTC');

%% ----------------------------- Simulation time ----------------------------
dt = 1;                                      % second
duration_sec = 1800;%floor(seconds(endTime - startTime));
N = duration_sec + 1;

t_EKF = startTime + seconds(0:duration_sec);
t = seconds(t_EKF - startTime)/60;           % minutes after startTime

t_start = startTime - t_epoch;               % duration
fprintf('EKF interval: %s to %s, N = %d\n', ...
    string(startTime), string(endTime), N);

%% ----------------------------- Preallocation ------------------------------
dopplerShift = nan(N,nStation);
azimuth      = nan(N,nStation);
elevation    = nan(N,nStation);

doppler_meas_match = nan(N,nStation);
doppler_error      = nan(N,nStation);

update_flag = false(N,1);
used_station_flag = false(N,nStation);

x_history = nan(N,7);

%% ----------------------------- EKF setup ----------------------------------
% x = [dM, dRAAN, dARGP, dINC, de, dn, dBstar]^T
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

%% ============================= EKF main loop ===============================
for j = 1:N

    t_model = t_EKF(j);

    % ---------- time update ----------
    x = F * x;
    P = F * P * F' + Q;

    % ---------- UTC / EOP ----------
    tsince = t(j);                              % minutes after startTime
    tsince_total = minutes(t_start) + tsince;   % minutes after TLE epoch

    MJD_UTC = MJD_Epoch + tsince_total/1440;

    [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = ... %#ok<ASGLU>
        IERS(eopdata,MJD_UTC,'l');

    [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC); %#ok<ASGLU>

    MJD_UT1 = MJD_UTC + UT1_UTC/86400;
    MJD_TT  = MJD_UTC + TT_UTC/86400;
    T       = (MJD_TT - const.MJD_J2000)/36525;

    % ---------- prediction before update ----------
    satdata_est = apply_tle_state_correction_local(satdata_err, x);

    [rteme_est, vteme_est] = sgp4(tsince_total, satdata_est);
    [recef_est, vecef_est] = teme2ecef( ...
        rteme_est, vteme_est, T, MJD_UT1+2400000.5, ...
        LOD, x_pole, y_pole, 2);

    for s = 1:nStation
        [dopplerShift(j,s), azimuth(j,s), elevation(j,s)] = ...
            predict_one_station_local( ...
                recef_est, vecef_est, ...
                rsite_ecef(:,s), vsite_ecef(:,s), ...
                gs(s).lat_deg, gs(s).lon_deg, fc, c);
    end

    % ---------- find matched measurements ----------
    z_all     = nan(nStation,1);
    zhat_all  = dopplerShift(j,:).';
    validObs  = false(nStation,1);

    for s = 1:nStation
        [dt_min, k_meas] = min(abs(t_real{s} - t_model));
        time_ok = seconds(dt_min) <= 0.6;

        if time_ok
            doppler_meas_match(j,s) = measured_data{s}.Doppler_Hz(k_meas);

            % Use this station only if CSV elevation is above horizon and
            % Doppler_Hz is available. In your CSV, Doppler_Hz is NaN below
            % the horizon.
            if measured_data{s}.Elevation_deg(k_meas) > 0 && ...
               ~isnan(measured_data{s}.Doppler_Hz(k_meas))

                z_all(s) = measured_data{s}.Doppler_Hz(k_meas);
                doppler_error(j,s) = z_all(s) - zhat_all(s);
                validObs(s) = true;
            end
        end
    end

    % ---------- EKF update ----------
    if any(validObs)

        z     = z_all(validObs);
        z_hat = zhat_all(validObs);

        H_full = numerical_H_doppler_tle_multi_local( ...
            x, satdata_err, tsince_total, ...
            T, MJD_UT1, LOD, x_pole, y_pole, ...
            rsite_ecef, vsite_ecef, gs, fc, c);

        H = H_full(validObs,:);
        R = R_full(validObs,validObs);

        innov = z - z_hat;

        S = H * P * H' + R;
        K = (P * H') / S;

        x = x + K * innov;

        % Joseph form for numerical stability.
        I7 = eye(7);
        P = (I7 - K*H) * P * (I7 - K*H)' + K * R * K';

        update_flag(j) = true;
        used_station_flag(j,validObs) = true;

        % ---------- recompute after update ----------
        satdata_est = apply_tle_state_correction_local(satdata_err, x);

        [rteme_est, vteme_est] = sgp4(tsince_total, satdata_est);
        [recef_est, vecef_est] = teme2ecef( ...
            rteme_est, vteme_est, T, MJD_UT1+2400000.5, ...
            LOD, x_pole, y_pole, 2);

        for s = 1:nStation
            [dopplerShift(j,s), azimuth(j,s), elevation(j,s)] = ...
                predict_one_station_local( ...
                    recef_est, vecef_est, ...
                    rsite_ecef(:,s), vsite_ecef(:,s), ...
                    gs(s).lat_deg, gs(s).lon_deg, fc, c);

            if validObs(s)
                doppler_error(j,s) = doppler_meas_match(j,s) - dopplerShift(j,s);
            end
        end
    end

    x_history(j,:) = x.';
end

%% ============================= Plot results ================================
figure(6);
set(gcf,'Name','Two-station Doppler EKF','Color','w');

for s = 1:nStation
    subplot(2,2,s);

    plot(t_real{s}(idx_doppler{s}), measured_data{s}.Doppler_Hz(idx_doppler{s}), ...
        'LineWidth', 1.0);
    grid on; hold on;

    idx_model = elevation(:,s) > 0;
    plot(t_EKF(idx_model), dopplerShift(idx_model,s), 'LineWidth', 1.2);

    ylabel('Doppler (Hz)');
    title(sprintf('%s Doppler', stationName(s)));
    legend('TLE epoch at 2026-07-02','TLE epoch at 2026-05-23','Location','best');
    hold off;
end
for s = 1:nStation
    subplot(2,2,2+s);

    idx_valid = find(isfinite(doppler_error(:,s)));
    idx_rmse = idx_valid(max(1, end-29):end);
    doppler_rmse = sqrt(mean(doppler_error(idx_rmse,s).^2));
    plot(t_EKF(idx_valid), doppler_error(idx_valid,s), 'LineWidth', 1.2);
    grid on;
    yline(0,'--');

    ylabel('Error (Hz)');
    xlabel('Time UTC');
    title(sprintf('%s Doppler Error (RMSE of last 30 valid errors = %.3f Hz)', ...
        stationName(s), doppler_rmse));
    hold off;
end

format_datetime_axes_local(gcf);
set_datatips_english_local(gcf);

figure(7);
set(gcf,'Name','Az El comparison and error','Color','w');

for s = 1:nStation
    idx_model = elevation(:,s) > 0;

    t_plot = t_EKF(idx_model);
    t_plot = t_plot(:);

    az_meas_on_ekf = interp1(t_real{s}(idx_visible{s}), ...
        measured_data{s}.Azimuth_deg(idx_visible{s}), ...
        t_plot, 'linear', NaN);

    el_meas_on_ekf = interp1(t_real{s}(idx_visible{s}), ...
        measured_data{s}.Elevation_deg(idx_visible{s}), ...
        t_plot, 'linear', NaN);

    az_meas_on_ekf = az_meas_on_ekf(:);
    el_meas_on_ekf = el_meas_on_ekf(:);

    az_ekf = azimuth(idx_model,s);
    el_ekf = elevation(idx_model,s);

    az_ekf = az_ekf(:);
    el_ekf = el_ekf(:);

    az_err = wrapTo180_local(az_meas_on_ekf - az_ekf);
    el_err = el_meas_on_ekf - el_ekf;

    idx_az_valid = find(isfinite(az_err));
    idx_el_valid = find(isfinite(el_err));
    idx_az_rmse = idx_az_valid(max(1, end-29):end);
    idx_el_rmse = idx_el_valid(max(1, end-29):end);
    az_rmse = sqrt(mean(az_err(idx_az_rmse).^2));
    el_rmse = sqrt(mean(el_err(idx_el_rmse).^2));

    % -------- Azimuth comparison --------
    subplot(4,nStation,s);
    plot(t_real{s}(idx_visible{s}), measured_data{s}.Azimuth_deg(idx_visible{s}), ...
        'LineWidth', 1.0);
    hold on; grid on;
    plot(t_EKF(idx_model), azimuth(idx_model,s), 'LineWidth', 1.2);
    ylabel('Azimuth (deg)');
    title(sprintf('%s Azimuth', stationName(s)));
    legend('TLE epoch at 2026-07-02','TLE epoch at 2026-05-23','Location','best');
    hold off;

    % -------- Azimuth error --------
    subplot(4,nStation,nStation+s);
    plot(t_plot, az_err, 'LineWidth', 1.2);
    grid on; yline(0,'--');
    ylabel('Az Error (deg)');
    title(sprintf('%s Azimuth Error (RMSE of last 30 valid errors = %.3f deg)', ...
        stationName(s), az_rmse));

    % -------- Elevation comparison --------
    subplot(4,nStation,2*nStation+s);
    plot(t_real{s}(idx_visible{s}), measured_data{s}.Elevation_deg(idx_visible{s}), ...
        'LineWidth', 1.0);
    hold on; grid on;
    plot(t_EKF(idx_model), elevation(idx_model,s), 'LineWidth', 1.2);
    ylabel('Elevation (deg)');
    title(sprintf('%s Elevation', stationName(s)));
    legend('TLE epoch at 2026-07-02','TLE epoch at 2026-05-23','Location','best');
    hold off;

    % -------- Elevation error --------
    subplot(4,nStation,3*nStation+s);
    plot(t_plot, el_err, 'LineWidth', 1.2);
    grid on; yline(0,'--');
    ylabel('El Error (deg)');
    xlabel('Time UTC');
    title(sprintf('%s Elevation Error (RMSE of last 30 valid errors  = %.3f deg)', ...
        stationName(s), el_rmse));
end

format_datetime_axes_local(gcf);
set_datatips_english_local(gcf);

figure(3);
set(gcf,'Name','EKF correction states','Color','w');

subplot(2,1,1);
plot(t_EKF, rad2deg(x_history(:,1)), 'LineWidth', 1.0); hold on;
plot(t_EKF, rad2deg(x_history(:,2)), 'LineWidth', 1.0);
plot(t_EKF, rad2deg(x_history(:,3)), 'LineWidth', 1.0);
plot(t_EKF, rad2deg(x_history(:,4)), 'LineWidth', 1.0);
grid on;
ylabel('Correction (deg)');
title('Angular TLE corrections');
legend('dM','dRAAN','dARGP','dINC','Location','best');
hold off;

subplot(2,1,2);
plot(t_EKF, x_history(:,5), 'LineWidth', 1.0); hold on;
plot(t_EKF, x_history(:,6), 'LineWidth', 1.0);
plot(t_EKF, x_history(:,7), 'LineWidth', 1.0);
grid on;
ylabel('Correction');
xlabel('Time UTC');
title('Non-angular TLE corrections');
legend('de','dn rad/min','dBstar','Location','best');
hold off;

format_datetime_axes_local(gcf);
set_datatips_english_local(gcf);

fprintf('\nFinal correction x:\n');
fprintf('dM     = %.8g deg\n', rad2deg(x(1)));
fprintf('dRAAN  = %.8g deg\n', rad2deg(x(2)));
fprintf('dARGP  = %.8g deg\n', rad2deg(x(3)));
fprintf('dINC   = %.8g deg\n', rad2deg(x(4)));
fprintf('de     = %.8g\n', x(5));
fprintf('dn     = %.8g rad/min\n', x(6));
fprintf('dBstar = %.8g\n', x(7));

%% ========================================================================
%  Local functions
%% ========================================================================

function t = parse_utc_time_local(x)
% Parse UTC time column from CSV.
% Supports string/cellstr/datetime input.

    if isdatetime(x)
        t = x;
        if isempty(t.TimeZone)
            t.TimeZone = 'UTC';
        else
            t.TimeZone = 'UTC';
        end
        return;
    end

    if iscell(x)
        x = string(x);
    end

    x = string(x);
    t = datetime(x, ...
        'InputFormat','yyyy-MM-dd HH:mm:ss', ...
        'TimeZone','UTC');
end


function satdata_new = apply_tle_state_correction_local(satdata_base, x)
% Apply EKF correction to TLE/SGP4 mean elements.
%
% x = [dM, dRAAN, dARGP, dINC, de, dn, dBstar]^T

    satdata_new = satdata_base;

    satdata_new.xmo    = wrapTo2Pi_local(satdata_base.xmo    + x(1));
    satdata_new.xnodeo = wrapTo2Pi_local(satdata_base.xnodeo + x(2));
    satdata_new.omegao = wrapTo2Pi_local(satdata_base.omegao + x(3));

    satdata_new.xincl = satdata_base.xincl + x(4);
    satdata_new.xincl = max(min(satdata_new.xincl, pi), 0);

    satdata_new.eo = satdata_base.eo + x(5);
    satdata_new.eo = max(min(satdata_new.eo, 0.99), 1e-8);

    satdata_new.xno = satdata_base.xno + x(6);
    satdata_new.xno = max(satdata_new.xno, 1e-10);

    satdata_new.bstar = satdata_base.bstar + x(7);
end


function [doppler_hz, az_deg, el_deg] = predict_one_station_local( ...
        recef, vecef, rsite_ecef, vsite_ecef, lat_deg, lon_deg, fc, c)

    rho_ecef    = recef - rsite_ecef;
    rhodot_ecef = vecef - vsite_ecef;

    [enu, ~] = ecef2enu_vector_local(rho_ecef, rhodot_ecef, lat_deg, lon_deg);

    E  = enu(1);
    Nn = enu(2);
    U  = enu(3);

    range_km = norm(enu);

    el = asin(U / range_km);
    az = atan2(E, Nn);
    if az < 0
        az = az + 2*pi;
    end

    az_deg = rad2deg(az);
    el_deg = rad2deg(el);

    rho_unit = rho_ecef / norm(rho_ecef);
    range_rate = dot(rhodot_ecef, rho_unit);       % km/s

    doppler_hz = - (range_rate / c) * fc;          % Hz
end


function z_hat = meas_model_doppler_tle_multi_local( ...
        x, satdata_base, tsince_total, ...
        T, MJD_UT1, LOD, x_pole, y_pole, ...
        rsite_ecef, vsite_ecef, gs, fc, c)

    nStation = numel(gs);
    z_hat = zeros(nStation,1);

    satdata_est = apply_tle_state_correction_local(satdata_base, x);

    [rteme, vteme] = sgp4(tsince_total, satdata_est);
    [recef, vecef] = teme2ecef( ...
        rteme, vteme, T, MJD_UT1+2400000.5, ...
        LOD, x_pole, y_pole, 2);

    for s = 1:nStation
        [doppler_hz, ~, ~] = predict_one_station_local( ...
            recef, vecef, ...
            rsite_ecef(:,s), vsite_ecef(:,s), ...
            gs(s).lat_deg, gs(s).lon_deg, fc, c);

        z_hat(s) = doppler_hz;
    end
end


function H = numerical_H_doppler_tle_multi_local( ...
        x, satdata_base, tsince_total, ...
        T, MJD_UT1, LOD, x_pole, y_pole, ...
        rsite_ecef, vsite_ecef, gs, fc, c)

    h0 = meas_model_doppler_tle_multi_local( ...
        x, satdata_base, tsince_total, ...
        T, MJD_UT1, LOD, x_pole, y_pole, ...
        rsite_ecef, vsite_ecef, gs, fc, c);

    n_state = length(x);
    n_obs = length(h0);
    H = zeros(n_obs,n_state);

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

        hp = meas_model_doppler_tle_multi_local( ...
            xp, satdata_base, tsince_total, ...
            T, MJD_UT1, LOD, x_pole, y_pole, ...
            rsite_ecef, vsite_ecef, gs, fc, c);

        H(:,ii) = (hp - h0) / dx_step(ii);
    end
end


function [rsite_ecef, vsite_ecef] = groundstation_ecef_local(lat_deg, lon_deg, alt_km)
% WGS-84 ground station ECEF position and Earth-rotation velocity.
% Output:
%   rsite_ecef km
%   vsite_ecef km/s

    a = 6378.137;                       % km
    f = 1/298.257223563;
    e2 = f*(2-f);

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    N = a / sqrt(1 - e2*sin(lat)^2);

    x = (N + alt_km) * cos(lat) * cos(lon);
    y = (N + alt_km) * cos(lat) * sin(lon);
    z = (N*(1-e2) + alt_km) * sin(lat);

    rsite_ecef = [x; y; z];

    omega_earth = 7.2921150e-5;         % rad/s
    vsite_ecef = cross([0;0;omega_earth], rsite_ecef);
end


function [enu, enu_dot] = ecef2enu_vector_local(rho_ecef, rhodot_ecef, lat_deg, lon_deg)
% Convert ECEF relative vector to local ENU vector.

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    R = [ ...
        -sin(lon),              cos(lon),             0; ...
        -sin(lat)*cos(lon),    -sin(lat)*sin(lon),    cos(lat); ...
         cos(lat)*cos(lon),     cos(lat)*sin(lon),    sin(lat) ...
    ];

    enu = R * rho_ecef;
    enu_dot = R * rhodot_ecef;
end


function angle_out = wrapTo2Pi_local(angle_in)
% Wrap angle in radians to [0, 2*pi).

    angle_out = mod(angle_in, 2*pi);
    if angle_out < 0
        angle_out = angle_out + 2*pi;
    end
end


function a = wrapTo180_local(a)
% Wrap angle in degrees to [-180, 180).

    a = mod(a + 180, 360) - 180;
end


function format_datetime_axes_local(figHandle)
% Make datetime axis labels stable and English-style.

    axs = findall(figHandle, 'Type', 'axes');

    for k = 1:numel(axs)
        ax = axs(k);

        if isa(ax.XAxis, 'matlab.graphics.axis.decorator.DatetimeRuler')
            ax.XAxis.TickLabelFormat = 'HH:mm:ss';
            ax.XAxis.SecondaryLabel.String = '';
        end
    end
end


function set_datatips_english_local(figHandle)
% Force line data tips to English labels.

    lines = findall(figHandle, 'Type', 'Line');

    for k = 1:numel(lines)
        dt = lines(k).DataTipTemplate;

        if numel(dt.DataTipRows) < 2
            continue;
        end

        dt.DataTipRows(1).Label = 'Time';
        dt.DataTipRows(1).Format = 'dd-MM-yyyy HH:mm:ss';

        ax = ancestor(lines(k), 'axes');
        ylab = ax.YLabel.String;

        if contains(ylab, 'Doppler')
            dt.DataTipRows(2).Label = 'Doppler';
            dt.DataTipRows(2).Format = '%.6f Hz';

        elseif contains(ylab, 'Az Error')
            dt.DataTipRows(2).Label = 'Azimuth Error';
            dt.DataTipRows(2).Format = '%.6f deg';

        elseif contains(ylab, 'Azimuth')
            dt.DataTipRows(2).Label = 'Azimuth';
            dt.DataTipRows(2).Format = '%.6f deg';

        elseif contains(ylab, 'El Error')
            dt.DataTipRows(2).Label = 'Elevation Error';
            dt.DataTipRows(2).Format = '%.6f deg';

        elseif contains(ylab, 'Elevation')
            dt.DataTipRows(2).Label = 'Elevation';
            dt.DataTipRows(2).Format = '%.6f deg';

        elseif contains(ylab, 'Correction')
            dt.DataTipRows(2).Label = 'Correction';
            dt.DataTipRows(2).Format = '%.10g';

        else
            dt.DataTipRows(2).Label = 'Value';
            dt.DataTipRows(2).Format = '%.6f';
        end
    end
end
