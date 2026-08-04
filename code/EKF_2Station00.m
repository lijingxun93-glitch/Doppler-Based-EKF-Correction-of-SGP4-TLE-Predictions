%% ========================================================================
%  Two-ground-station Doppler EKF for SGP4/TLE correction
%  State:
%  x = [dM, dRAAN, dARGP, dINC, de, dn, dBstar]^T
%
%  Measurement:
%  z = [Doppler_station_1;
%       Doppler_station_2]
%
%  If one station is below horizon or has no matched measurement, only the
%  available station is used in EKF update.
%% ========================================================================

clear; clc; close all;
format long g

global const
SAT_Const
%% ----------------------------- Ground stations ----------------------------
% Station 1: Tehran
gs(1).lat_deg = 35.6833;
gs(1).lon_deg = 51.35;
gs(1).alt_m   = 1191;

% Station 2: Tabriz
gs(2).lat_deg = 38.0833;
gs(2).lon_deg = 46.2833;
gs(2).alt_m   = 1361;
% gs(2).lat_deg = 35.6833;
% gs(2).lon_deg = 51.35;
% gs(2).alt_m   = 1191;
nStation = 2;

rsite_ecef = zeros(3,nStation);
vsite_ecef = zeros(3,nStation);

for s = 1:nStation
    [rsite_ecef(:,s), vsite_ecef(:,s)] = groundstation_ecef_local( ...
        gs(s).lat_deg, gs(s).lon_deg, gs(s).alt_m/1000);
end
%% ----------------------------- Measured data -----------------------------
% % 注意：
% % 两个 CSV 最好都有这些列：
% %   t_video_sec_
% %   azimuth_deg
% %   elevation_deg
% %   doppler_hz
% %
% % 如果第二个站现在还没有实测数据，可以先用同一个文件测试代码结构，
% % 但真实双站效果必须用两个不同观测站的数据。
% 
% measFile(1) = "C:\Users\22743\Desktop\Master project\satellite data\20251101\az_el_doppler_2stations_v3.csv";
% measFile(2) = "C:\Users\22743\Desktop\Master project\satellite data\20251101\az_el_doppler_2stations_v3_2.csv";
% 
% measured_data = cell(2,1);
% t_real        = cell(2,1);
% idx_real      = cell(2,1);
% 
% % 两个录像/数据的 UTC 起始时间
% % 如果两个站录像起始时间不同，分别修改这里
% measStartTime(1) = datetime(2026,5,15,14,5,24,'TimeZone','UTC');
% % measStartTime(2) = datetime(2026,5,15,14,5,24,'TimeZone','UTC');
% measStartTime(2) = measStartTime(1);
% for s = 1:2
%     measured_data{s} = readtable(measFile(s));
%     t_real{s} = measStartTime(s) + seconds(measured_data{s}.t_video_sec_);
%     idx_real{s} = find(measured_data{s}.elevation_deg > 0);
% end

% persudo measured data
measured_data = cell(2,1);
t_real        = cell(2,1);
idx_real      = cell(2,1);
measStartTime(1) = datetime(2025,10,31,13,19,11,'TimeZone','UTC','Locale','en_US');
% measStartTime(2) = datetime(2026,5,15,14,5,24,'TimeZone','UTC','Locale','en_US');
measStartTime(2) = measStartTime(1);
duration_min = 16;
for s = 1:2
    [measured_data{s}.doppler_hz, measured_data{s}.azimuth_deg, measured_data{s}.elevation_deg, satData] = Doppler_graph(measStartTime(s), duration_min, gs(s).lat_deg, gs(s).lon_deg, gs(s).alt_m);
    t_real{s} = measStartTime(s) + seconds(0:duration_min*60);
    idx_real{s} = find(measured_data{s}.elevation_deg > 0);
end
%% ----------------------------- TLE initialization -------------------------
ge = 398600.8;                 % Earth gravitational constant [km^3/s^2]
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = MINUTES_PER_DAY^2;
MINUTES_PER_DAY_CUBED = MINUTES_PER_DAY^3;

fname = 'COSMOS_2428.txt';

fid = fopen(fname, 'r');
if fid < 0
    error('Cannot open TLE file: %s', fname);
end

% read first line
tline = fgetl(fid);
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
Enum  = str2double(tline(65:end));

% read second line
tline = fgetl(fid);
i     = str2double(tline(9:16));
raan  = str2double(tline(18:25));
e     = str2double(strcat('0.',tline(27:33)));
omega = str2double(tline(35:42));
M     = str2double(tline(44:51));
no    = str2double(tline(53:63));
a     = ( ge/(no*2*pi/86400)^2 )^(1/3);
rNo   = str2double(tline(65:end));

fclose(fid);

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

%% ----------------------------- Inject TLE error ---------------------------
satdata_err = satdata;

% satdata_err.xmo    = satdata.xmo    + deg2rad(0.01);
% satdata_err.xnodeo = satdata.xnodeo + deg2rad(0.01);
% satdata_err.omegao = satdata.omegao + deg2rad(0.01);
% satdata_err.xincl  = satdata.xincl  + deg2rad(0.005);
% satdata_err.eo     = satdata.eo     + 1e-5;
% satdata_err.xno    = satdata.xno    * (1 + 1e-6);
% satdata_err.bstar  = satdata.bstar  * (1 + 0.1);
satdata_err.xmo    = satdata.xmo    + deg2rad(5);   % Mean anomaly: orbital phase
satdata_err.xnodeo = satdata.xnodeo + deg2rad(2);   % RAAN: orbit plane orientation
satdata_err.omegao = satdata.omegao + deg2rad(2);   % Argument of perigee: ellipse orientation in orbital plane
satdata_err.xincl  = satdata.xincl  + deg2rad(0.5);   % Inclination: orbit tilt
satdata_err.eo     = satdata.eo     + 1e-3;           % Eccentricity: orbit shape
satdata_err.xno    = satdata.xno    * (1 + 5e-5);     % Mean motion: average angular rate
satdata_err.bstar  = satdata.bstar  * (1 + 1.0);      % B*: atmospheric drag parameter

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
        tline = fgetl(fid);
        for ii = 1:numrecsobs
            eopdata(:,ii) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        for ii = 1:4
            tline = fgetl(fid);
        end
        numrecspred = str2double(tline(22:end));
        tline = fgetl(fid);
        for ii = numrecsobs+1:numrecsobs+numrecspred
            eopdata(:,ii) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
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

[mon,day,hr,minute,sec] = days2mdh(year_full,doy);
MJD_Epoch = Mjday(year_full,mon,day,hr,minute,sec);



%% ----------------------------- Simulation time ----------------------------
startTime = measStartTime(1);
t_epoch   = datetime(year_full,mon,day,hr,minute,sec,'TimeZone','UTC');
t_start   = startTime - t_epoch;

% 用两个数据文件中最短的结束时间，避免索引越界
endTime = min(t_real{1}(end), t_real{2}(end));

duration_min = floor(seconds(endTime - startTime)/60);
dt = 1;                         % second
t = (0:duration_min*60)/60;      % minutes after startTime

fc = 100e6;                     % Hz
c  = 299792.458;                % km/s

N = length(t);
t_EKF = startTime + seconds((0:N-1)*dt);

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
    deg2rad(1.0), ...
    deg2rad(0.5), ...
    deg2rad(0.5), ...
    deg2rad(0.2), ...
    5e-4, ...
    satdata.xno*1e-4, ...
    max(abs(satdata.bstar),1e-8) ...
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

% 两个站的 Doppler measurement noise
% 如果第二站接收机更差，可以设成 [1^2, 5^2]
R_full = diag([20^2, 20^2]);      % Hz^2

%% ============================= EKF main loop ===============================
for j = 1:N

    t_model = t_EKF(j);

    % ---------- time update ----------
    x = F * x;
    P = F * P * F' + Q;

    % ---------- UTC / EOP ----------
    tsince = t(j);
    tsince_total = minutes(t_start) + tsince;

    MJD_UTC = MJD_Epoch + tsince_total/1440;

    [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = ...
        IERS(eopdata,MJD_UTC,'l');

    [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);

    MJD_UT1 = MJD_UTC + UT1_UTC/86400;
    MJD_TT  = MJD_UTC + TT_UTC/86400;
    T = (MJD_TT - const.MJD_J2000)/36525;

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

        % 允许最大时间差。你的数据是 1 s，设 0.6 s 比较严格。
        time_ok = seconds(dt_min) <= 0.6;

        if time_ok
            doppler_meas_match(j,s) = measured_data{s}.doppler_hz(k_meas);

            % 用实测 elevation 判断该站是否有效
            if measured_data{s}.elevation_deg(k_meas) > 0
                z_all(s) = measured_data{s}.doppler_hz(k_meas);
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

        % Joseph form: 比 P=(I-KH)P 数值稳定
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
figure(1);
set(gcf,'Name','Two-station Doppler EKF','Color','w');

for s = 1:nStation
    % subplot(3,2,s);
    subplot(2,2,s);
    plot(t_real{s}(idx_real{s}), measured_data{s}.doppler_hz(idx_real{s}), 'LineWidth', 1.0);
    grid on; hold on;
    idx_model = elevation(:,s) > 0;
    plot(t_EKF(idx_model), dopplerShift(idx_model,s), 'LineWidth', 1.2);
    ylabel('Doppler (Hz)');
    title(sprintf('Station %d Doppler',s));
    legend('Measured','SGP4 + EKF','Location','best');
end

for s = 1:nStation
    % subplot(3,2,2+s);
    subplot(2,2,2+s);
    idx_valid = ~isnan(doppler_error(:,s));
    plot(t_EKF(idx_valid), doppler_error(idx_valid,s), 'LineWidth', 1.2);
    grid on;
    yline(0,'--');
    ylabel('Error (Hz)');
    title(sprintf('Station %d Doppler Error',s));
end

% subplot(3,2,5);
% plot(t_EKF, used_station_flag(:,1), 'LineWidth', 1.2); hold on;
% plot(t_EKF, used_station_flag(:,2), 'LineWidth', 1.2);
% grid on;
% ylim([-0.1 1.1]);
% ylabel('Used in EKF');
% xlabel('Time UTC');
% title('Station usage');
% legend('Station 1','Station 2','Location','best');
% 
% subplot(3,2,6);
% plot(t_EKF, x_history(:,1)*180/pi, 'LineWidth', 1.0); hold on;
% plot(t_EKF, x_history(:,2)*180/pi, 'LineWidth', 1.0);
% plot(t_EKF, x_history(:,3)*180/pi, 'LineWidth', 1.0);
% plot(t_EKF, x_history(:,4)*180/pi, 'LineWidth', 1.0);
% grid on;
% ylabel('Correction (deg)');
% xlabel('Time UTC');
% title('Angular state corrections');
% legend('dM','dRAAN','dARGP','dINC','Location','best');

axs = findall(gcf, 'Type', 'axes');
for k = 1:numel(axs)
    ax = axs(k);

    if isa(ax.XAxis, 'matlab.graphics.axis.decorator.DatetimeRuler')
        ax.XAxis.TickLabelFormat = 'HH:mm';
        ax.XAxis.SecondaryLabel.String = 'May 15, 2026';
    end
end
figure(2);
set(gcf,'Name','Az El comparison and error','Color','w');

% Figure 2 layout:
%   Row 1: azimuth comparison + azimuth error
%   Row 2: elevation comparison + elevation error
% Each station uses one column block. For two stations, this gives 4 rows x 2 columns.

for s = 1:nStation
    idx_model = elevation(:,s) > 0;

    % Interpolate measured az/el onto EKF time base for error calculation.
    % This is safer than subtracting directly because measured time and EKF time
    % may not always have exactly the same samples.
    t_plot = t_EKF(idx_model);
    t_plot = t_plot(:);
    
    az_meas_on_ekf = interp1(t_real{s}(idx_real{s}), ...
        measured_data{s}.azimuth_deg(idx_real{s}), ...
        t_plot, 'linear', NaN);
    
    el_meas_on_ekf = interp1(t_real{s}(idx_real{s}), ...
        measured_data{s}.elevation_deg(idx_real{s}), ...
        t_plot, 'linear', NaN);
    
    % Force all vectors to be column vectors
    az_meas_on_ekf = az_meas_on_ekf(:);
    el_meas_on_ekf = el_meas_on_ekf(:);
    
    az_ekf = azimuth(idx_model,s);
    el_ekf = elevation(idx_model,s);
    
    az_ekf = az_ekf(:);
    el_ekf = el_ekf(:);
    
    az_err = wrapTo180_local(az_meas_on_ekf - az_ekf);
    el_err = el_meas_on_ekf - el_ekf;

    % -------- Azimuth comparison --------
    subplot(4,nStation,s);
    plot(t_real{s}(idx_real{s}), measured_data{s}.azimuth_deg(idx_real{s}), 'LineWidth', 1.0);
    hold on; grid on;
    plot(t_EKF(idx_model), azimuth(idx_model,s), 'LineWidth', 1.2);
    ylabel('Azimuth (deg)');
    title(sprintf('Station %d Azimuth',s));
    legend('Measured','SGP4 + EKF','Location','best');

    % -------- Azimuth error --------
    subplot(4,nStation,nStation+s);
    plot(t_plot, az_err, 'LineWidth', 1.2);
    grid on; yline(0,'--');
    ylabel('Az Error (deg)');
    title(sprintf('Station %d Azimuth Error',s));

    % -------- Elevation comparison --------
    subplot(4,nStation,2*nStation+s);
    plot(t_real{s}(idx_real{s}), measured_data{s}.elevation_deg(idx_real{s}), 'LineWidth', 1.0);
    hold on; grid on;
    plot(t_EKF(idx_model), elevation(idx_model,s), 'LineWidth', 1.2);
    ylabel('Elevation (deg)');
    title(sprintf('Station %d Elevation',s));
    legend('Measured','SGP4 + EKF','Location','best');

    % -------- Elevation error --------
    subplot(4,nStation,3*nStation+s);
    plot(t_plot, el_err, 'LineWidth', 1.2);
    grid on; yline(0,'--');
    ylabel('El Error (deg)');
    xlabel('Time UTC');
    title(sprintf('Station %d Elevation Error',s));
end

axs = findall(gcf, 'Type', 'axes');
for k = 1:numel(axs)
    ax = axs(k);

    if isa(ax.XAxis, 'matlab.graphics.axis.decorator.DatetimeRuler')
        ax.XAxis.TickLabelFormat = 'HH:mm';
        ax.XAxis.SecondaryLabel.String = 'May 15, 2026';
    end
end
% -------- Force data tips to English --------
lines = findall(gcf, 'Type', 'Line');

for k = 1:numel(lines)
    dt = lines(k).DataTipTemplate;

    % X value: datetime shown in English-like numeric format
    dt.DataTipRows(1).Label = 'Time';
    dt.DataTipRows(1).Format = 'dd-MM-yyyy HH:mm:ss';

    % Y value label depends on subplot ylabel
    ax = ancestor(lines(k), 'axes');
    ylab = ax.YLabel.String;

    if contains(ylab, 'Azimuth')
        dt.DataTipRows(2).Label = 'Azimuth';
        dt.DataTipRows(2).Format = '%.6f deg';

    elseif contains(ylab, 'Az Error')
        dt.DataTipRows(2).Label = 'Azimuth Error';
        dt.DataTipRows(2).Format = '%.6f deg';

    elseif contains(ylab, 'Elevation')
        dt.DataTipRows(2).Label = 'Elevation';
        dt.DataTipRows(2).Format = '%.6f deg';

    elseif contains(ylab, 'El Error')
        dt.DataTipRows(2).Label = 'Elevation Error';
        dt.DataTipRows(2).Format = '%.6f deg';

    else
        dt.DataTipRows(2).Label = 'Value';
        dt.DataTipRows(2).Format = '%.6f';
    end
end

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

function angle_out = wrapTo2Pi_local(angle_in) %[0, 2*pi]
    angle_out = mod(angle_in, 2*pi);
    if angle_out < 0
        angle_out = angle_out + 2*pi;
    end
end

function a = wrapTo180_local(a)
    % Wrap angle in degrees to [-180, 180].
    a = mod(a + 180, 360) - 180;
end