measured_data = readtable('C:\Users\22743\Desktop\Master project\satellite data\20251101\az_el_doppler_1s.csv'); % Prealloc 
t_real = datetime(2025,10,31,13,19,11,'TimeZone','UTC') + seconds(measured_data.t_video_sec_); 
fr_real = [measured_data.doppler_hz;0] - [0;measured_data.doppler_hz]; 
fr_real(end) = []; fr_real(1) = 0; 
idx_real = find(measured_data.elevation_deg>0); 


figure(1);
subplot(2,2,1);
plot(t_real(idx_real),measured_data.doppler_hz(idx_real)); grid on; hold on;
ylabel('Doppler (Hz)'); title('Doppler shift'); 
subplot(2,2,2); 
plot(t_real(idx_real),movmean(fr_real(idx_real),5)); grid on; hold on;
ylabel('Doppler rate (Hz/s)'); title('Doppler rate'); 
subplot(2,2,3);
ylabel('Angles (deg)'); 
plot(t_real(idx_real),measured_data.azimuth_deg(idx_real));  grid on; hold on;
title('Az angles'); 
subplot(2,2,4);
ylabel('Angles (deg)'); 
plot(t_real(idx_real),measured_data.elevation_deg(idx_real));grid on; hold on;
title('El angles'); % xlabel('Time (UTC)'); 

%% ------------------------------------Model initialization-----------------------

format long g

global const
SAT_Const

ge = 398600.8; % Earth gravitational constant [km3/s2]
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);

% TLE file name
fname = 'COSMOS_2428.txt';

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

% satdata_err.xmo    = satdata.xmo    + deg2rad(0.01);
% satdata_err.xnodeo = satdata.xnodeo + deg2rad(0.01);
% satdata_err.omegao = satdata.omegao + deg2rad(0.01);
% satdata_err.xincl  = satdata.xincl  + deg2rad(0.005);
% satdata_err.eo     = satdata.eo     + 1e-5;
% satdata_err.xno    = satdata.xno    * (1 + 1e-6);
% satdata_err.bstar  = satdata.bstar  * (1 + 0.1);

satdata_err.xmo    = satdata.xmo    + deg2rad(0.5);   % Mean anomaly: orbital phase
satdata_err.xnodeo = satdata.xnodeo + deg2rad(0.2);   % RAAN: orbit plane orientation
satdata_err.omegao = satdata.omegao + deg2rad(0.2);   % Argument of perigee: ellipse orientation in orbital plane
satdata_err.xincl  = satdata.xincl  + deg2rad(0.1);   % Inclination: orbit tilt
satdata_err.eo     = satdata.eo     + 1e-4;           % Eccentricity: orbit shape
satdata_err.xno    = satdata.xno    * (1 + 5e-5);     % Mean motion: average angular rate
satdata_err.bstar  = satdata.bstar  * (1 + 1.0);      % B*: atmospheric drag parameter

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
gs.lat_deg  = 35.6833;
gs.lon_deg  = 51.35;
gs.alt_m    = 1191;

% ------ choose model start time
startTime = datetime(2025,10,31,13,19,11,'TimeZone','UTC');  
t_epoch   = datetime(year,mon,day,hr,minute,sec,'TimeZone','UTC');
t_start   = startTime - t_epoch;   % duration

duration_min = floor(seconds(t_real(end) - startTime)/60);
dt = 1;   % second
t = (0:duration_min*60)/60;   % in minutes for sgp4 tsince

fc = 100e6;
c = 299792.458;   % km/s

N = length(t);
%% -------------------- Compare nominal SGP4 and disturbed SGP4 --------------------

N = length(t);

doppler_nom = zeros(N,1);
az_nom      = zeros(N,1);
el_nom      = zeros(N,1);

doppler_err = zeros(N,1);
az_err      = zeros(N,1);
el_err      = zeros(N,1);

j = 1;

for tsince = t

    % current model time
    t_model = startTime + seconds((j-1)*dt);

    % ---------- UTC / EOP ----------
    MJD_UTC = MJD_Epoch + (minutes(t_start) + tsince)/1440;

    [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
    [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    MJD_UT1 = MJD_UTC + UT1_UTC/86400;
    MJD_TT  = MJD_UTC + TT_UTC/86400;
    T = (MJD_TT-const.MJD_J2000)/36525;

    % ---------- ground station ----------
    [rsite_ecef, vsite_ecef] = groundstation_ecef(gs.lat_deg, gs.lon_deg, gs.alt_m/1000);

    %% ---------- 1. Nominal SGP4 ----------
    [rteme_nom, vteme_nom] = sgp4(minutes(t_start) + tsince, satdata);
    [recef_nom, vecef_nom] = teme2ecef(rteme_nom, vteme_nom, T, MJD_UT1+2400000.5, LOD, x_pole, y_pole, 2);

    [az_nom(j), el_nom(j), doppler_nom(j)] = calc_az_el_doppler( ...
        recef_nom, vecef_nom, rsite_ecef, vsite_ecef, gs, fc, c);

    %% ---------- 2. Disturbed SGP4 ----------
    [rteme_err, vteme_err] = sgp4(minutes(t_start) + tsince, satdata_err);
    [recef_err, vecef_err] = teme2ecef(rteme_err, vteme_err, T, MJD_UT1+2400000.5, LOD, x_pole, y_pole, 2);

    [az_err(j), el_err(j), doppler_err(j)] = calc_az_el_doppler( ...
        recef_err, vecef_err, rsite_ecef, vsite_ecef, gs, fc, c);

    j = j + 1;
end

%% ---------- Doppler rate ----------
dopplerRate_nom = [doppler_nom;0] - [0;doppler_nom];
dopplerRate_nom(1) = [];
dopplerRate_nom(end) = 0;

dopplerRate_err = [doppler_err;0] - [0;doppler_err];
dopplerRate_err(1) = [];
dopplerRate_err(end) = 0;

%% ---------- Differences caused by manually injected TLE error ----------
dAz = wrapTo180(az_err - az_nom);
dEl = el_err - el_nom;
dDoppler = doppler_err - doppler_nom;
dDopplerRate = dopplerRate_err - dopplerRate_nom;

fprintf('\nInjected TLE error effect:\n');
fprintf('Az RMS error       = %.6f deg\n', rms(dAz));
fprintf('El RMS error       = %.6f deg\n', rms(dEl));
fprintf('Doppler RMS error  = %.6f Hz\n',  rms(dDoppler));
fprintf('Max |Az error|     = %.6f deg\n', max(abs(dAz)));
fprintf('Max |El error|     = %.6f deg\n', max(abs(dEl)));
fprintf('Max |Doppler err|  = %.6f Hz\n',  max(abs(dDoppler)));

%% ---------- Plot against measured data ----------
time_model = minutes(t) + startTime;

figure(1);

subplot(2,2,1);
plot(time_model, doppler_nom, 'LineWidth', 1.1);
plot(time_model, doppler_err, '--', 'LineWidth', 1.1);
legend('Doppler real','Doppler nominal SGP4','Doppler disturbed SGP4','Location','best');
ylabel('Doppler (Hz)');
title('Doppler shift');
grid on;

subplot(2,2,2);
plot(time_model, dopplerRate_nom, 'LineWidth', 1.1);
plot(time_model, dopplerRate_err, '--', 'LineWidth', 1.1);
legend('Doppler rate real','Doppler rate nominal SGP4','Doppler rate disturbed SGP4','Location','best');
ylabel('Doppler rate (Hz/s)');
title('Doppler rate');
grid on;

subplot(2,2,3);
plot(time_model, az_nom, 'LineWidth', 1.1);
plot(time_model, az_err, '--', 'LineWidth', 1.1);
legend('Az real','Az nominal SGP4','Az disturbed SGP4','Location','best');
ylabel('Azimuth (deg)');
title('Az angles');
grid on;

subplot(2,2,4);
plot(time_model, el_nom, 'LineWidth', 1.1);
plot(time_model, el_err, '--', 'LineWidth', 1.1);
legend('El real','El nominal SGP4','El disturbed SGP4','Location','best');
ylabel('Elevation (deg)');
title('El angles');
grid on;

%% ---------- Plot only injected error effect ----------
figure(2);

subplot(3,1,1);
plot(time_model, dAz);
ylabel('\Delta Az (deg)');
title('Effect of injected TLE error');
grid on;

subplot(3,1,2);
plot(time_model, dEl);
ylabel('\Delta El (deg)');
grid on;

subplot(3,1,3);
plot(time_model, dDoppler);
ylabel('\Delta Doppler (Hz)');
xlabel('Time UTC');
grid on;

function [az_deg, el_deg, doppler_hz] = calc_az_el_doppler( ...
    recef, vecef, rsite_ecef, vsite_ecef, gs, fc, c)

    rho_ecef    = recef - rsite_ecef;
    rhodot_ecef = vecef - vsite_ecef;

    [enu, ~] = ecef2enu_vector(rho_ecef, rhodot_ecef, gs.lat_deg, gs.lon_deg);

    E  = enu(1);
    Nn = enu(2);
    U  = enu(3);

    range_km = norm(enu);

    el = asin(U / range_km);
    az = atan2(E, Nn);

    if az < 0
        az = az + 2*pi;
    end

    rho_unit = rho_ecef / norm(rho_ecef);
    range_rate = dot(rhodot_ecef, rho_unit);   % km/s

    doppler_hz = - (range_rate / c) * fc;
    az_deg = rad2deg(az);
    el_deg = rad2deg(el);
end