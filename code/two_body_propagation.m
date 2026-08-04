% DOPPLER_FROM_TLE_DEMO_FIXED 
% % Two-body propagation (Kepler) with key fixes: 
% - PQW state uses consistent E-parameterization efor r and v 
% - WGS-84 ellipsoid for GS ECEF 
% - Optional GAST (GMST + EqEq) rotation (simple truncated model) 
% % NOTE: This is still NOT SGP4. Expect residual bias vs Gpredict. 
%% -------- USER INPUTS -------- 
L1 = '1 25338U 98030A   25304.26992000  .00000224  00000+0  10966-3 0  9995';
L2 = '2 25338  98.5267 325.8917 0009236 301.7711  58.2570 14.27044591428643';
gs.lat_deg = 35.6833; % deg (+N) 
gs.lon_deg = 51.35; % deg (+E) 
gs.alt_m = 1191; % m 
fc = 100e6; % Hz 
duration_min= 20; 
dt_sec = 1; 
useGAST_approx = true; % <- 可切换：true 用 GMST+EqEq 近似；false 用 GMST 直接 
%% -------- CONSTANTS -------- 
mu = 3.986004418e14; % m^3/s^2 
Re = 6378137.0; % WGS-84 a 
c = 299792458; % m/s 
wE = 7.2921150e-5; % rad/s earth 
%% -------- PARSE TLE -------- 
tle = parseTLE(L1,L2); 
n0 = tle.n_revday * 2*pi/86400; % rad/s (mean motion from TLE) 
a = (mu / n0^2)^(1/3); % semi-major (二体等效)Kepler’s third law 
t_epoch = tle.epoch; % datetime UTC 
i = deg2rad(tle.inc_deg); % inclination 
raan = deg2rad(tle.raan_deg); % Right Ascension of Ascending Node 升交点赤经 
e = tle.ecc; % eccentricity 偏心率 
argp = deg2rad(tle.argp_deg); % argument of perigee 
M0 = deg2rad(tle.M_deg); % mean anomaly at epoch 
%% -------- TIME VECTOR -------- 
startTime = datetime(2025,10,31,13,11,11,'TimeZone','UTC'); 
% startTime = t_epoch; 
tvec = startTime : seconds(dt_sec) : (startTime + minutes(duration_min)); 
Nt = numel(tvec); 
measured_data = readtable('C:\Users\22743\Desktop\Master project\satellite data\20251101\az_el_doppler_1s.csv'); % Prealloc 
r_teme = zeros(3,Nt); 
v_teme = zeros(3,Nt); 
r_ecef = zeros(3,Nt); 
v_ecef = zeros(3,Nt); 
range_m = zeros(1,Nt); 
rrate = zeros(1,Nt); 
fd = zeros(1,Nt); 
az_deg = zeros(1,Nt); 
el_deg = zeros(1,Nt); % GS ECEF (WGS-84) 
gs_ecef = llh2ecef_wgs84([deg2rad(gs.lat_deg); deg2rad(gs.lon_deg); gs.alt_m]);
%% -------- PROPAGATE & DOPPLER -------- 
R3_raan = rot3(raan); 
R1_i = rot1(i); 
R3_argp = rot3(argp); 
Q_pqw2teme = R3_raan * R1_i * R3_argp; % treat as TEME-like inertial 
%n_mean = sqrt(mu/a^3); % 二体角速度（与 n0 数值相近但来源不同(存疑)） 
for k = 1:Nt 
    t = tvec(k); 
    dt = seconds(t - t_epoch); 
    % ---- Mean anomaly & Kepler E ---- 
    M = M0 + n0*dt; % 用 TLE n0 匀速推进 
    E = keplerE(M, e); 
    % ---- Consistent E-parameterization for r,v in PQW ---- Kepler’s ellipse equation r = a(1−ecosE) 
    r_peri = [ a*(cos(E)-e); % rx​ = a(cosE−e) 
        a*sqrt(1-e^2)*sin(E); % ry​ = a*sqrt(1−e2)​sinE 
        0 ]; % rz = 0 
    v_peri = [ -a*n0*sin(E); 
        +a*n0*sqrt(1-e^2)*cos(E); 
        0 ]; 
     % ---- PQW -> TEME-like inertial ---- 
     r_teme(:,k) = Q_pqw2teme * r_peri; 
     v_teme(:,k) = Q_pqw2teme * v_peri; 
     % ---- Inertial -> ECEF (GMST or GAST approx) ---- % 
     th = wE*dt;%???????????????????????？？？？？？？？？？？？ 
     % th_gmst = gmst_angle(t); 
     % if useGAST_approx 
     %    eqe = equation_of_equinoxes_trunc(t); % radians (Simplified approximation) 
     %    th = th_gmst + eqe; 
     % else 
     %     th = th_gmst; 
     % end 
     R3_g = rot3(th); 
     r_ecef(:,k) = R3_g * r_teme(:,k); 
     v_ecef(:,k) = R3_g * (v_teme(:,k) - cross([0;0;wE], r_teme(:,k))); 
     % ---- Topocentric & Doppler ---- 
     rel = r_ecef(:,k) - gs_ecef; 
     range_m(k) = norm(rel); 
     vrel = v_ecef(:,k); % GS in ECEF is stationary 
     rrate(k) = dot(vrel, rel)/range_m(k); %rangerate 
     fd(k) = -(fc/c)*rrate(k); 
     [az_deg(k), el_deg(k)] = ecef2azel(rel, gs); 
end 
fr = [fd,0] - [0,fd]; 
fr(end) = []; fr(1) = 0; 
t_real = datetime(2025,10,31,13,19,11,'TimeZone','UTC') + seconds(measured_data.t_video_sec_); 
tc = seconds(t_real(1) - startTime); 
fr_real = [measured_data.doppler_hz;0] - [0;measured_data.doppler_hz]; 
fr_real(end) = []; fr_real(1) = 0; 
error = 1; 
Delta_t = 0; 
%% -------- PLOTS -------- 
%——————SGP4 
[dopplerShift, dopplerRate, azimuth, elevation] = Doppler_graph(startTime - t_epoch, duration_min); 
idx_sgp4 = find(elevation > 0); 
figure(1); 
subplot(2,2,1); 
plot(tvec(idx_sgp4),dopplerShift(idx_sgp4));hold on; 
subplot(2,2,2);  plot(tvec(idx_sgp4),dopplerRate(idx_sgp4));hold on; 
subplot(2,2,3);  plot(tvec(idx_sgp4),azimuth(idx_sgp4));hold on; 
subplot(2,2,4);  plot(tvec(idx_sgp4),elevation(idx_sgp4));hold on; 
% figure(1); 
% subplot(2,2,1); 
% plot(tvec,dopplerShift);hold on; 
% subplot(2,2,2); plot(tvec,dopplerRate);hold on; 
% subplot(2,2,3); plot(tvec,azimuth);hold on; 
% subplot(2,2,4); plot(tvec,elevation);hold on; 
%___________________________ 
idx = find(el_deg > 0); 
idx_real = find(measured_data.elevation_deg>0); 
subplot(2,2,1); %plot(tvec(idx), fd(idx)); 
plot(t_real(idx_real),measured_data.doppler_hz(idx_real)); grid on; 
hold off;
ylabel('Doppler (Hz)'); title('Doppler angles'); 
legend('Doppler sgp4','Doppler Gpredict','Doppler real','Location','best'); % 
subplot(2,2,2); %plot(tvec(idx), fr(idx)); 
plot(t_real(idx_real),movmean(fr_real(idx_real),5)); grid on; 
hold off;
ylabel('Doppler rate (Hz/s)'); title('Doppler rate'); 
legend( 'Doppler rate sgp4','Doppler rate Gpredict','Doppler rate real','Location','best'); 
subplot(2,2,3); %plot(tvec(idx), az_deg(idx));
ylabel('Angles (deg)'); 
plot(t_real(idx_real),measured_data.azimuth_deg(idx_real));  grid on; 
hold off;
title('Az angles'); 
legend('Az sgp4','Az Gpredict','Az real', 'Location','best'); % 
subplot(2,2,4); %plot(tvec(idx), el_deg(idx)); 
ylabel('Angles (deg)'); 
plot(t_real(idx_real),measured_data.elevation_deg(idx_real));grid on; hold off;
title('El angles'); % xlabel('Time (UTC)'); 
legend('El sgp4','El Gpredict','El real', 'Location','best'); 
% subplot(2,2,1); %plot(tvec, fd); grid on; hold off;
% plot(t_real(idx_real),measured_data.doppler_hz(idx_real)); grid on; hold off;
% ylabel('Doppler (Hz)'); title('Doppler angles'); 
% legend('Doppler sgp4','Doppler Gpredict','Location','best'); 
% subplot(2,2,2); %plot(tvec, fr); grid on; hold off;
% plot(t_real(idx_real),movmean(fr_real(idx_real),5)); grid on; hold off;
% ylabel('Doppler rate (Hz/s)'); title('Doppler rate'); 
% legend( 'Doppler rate sgp4','Doppler rate Gpredict','Location','best'); 
% subplot(2,2,3); %plot(tvec, az_deg); grid on; 
% ylabel('Angles (deg)'); hold off;
% plot(t_real(idx_real),measured_data.azimuth_deg(idx_real)); grid on; hold off;
% title('Az angles'); 
% legend('Az sgp4','Az Gpredict', 'Location','best'); 
% subplot(2,2,4); %plot(tvec, el_deg); grid on; 
% ylabel('Angles (deg)'); hold off;
% plot(t_real(idx_real),measured_data.elevation_deg(idx_real)); grid on; hold off;
% title('El angles'); xlabel('Time (UTC)'); 
% legend('El sgp4','El Gpredict', 'Location','best'); 

% Visibility quick check 
% vis = el_deg > 0; 
% if any(vis) 
% t_on = tvec(find(vis,1,'first')); 
% t_off = tvec(find(vis,1,'last')); 
% fprintf('Visible %s to %s (UTC). Max El=%.1f deg\n',... % datestr(t_on), datestr(t_off), max(el_deg)); 
% else % fprintf('No visibility (El<=0).\n'); 
% end 
%% ========== helpers ========== 
function tle = parseTLE(L1,L2) 
    epoch_str = strtrim(L1(19:32)); 
    yy = str2double(epoch_str(1:2)); 
    doy = str2double(epoch_str(3:end)); 
    year = (yy<57) * (2000+yy) + (yy>=57) * (1900+yy); 
    doy_floor = floor(doy); 
    frac = doy - doy_floor; 
    epoch_dt = datetime(year,1,1,0,0,0,'TimeZone','UTC') + days(doy_floor-1) + seconds(frac*86400); 
    inc_deg = str2double(L2(9:16)); 
    raan_deg = str2double(L2(18:25)); 
    ecc = str2double(['0.' L2(27:33)]); 
    argp_deg = str2double(L2(35:42)); 
    M_deg = str2double(L2(44:51)); 
    n_revday = str2double(L2(53:63)); 
    tle = struct('epoch',epoch_dt,'inc_deg',inc_deg,'raan_deg',raan_deg,... 
        'ecc',ecc,'argp_deg',argp_deg,'M_deg',M_deg,... 
        'n_revday',n_revday); 
end 
function E = keplerE(M,e) % Kepler’s Equation M=E−esin(E) 
    E = M; 
    for it=1:25 
        f = E - e*sin(E) - M; 
        fp= 1 - e*cos(E); 
        dE= -f/fp; 
        E = E + dE; 
        if abs(dE) < 1e-13, 
            break; 
        end 
    end 
    E = atan2(sin(E), cos(E)); 
end 
function R = rot1(a), R=[1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)]; end 
function R = rot3(a), R=[cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]; end 
% Earth’s rotation relative to the fixed stars % tells how much the Earth has rotated since 0h UT 
function th = gmst_angle(tUTC) % GMST (rad), Vallado simplified 
    [JD, Tu] = julianDateAndTu(tUTC); 
    GMST_sec = 67310.54841 + (876600*3600 + 8640184.812866)*Tu + 0.093104*Tu^2 - 6.2e-6*Tu^3; 
    GMST_sec = mod(GMST_sec,86400); 
    if GMST_sec<0, 
        GMST_sec=GMST_sec+86400; 
    end 
    th = 2*pi*(GMST_sec/86400); 
end 
function [JD, Tu] = julianDateAndTu(tUTC)
    tUTC.TimeZone = 'UTC'; 
    y=year(tUTC); 
    m=month(tUTC); 
    d=day(tUTC); 
    h=hour(tUTC); 
    mi=minute(tUTC); 
    s=second(tUTC); 
    if m<=2, 
        y=y-1; 
        m=m+12; 
    end 
    A=floor(y/100); 
    B=2-A+floor(A/4); 
    JD0=floor(365.25*(y+4716))+floor(30.6001*(m+1))+d+B-1524.5; 
    frac=(h+mi/60+s/3600)/24; 
    JD=JD0+frac; 
    Tu=(JD-2451545.0)/36525; 
end 
function [az_deg, el_deg] = ecef2azel(rel_ecef, gs) 
    lat=deg2rad(gs.lat_deg); 
    lon=deg2rad(gs.lon_deg); 
    R = [-sin(lon) cos(lon) 0; 
        -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); 
        cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)]; 
    enu = R*rel_ecef; 
    E=enu(1); 
    N=enu(2); 
    U=enu(3); 
    az=atan2(E,N); 
    if az<0, 
        az=az+2*pi; 
    end 
    el=atan2(U, sqrt(E^2+N^2)); 
    az_deg=rad2deg(az); 
    el_deg=rad2deg(el); 
end 
function eqe = equation_of_equinoxes_trunc(tUTC) % Very truncated EqEq (IAU-1980–style) for demo; returns radians. 
    [JD, Tu] = julianDateAndTu(tUTC); 
    T = (JD - 2451545.0)/36525; % centuries TT ~ UT1 here (demo) % Mean longitudes (rad) 
    L = deg2rad(280.4665 + 36000.7698*T); % Sun mean long. 
    l = deg2rad(218.3165 + 481267.8813*T); % Moon mean long. (≈ argument) 
    Om = deg2rad(125.04452 - 1934.136261*T); % Ascending node of Moon % Nutation in longitude (arcsec) very truncated: 
    dPsi_as = -17.20*sin(Om) - 1.32*sin(2*L) - 0.23*sin(2*l) + 0.21*sin(2*Om); % Mean obliquity (arcsec) 
    eps0_as = 84381.448 - 46.8150*T - 0.00059*T^2 + 0.001813*T^3; eqe = deg2rad( (dPsi_as*cosd(eps0_as))/3600 ); % radians 
end
