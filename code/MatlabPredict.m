clc; clear;
%% ============================================================
%  MATLAB SGP4 + Ground Station + Azimuth/Elevation/Doppler
%  Replacement for Gpredict predicted data
%  Required toolbox: Satellite Communications Toolbox
%% ============================================================

%% ---------- 1. TLE ----------
% tleLine1 = "1 25338U 98030A   26134.95357256  .00000056  00000-0  40526-4 0  9991";
% tleLine2 = "2 25338  98.5085 156.4816 0011503  85.8931 274.3565 14.27135756456551";

tleFile = "ISS_20260727.tle";
% fid = fopen(tleFile, "w");
% fprintf(fid, "NOAA 15\n");
% fprintf(fid, "%s\n", tleLine1);
% fprintf(fid, "%s\n", tleLine2);
% fclose(fid);

%% ---------- 2. Simulation time ----------
% 注意：这里时间必须覆盖你想看的 pass
startTime  = datetime(2026,7,4,7,19,53,"TimeZone","UTC");
stopTime   = startTime + minutes(1440);
sampleTime = 1;     % seconds

sc = satelliteScenario(startTime, stopTime, sampleTime);

%% ---------- 3. Satellite using SGP4 ----------
sat = satellite(sc, tleFile, ...
    "Name", "NOAA 15", ...
    "OrbitPropagator", "sgp4");

%% ---------- 4. Ground station ----------
% Tehran
gsLat = 35.6833;      % deg
gsLon = 51.35;      % deg
gsAlt = 1191;            % metre, altitude above WGS84 ellipsoid
% gsLat = 38.0833;
% gsLon = 46.2833;
% gsAlt = 1361;
gs = groundStation(sc, ...
    "Name", "Ground Station 1", ...
    "Latitude", gsLat, ...
    "Longitude", gsLon, ...
    "Altitude", gsAlt);

%% ---------- 5. Azimuth, elevation, range ----------
[az, el, range, t] = aer(gs, sat);

az = az(:);
el = el(:);
range = range(:);
t = t(:);

%% ---------- 6. Doppler shift ----------
fc = 100e6;            % Hz, carrier frequency
c  = physconst("LightSpeed");

% Method A: MATLAB built-in dopplershift
try
    [fd_builtin, tDop] = dopplershift(sat, gs, "Frequency", fc);
    fd = fd_builtin(:);
catch
    warning("dopplershift function not available. Using manual Doppler calculation.");

    % Method B: manual Doppler from position and velocity
    [rSat, vSat] = states(sat, "CoordinateFrame", "ecef");
    [rGs,  vGs]  = states(gs,  "CoordinateFrame", "ecef");

    rSat = squeeze(rSat);   % 3 x N, m
    vSat = squeeze(vSat);   % 3 x N, m/s
    rGs  = squeeze(rGs);    % 3 x N, m
    vGs  = squeeze(vGs);    % 3 x N, m/s

    rho    = rSat - rGs;
    rhoDot = vSat - vGs;

    rangeManual = vecnorm(rho, 2, 1);
    rangeRate = sum(rho .* rhoDot, 1) ./ rangeManual;

    fd = -(fc/c) .* rangeRate(:);
end

%% ---------- 7. Visible pass only ----------
visibleIdx = el > 0;

t_visible     = t(visibleIdx);
az_visible    = az(visibleIdx);
el_visible    = el(visibleIdx);
range_visible = range(visibleIdx);
fd_visible    = fd(visibleIdx);

% NaN-masked data for plotting. Keep the complete time vector so that
% MATLAB does not connect separate visible passes with a straight line.
az_visible_plot = az;
el_visible_plot = el;
fd_visible_plot = fd;

az_visible_plot(~visibleIdx) = NaN;
el_visible_plot(~visibleIdx) = NaN;
fd_visible_plot(~visibleIdx) = NaN;

%% ---------- 8. Save all data ----------
T_all = table( ...
    t, az, el, range, fd, ...
    'VariableNames', {'Time_UTC','Azimuth_deg','Elevation_deg','Range_m','Doppler_Hz'} ...
);

writetable(T_all, "SGP4_Az_El_Doppler_All.csv");

T_visible = table( ...
    t_visible, az_visible, el_visible, range_visible, fd_visible, ...
    'VariableNames', {'Time_UTC','Azimuth_deg','Elevation_deg','Range_m','Doppler_Hz'} ...
);

writetable(T_visible, "SGP4_Az_El_Doppler_Visible.csv");

%% ---------- 9. Display pass summary ----------
fprintf("\nTotal samples: %d\n", length(t));
fprintf("Visible samples, El > 0: %d\n", length(t_visible));

if ~isempty(t_visible)
    fprintf("First visible time: %s UTC\n", string(t_visible(1)));
    fprintf("Last visible time : %s UTC\n", string(t_visible(end)));
    fprintf("Max elevation     : %.3f deg\n", max(el_visible));
else
    fprintf("No visible pass found in this time interval.\n");
end

%% ---------- 10. Plot ----------
figure;

subplot(3,1,1);
plot(t, az, 'LineWidth', 1.2);
grid on;
ylabel("Azimuth (deg)");
title("SGP4 Azimuth, Elevation, Doppler");

subplot(3,1,2);
plot(t, el, 'LineWidth', 1.2);
grid on;
ylabel("Elevation (deg)");
yline(0, "--");
    
subplot(3,1,3);
plot(t, fd, 'LineWidth', 1.2);
grid on;
ylabel("Doppler (Hz)");
xlabel("UTC Time");

%% ---------- 11. Plot visible pass only ----------
if ~isempty(t_visible)
    figure;

    subplot(3,1,1);
    plot(t, az_visible_plot, 'LineWidth', 1.2);
    grid on;
    ylabel("Azimuth (deg)");
    title("Visible Pass Only: El > 0");

    subplot(3,1,2);
    plot(t, el_visible_plot, 'LineWidth', 1.2);
    grid on;
    ylabel("Elevation (deg)");

    subplot(3,1,3);
    plot(t, fd_visible_plot, 'LineWidth', 1.2);
    grid on;
    ylabel("Doppler (Hz)");
    xlabel("UTC Time");
end
