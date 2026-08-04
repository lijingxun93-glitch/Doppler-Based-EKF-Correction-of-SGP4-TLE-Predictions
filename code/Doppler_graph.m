%--------------------------------------------------------------------------
%                   SGP4 Doppler generator
%

function [dopplerShift, azimuth, elevation, satdata] = Doppler_graph(startTime, dur, lat_deg, lon_deg, alt_m)

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
    Cnum = tline(3:7);      			        % Catalog Number (NORAD)
    SC   = tline(8);					        % Security Classification
    ID   = tline(10:17);			            % Identification Number
    year = str2num(tline(19:20));               % Year
    doy  = str2num(tline(21:32));               % Day of year
    epoch = str2num(tline(19:32));              % Epoch
    TD1   = str2num(tline(34:43));              % first time derivative
    TD2   = str2num(tline(45:50));              % 2nd Time Derivative
    ExTD2 = str2num(tline(51:52));              % Exponent of 2nd Time Derivative
    BStar = str2num(tline(54:59));              % Bstar/drag Term
    ExBStar = str2num(tline(60:61));            % Exponent of Bstar/drag Term
    BStar = BStar*1e-5*10^ExBStar;
    Etype = tline(63);                          % Ephemeris Type
    Enum  = str2num(tline(65:end));             % Element Number
    
    % read second line
    tline = fgetl(fid);
    i = str2num(tline(9:16));                   % Orbit Inclination (degrees)
    raan = str2num(tline(18:25));               % Right Ascension of Ascending Node (degrees)
    e = str2num(strcat('0.',tline(27:33)));     % Eccentricity
    omega = str2num(tline(35:42));              % Argument of Perigee (degrees)
    M = str2num(tline(44:51));                  % Mean Anomaly (degrees)
    no = str2num(tline(53:63));                 % Mean Motion
    a = ( ge/(no*2*pi/86400)^2 )^(1/3);         % semi major axis (km)
    rNo = str2num(tline(65:end));               % Revolution Number at Epoch
    
    fclose(fid);
    
    satdata.epoch = epoch;
    satdata.norad_number = Cnum;
    satdata.bulletin_number = ID;
    satdata.classification = SC; % almost always 'U'
    satdata.revolution_number = rNo;
    satdata.ephemeris_type = Etype;
    satdata.xmo = M * (pi/180);
    satdata.xnodeo = raan * (pi/180);
    satdata.omegao = omega * (pi/180);
    satdata.xincl = i * (pi/180);
    satdata.eo = e;
    satdata.xno = no * TWOPI / MINUTES_PER_DAY;
    satdata.xndt2o = TD1 * TWOPI / MINUTES_PER_DAY_SQUARED;
    satdata.xndd6o = TD2 * 10^ExTD2 * TWOPI / MINUTES_PER_DAY_CUBED;
    satdata.bstar = BStar;
    
    % read Earth orientation parameters
    fid = fopen('EOP-All.txt','r');
    %  ----------------------------------------------------------------------------------------------------
    % |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
    % |(0h UTC)           "         "          s          s          "        "          "         "     s 
    %  ----------------------------------------------------------------------------------------------------
    while ~feof(fid)
        tline = fgetl(fid);
        k = strfind(tline,'NUM_OBSERVED_POINTS');
        if (k == 1)
            numrecsobs = str2num(tline(21:end));
            tline = fgetl(fid);
            for i=1:numrecsobs
                eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
            end
            for i=1:4
                tline = fgetl(fid);
            end
            numrecspred = str2num(tline(22:end));
            tline = fgetl(fid);
            for i=numrecsobs+1:numrecsobs+numrecspred
                eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
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
    
    %------
    gs.lat_deg  = lat_deg;    % deg (+N)
    gs.lon_deg  = lon_deg;    % deg (+E)
    gs.alt_m    = alt_m;         % m
    gs_ecef = llh2ecef_wgs84([deg2rad(gs.lat_deg); deg2rad(gs.lon_deg); gs.alt_m])/1000;
    
    duration_min = dur;
    % startTime = datetime(2025,10,31,13,15,11,'TimeZone','UTC');
    % startTime = datetime(2026,04,19,4,12,00,'TimeZone','UTC');
    t_epoch = datetime(year,mon,day,hr,minute,sec,'TimeZone','UTC');
    t_start = startTime-t_epoch;
    
    t = (0:duration_min*60)/60;
    fc = 100e6; c = 299792.458;
    dopplerShift = zeros(length(t),1);
    azimuth = zeros(length(t),1);
    elevation = zeros(length(t),1);
    j = 1;
    for tsince = t % amount of time in which you are going to propagate satellite's state vector forward (+) or backward (-) [minutes] 
    
        MJD_UTC = MJD_Epoch + (minutes(t_start) + tsince)/1440;
        
        % Earth Orientation Parameters
        [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
        [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
        MJD_UT1 = MJD_UTC + UT1_UTC/86400;
        MJD_TT  = MJD_UTC + TT_UTC/86400;
        T = (MJD_TT-const.MJD_J2000)/36525;
        
        [rteme, vteme] = sgp4(minutes(t_start) + tsince, satdata);
        [recef,vecef] = teme2ecef(rteme,vteme,T,MJD_UT1+2400000.5,LOD,x_pole,y_pole,2);
    %-------------------------doppler-------------------------
        % rho = recef - gs_ecef;
        % rho_unit = rho / norm(rho);
        % range_rate = dot(vecef, rho_unit); %径向速度
        % dopplerShift(j) = - (range_rate / c) * fc;
    %-----------------------azimuth and elevation-------------------------
        % Ground station ECEF position/velocity
        [rsite_ecef, vsite_ecef] = groundstation_ecef(gs.lat_deg, gs.lon_deg, gs.alt_m/1000);
        
        % Relative position/velocity in ECEF
        rho_ecef    = recef - rsite_ecef; % exactly the same as rho
        rhodot_ecef = vecef - vsite_ecef;
        
        % Convert relative vector to local ENU frame
        [enu, enu_dot] = ecef2enu_vector(rho_ecef, rhodot_ecef, gs.lat_deg, gs.lon_deg);
        
        E  = enu(1);
        N  = enu(2);
        U  = enu(3);
        
        Ed = enu_dot(1);
        Nd = enu_dot(2);
        Ud = enu_dot(3);
        
        % Range
        range_km = norm(enu);
        
        % Elevation
        el = asin(U / range_km);
        elevation(j) = el * 180/pi;
        
        % Azimuth (clockwise from north)
        az = atan2(E, N);
        if az < 0
            az = az + 2*pi;
        end
        azimuth(j) = az * 180/pi;
     %-------------------------doppler-------------------------
        rho_unit = rho_ecef / norm(rho_ecef);
        range_rate = dot(rhodot_ecef, rho_unit); %径向速度
        dopplerShift(j) = - (range_rate / c) * fc;
    
        j = j+1;
    end
    % % dopplerRate = [dopplerShift;0] - [0;dopplerShift];
    % % dopplerRate(1) = []; dopplerRate(end) = 0;
    % tplot = minutes(t)+startTime;
    % idx = find(elevation > 0); 
    % figure;
    % subplot(2,1,1);
    % plot(tplot(idx),dopplerShift(idx));
    % legend('Gpredict','SGP4');
    % % subplot(2,2,2);
    % % plot(tplot(idx),dopplerRate(idx));
    % % legend('Gpredict','SGP4');
    % subplot(2,2,3);
    % plot(tplot(idx),azimuth(idx));
    % legend('Gpredict','SGP4');
    % subplot(2,2,4);
    % plot(tplot(idx),elevation(idx));
    % legend('Gpredict','SGP4');
end

