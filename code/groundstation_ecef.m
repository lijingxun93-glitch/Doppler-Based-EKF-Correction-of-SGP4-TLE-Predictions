function [r_ecef, v_ecef] = groundstation_ecef(lat_deg, lon_deg, h_km)
    % WGS-84
    a = 6378.137;                  % km
    f = 1/298.257223563;
    e2 = f*(2-f);

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    sinlat = sin(lat);
    coslat = cos(lat);
    sinlon = sin(lon);
    coslon = cos(lon);

    N = a / sqrt(1 - e2*sinlat^2);

    x = (N + h_km) * coslat * coslon;
    y = (N + h_km) * coslat * sinlon;
    z = (N*(1 - e2) + h_km) * sinlat;

    r_ecef = [x; y; z];

    % Ground station velocity in ECEF due to Earth rotation
    omega_earth = 7.2921150e-5;    % rad/s
    omega_vec = [0; 0; omega_earth];

    v_ecef = cross(omega_vec, r_ecef);
end