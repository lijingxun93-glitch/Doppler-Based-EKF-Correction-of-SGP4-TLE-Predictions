function [enu, enu_dot] = ecef2enu_vector(rho_ecef, rhodot_ecef, lat_deg, lon_deg)
    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    sinlat = sin(lat);
    coslat = cos(lat);
    sinlon = sin(lon);
    coslon = cos(lon);

    % ECEF -> ENU rotation
    R = [ -sinlon,            coslon,           0;
          -sinlat*coslon, -sinlat*sinlon,  coslat;
           coslat*coslon,  coslat*sinlon,  sinlat ];

    enu     = R * rho_ecef;
    enu_dot = R * rhodot_ecef;
end