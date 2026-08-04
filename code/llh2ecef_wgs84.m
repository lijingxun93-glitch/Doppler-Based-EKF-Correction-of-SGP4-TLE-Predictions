function ecef = llh2ecef_wgs84(llh)
% llh=[lat;lon;h], rad, rad, m
a = 6378137.0; f = 1/298.257223563; e2 = f*(2-f);
lat=llh(1); lon=llh(2); h=llh(3);
N = a/sqrt(1 - e2*sin(lat)^2);
ecef = [ (N+h)*cos(lat)*cos(lon);
         (N+h)*cos(lat)*sin(lon);
         (N*(1-e2)+h)*sin(lat) ];
end