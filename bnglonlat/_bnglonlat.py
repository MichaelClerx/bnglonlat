#!/usr/bin/env python3
#
# Pure Python implementation of the convertbng.util.convert_lonlat method.
#
from numpy import pi, sin, cos, tan, arctan2, sqrt, abs, any

#/// assert_eq!((-0.328248, 51.44534), convert_lonlat(&516276, &173141));


def _curvature(a, f0, e2, lat):
    """
    Calculate the meridional radius of curvature"""
    #a * f0 * (1. - e2) * (1. - e2 * lat.sin().powi(2)).powf(-1.5)
    return a * f0 * (1 - e2) * (1 - e2 * sin(lat)**2)**(-1.5)


a = 6377563.396
b = 6356256.909
a_2 = 6378137.000
b_2 = 6356752.3141
E0 = 400000
N0 = -100000
F0 = 0.9996012717
TX = -446.448
TY = 125.157
TZ = -542.060
RXS = -0.1502
RYS = -0.2470
RZS = -0.8421
CS = 20.4894 * 0.000001
minus_s = -CS
pio = pi / 180

# let lat0 = 49. * PI / 180.;
lat0 = 49 * pio
# let lon0 = -2. * PI / 180.;
lon0 = -2 * pio

# let e2 = 1. - b.powi(2) / a.powi(2);
e2 = 1 - (b / a)**2

# let n = (a - b) / (a + b);
n = (a - b) / (a + b)
n2 = n**2
n3 = n**3


def bnglonlat(easting, northing, tol1=1e-5, tol2=1e-16):
    """
    Converts a BNG easting and northin (in meters) to longitude and lattitude.
    """
    # let mut lat = lat0;
    lat = lat0
    # let mut M: f64 = 0.0;
    M = 0
    # while (northing - N0 - M) >= 0.00001 {
    while any(northing - N0 - M >= 1e-5):
        # lat += (northing - N0 - M) / (a * F0);
        lat += (northing - N0 - M) / (a * F0)

        # let M1 = (1. + n + (5. / 4.) * n.powi(3) + (5. / 4.) * n.powi(3))
        # * (lat - lat0);
        M1 = (1 + n + (5 / 4) * (n2 + n3)) * (lat - lat0)

        # let M2 = (3. * n + 3. * n.powi(2) + (21. / 8.) * n.powi(3))
        #     * ((lat.sin() * lat0.cos()) - (lat.cos() * lat0.sin()))
        #         .ln_1p()
        #         .exp_m1() * (lat + lat0).cos();
        M2 = (3 * n + 3 * n2 + (21 / 8) * n3) * cos(lat + lat0) * (
            sin(lat) * cos(lat0) - cos(lat) * sin(lat0))

        # let M3 = ((15. / 8.) * n.powi(2) + (15. / 8.) * n.powi(3))
        #     * (2. * (lat - lat0)).sin()
        #     * (2. * (lat + lat0)).cos();
        M3 = ((15 / 8) * (n2 + n3) * sin(2 * (lat - lat0))
              * cos(2 * (lat + lat0)))

        # let M4 = (35. / 24.) * n.powi(3) * (3. * (lat - lat0)).sin()
        # * (3. * (lat + lat0)).cos();
        M4 = (35 / 24) * n3 * sin(3 * (lat - lat0)) * cos(3 * (lat + lat0))

        # M = b * F0 * (M1 - M2 + M3 - M4);
        M = b * F0 * (M1 - M2 + M3 - M4)

    # let nu = a * F0 / (1. - e2 * lat.sin().powi(2)).sqrt();
    nu = a * F0 / sqrt(1 - e2 * sin(lat)**2)

    # Meridional radius of curvature
    # let rho = curvature(a, F0, e2, lat);
    rho = _curvature(a, F0, e2, lat)
    # let eta2 = nu / rho - 1.;
    eta2 = nu / rho - 1

    ta = tan(lat)
    ta2 = ta**2
    ta4 = ta**4
    nu3 = nu**3
    nu5 = nu**5
    # let secLat = 1. / lat.cos();
    se = 1 / cos(lat)

    # let VII = lat.tan() / (2. * rho * nu);
    VII = ta / (2 * rho * nu)
    # let VIII = lat.tan() / (24. * rho * nu.powi(3))
    #   * (5. + 3. * lat.tan().powi(2) + eta2 - 9. * lat.tan().powi(2) * eta2);
    VIII = ta / (24 * rho * nu3) * (5 + 3 * ta2 + eta2 - 9 * ta2 * eta2)
    # let IX = lat.tan() / (720. * rho * nu.powi(5))
    #              * (61. + 90. * lat.tan().powi(2) + 45. * lat.tan().powi(4));
    IX = ta / (720 * rho * nu5) * (61 + 90 * ta2 + 45 * ta4)
    # let X = secLat / nu;
    X = se / nu
    # let XI = secLat / (6. * nu.powi(3)) * (nu / rho + 2. * lat.tan().powi(2))
    XI = se / (6 * nu3) * (nu / rho + 2 * ta2)
    # let XII = secLat / (120. * nu.powi(5)) * (5. + 28. * lat.tan().powi(2)
    #                                        + 24. * lat.tan().powi(4));
    XII = se / (120 * nu5) * (5 + 28 * ta2 + 24 * ta4)
    # let XIIA = secLat / (5040. * nu.powi(7))
    #     * (61. + 662. * lat.tan().powi(2)
    #        + 1320. * lat.tan().powi(4) + 720. * lat.tan().powi(6));
    XIIA = se / (5040 * nu**7) * (61 + 662 * ta2 + 1320 * ta4 + 720 * ta**6)
    # let dE = easting - E0;
    dE = easting - E0

    # let lat_1 = lat - VII * dE.powi(2) + VIII * dE.powi(4) - IX * dE.powi(6);
    lat_1 = lat - VII * dE**2 + VIII * dE**4 - IX * dE**6
    # let lon_1 = lon0 + X * dE - XI * dE.powi(3) + XII * dE.powi(5)
    #                                                - XIIA * dE.powi(7);
    lon_1 = lon0 + X * dE - XI * dE**3 + XII * dE**5 - XIIA * dE**7

    # let H = 0.;
    H = 0
    # let x_1 = (nu / F0 + H) * lat_1.cos() * lon_1.cos();
    x_1 = (nu / F0 + H) * cos(lat_1) * cos(lon_1)
    # let y_1 = (nu / F0 + H) * lat_1.cos() * lon_1.sin();
    y_1 = (nu / F0 + H) * cos(lat_1) * sin(lon_1)
    # let z_1 = ((1. - e2) * nu / F0 + H) * lat_1.sin();
    z_1 = ((1 - e2) * nu / F0 + H) * sin(lat_1)

    # let tx = TX.abs();
    tx = -TX
    # let ty = TY * -1.;
    ty = -TY
    # let tz = TZ.abs();
    tz = -TZ

    # let rxs = RXS * -1.;
    rxs = -RXS
    # let rys = RYS * -1.;
    rys = -RYS
    # let rzs = RZS * -1.;
    rzs = -RZS

    # let rx = rxs * PI / (180. * 3600.);
    rx = rxs * pio / 3600
    # let ry = rys * PI / (180. * 3600.);
    ry = rys * pio / 3600
    # let rz = rzs * PI / (180. * 3600.); // In radians
    rz = rzs * pio / 3600
    # let x_2 = tx + (1. + minus_s) * x_1 + (-rz) * y_1 + (ry) * z_1;
    x_2 = tx + (1 + minus_s) * x_1 - rz * y_1 + ry * z_1
    # let y_2 = ty + (rz) * x_1 + (1. + minus_s) * y_1 + (-rx) * z_1;
    y_2 = ty + rz * x_1 + (1 + minus_s) * y_1 - rx * z_1
    # let z_2 = tz + (-ry) * x_1 + (rx) * y_1 + (1. + minus_s) * z_1;
    z_2 = tz - ry * x_1 + rx * y_1 + (1 + minus_s) * z_1

    # let e2_2 = 1. - b_2.powi(2) / a_2.powi(2);
    e2_2 = 1 - (b_2 / a_2)**2

    # let p = (x_2.powi(2) + y_2.powi(2)).sqrt();
    p = sqrt(x_2**2 + y_2**2)

    # Lat is obtained by iterative procedure
    # Initial value
    # let mut lat = z_2.atan2(p * (1. - e2_2));
    lat = arctan2(z_2, p * (1 - e2_2))
    # let mut latold = 2. * PI;
    old = 2 * pi
    # while (lat - latold).abs() > (10. as f64).powi(-16) {
    while any(abs(lat - old) > tol2):
        # mem::swap(&mut lat, &mut latold);
        old = lat
        # nu_2 = a_2 / (1. - e2_2 * latold.sin().powi(2)).sqrt();
        nu_2 = a_2 / sqrt(1 - e2_2 * sin(old)**2)
        # lat = (z_2 + e2_2 * nu_2 * latold.sin()).atan2(p);
        lat = arctan2(z_2 + e2_2 * nu_2 * sin(old), p)

    # let mut lon = y_2.atan2(x_2);
    lon = arctan2(y_2, x_2)
    # lat = lat * 180. / PI;
    lat /= pio
    # lon = lon * 180. / PI;
    lon /= pio
    return lon, lat

