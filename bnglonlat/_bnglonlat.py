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

MAX_EASTING = 700000
MAX_NORTHING = 1250000



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


175396

'''
get_ostn_ref(x, y):
    """ Try to get OSTN15 shift parameters, and calculate offsets. """
    key = x + (y * 701) + 1

    result = ostn15_lookup(&key)
    Ok((result.0, result.1, result.2))
}


def ostn15_shifts(x, y):
    """ Calculate OSTN15 shifts for a given coordinate. """
    # let e_index = (x / 1000.) as i32;
    e_index = x // 1000
    # let n_index = (y / 1000.) as i32;
    n_index = y // 1000

    # let x0 = e_index as i32 * 1000;
    x0 = e_index * 1000
    # let y0 = n_index as i32 * 1000;
    y0 = n_index * 1000

    # let s0: (f64, f64, f64) = get_ostn_ref(e_index, n_index)?;
    let s0: (f64, f64, f64) = get_ostn_ref(e_index, n_index)?;
    # let s1: (f64, f64, f64) = get_ostn_ref(e_index + 1, n_index)?;
    let s1: (f64, f64, f64) = get_ostn_ref(e_index + 1, n_index)?;
    # let s2: (f64, f64, f64) = get_ostn_ref(e_index, n_index + 1)?;
    let s2: (f64, f64, f64) = get_ostn_ref(e_index, n_index + 1)?;
    # let s3: (f64, f64, f64) = get_ostn_ref(e_index + 1, n_index + 1)?;
    let s3: (f64, f64, f64) = get_ostn_ref(e_index + 1, n_index + 1)?;

    # let dx = x - f64::from(x0);
    let dx = x - f64::from(x0);
    # let dy = y - f64::from(y0);
    let dy = y - f64::from(y0);

    # let t = dx / 1000.;
    let t = dx / 1000.;
    # let u = dy / 1000.;
    let u = dy / 1000.;

    # let f0 = (1. - t) * (1. - u);
    let f0 = (1. - t) * (1. - u);
    # let f1 = t * (1. - u);
    let f1 = t * (1. - u);
    # let f2 = (1. - t) * u;
    let f2 = (1. - t) * u;
    # let f3 = t * u;
    let f3 = t * u;

    # let se = f0 * s0.0 + f1 * s1.0 + f2 * s2.0 + f3 * s3.0;
    let se = f0 * s0.0 + f1 * s1.0 + f2 * s2.0 + f3 * s3.0;
    # let sn = f0 * s0.1 + f1 * s1.1 + f2 * s2.1 + f3 * s3.1;
    let sn = f0 * s0.1 + f1 * s1.1 + f2 * s2.1 + f3 * s3.1;
    # let sg = f0 * s0.2 + f1 * s1.2 + f2 * s2.2 + f3 * s3.2;

    return se, sn


'''

'''



/// Convert OSGB36 coordinates to Lon, Lat using OSTN15 data
#[allow(non_snake_case)]
pub fn convert_osgb36_to_ll(E: f64, N: f64) -> Result<(f64, f64), ()> {
    // Apply reverse OSTN15 adustments
    let epsilon = 0.009;
    let (mut dx, mut dy, _) = ostn15_shifts(E, N)?;
    let (mut x, mut y) = (E - dx, N - dy);
    let (mut last_dx, mut last_dy) = (dx, dy);
    let mut res;
    loop {
        res = ostn15_shifts(x, y)?;
        dx = res.0;
        dy = res.1;
        x = E - dx;
        y = N - dy;
        // If the difference […] is more than 0.00010m (User Guide, p15)
        // TODO: invert this logic
        if (dx - last_dx).abs() < epsilon && (dy - last_dy).abs() < epsilon {
            break;
        }
        last_dx = dx;
        last_dy = dy;
    }
    let x = (E - dx).round_to_mm();
    let y = (N - dy).round_to_mm();
    // We've converted to ETRS89, so we need to use the WGS84/ GRS80 ellipsoid constants
    convert_to_ll(x, y, GRS80_SEMI_MAJOR, GRS80_SEMI_MINOR)
}


// Easting and Northing to Lon, Lat conversion using a Helmert transform
// Note that either GRS80 or Airy 1830 ellipsoids can be passed
#[allow(non_snake_case)]
fn convert_to_ll(eastings: f64, northings: f64, ell_a: f64, ell_b: f64) -> Result<(f64, f64), ()> {
    // ensure that we're within the boundaries
    check(eastings, (0.000, MAX_EASTING))?;
    check(northings, (0.000, MAX_NORTHING))?;
    // ellipsoid squared eccentricity constant
    let a = ell_a;
    let b = ell_b;
    let e2 = (a.powi(2) - b.powi(2)) / a.powi(2);
    let n = (a - b) / (a + b);

    let dN = northings - N0;
    let mut phi = PHI0 + dN / (a * F0);
    let mut m = compute_m(phi, b, n);
    while (dN - m) >= 0.00001 {
        m = compute_m(phi, b, n);
        phi += (dN - m) / (a * F0);
    }
    let sp2 = phi.sin().powi(2);
    let nu = a * F0 * (1. - e2 * sp2).powf(-0.5);
    let rho = a * F0 * (1. - e2) * (1. - e2 * sp2).powf(-1.5);
    let eta2 = nu / rho - 1.;

    let tp = phi.tan();
    let tp2 = tp.powi(2);
    let tp4 = tp.powi(4);

    let VII = tp / (2. * rho * nu);
    let VIII = tp / (24. * rho * nu.powi(3)) * (5. + 3. * tp2 + eta2 - 9. * tp2 * eta2);
    let IX = tp / (720. * rho * nu.powi(5)) * (61. + 90. * tp2 + 45. * tp4);

    let sp = 1.0 / phi.cos();
    let tp6 = tp4 * tp2;

    let X = sp / nu;
    let XI = sp / (6. * nu.powi(3)) * (nu / rho + 2. * tp2);
    let XII = sp / (120. * nu.powi(5)) * (5. + 28. * tp2 + 24. * tp4);
    let XIIA = sp / (5040. * nu.powi(7)) * (61. + 662. * tp2 + 1320. * tp4 + 720. * tp6);

    let e = eastings - E0;

    phi = phi - VII * e.powi(2) + VIII * e.powi(4) - IX * e.powi(6);
    let mut lambda = LAM0 + X * e - XI * e.powi(3) + XII * e.powi(5) - XIIA * e.powi(7);

    phi = phi.to_degrees();
    lambda = lambda.to_degrees();
    return lambda, phi


// Intermediate calculation used for lon, lat to ETRS89 and reverse conversion
fn compute_m(phi: f64, b: f64, n: f64) -> f64 {
    let p_plus = phi + PHI0;
    let p_minus = phi - PHI0;

    b * F0
        * ((1. + n * (1. + 5. / 4. * n * (1. + n))) * p_minus
            - 3. * n * (1. + n * (1. + 7. / 8. * n)) * p_minus.sin() * p_plus.cos()
            + (15. / 8. * n * (n * (1. + n))) * (2. * p_minus).sin() * (2. * p_plus).cos()
            - 35. / 24. * n.powi(3) * (3. * p_minus).sin() * (3. * p_plus).cos())
}
'''
