#!/usr/bin/env python3
#
# Pure Python implementation of the convertbng.util.convert_lonlat method.
#
# This itself is based on
#
# [1] https://docs.os.uk/more-than-maps/deep-dive/...
#     ...a-guide-to-coordinate-systems-in-great-britain/
#     ...converting-between-grid-eastings-and-northings-and-ellipsoidal-
#     ...latitude-and-longitude
#
# [2] https://docs.os.uk/more-than-maps/deep-dive/...
#     ...a-guide-to-coordinate-systems-in-great-britain/
#     ...datum-ellipsoid-and-projection-information
#
# [3] https://docs.os.uk/more-than-maps/deep-dive/...
#     ...a-guide-to-coordinate-systems-in-great-britain/converting-between-
#     ...3d-cartesian-and-ellipsoidal-latitude-longitude-and-height-coordinates
#
# [4] https://docs.os.uk/more-than-maps/deep-dive/...
#     ...a-guide-to-coordinate-systems-in-great-britain/
#     ...from-one-coordinate-system-to-another-geodetic-transformations/
#     ...approximate-wgs84-to-osgb36-odn-transformation
#
from numpy import pi, sin, cos, tan, arctan2, sqrt, abs, any


# Northing of true origin, from [1, 2]
N0 = -100000

# Easting of true origin, from [1, 2]
E0 = 400000

# Scale factor on central meridian, from [1]
F0 = 0.9996012717

# Latitude of true origin, from [2]
pio = pi / 180
phi0 = 49 * pio

# Longitude of true origina and central meridian, from [2]
lambda0 = -2 * pio

# Semi-major axis a (m), Airy 1830, from [2]
a = 6377563.396

# Semi major axis a (m), GRS80 aka WGS84 ellipsoid, from [2]
a_2 = 6378137

# Semi minor axis b (m), Airy 1830, from [2]
b = 6356256.909

# Semi minor axis b (m), GRS80 aka WGS84 ellipsoid, from [2]
b_2 = 6356752.3141

# The ellipsoid squared eccentricity constant e^2, from [3]
e2 = 1 - (b / a)**2
e2_2 = 1 - (b_2 / a_2)**2

# ETRS89 (WGS84) to OSGB36/ODN Helmert transformation parameters, from [4]
tx = -446.448
ty = 125.157
tz = -542.06
rxs = -0.1502
rys = -0.2470
rzs = -0.8421
cs = 20.4894 * 0.000001

# Equation n = (a - b) / (a + b), from [1]
n = (a - b) / (a + b)
n2 = n**2
n3 = n**3


def bnglonlat(easting, northing, tol1=1e-5, tol2=1e-16):
    """
    Converts a BNG easting and northing (in meters) to longitude and latitude.


    """
    # Part 1, following C.2 in the PDF: Converting eastings and northings to
    # (ellipsoidal) latitude and longitude

    # Iterative procedure to get phi
    M = 0
    phi = phi0
    while any(northing - N0 - M >= tol1):
        # First iteration: phi' = (N - N0) / (a * F0) + phi_0
        # Next iterations: phi_new = (N - N0 - M) / (a * F0) + phi'
        phi += (northing - N0 - M) / (a * F0)

        # Equations for M, from [1]
        psum = phi + phi0
        pdif = phi - phi0
        M = b * F0 * (
            + (1 + n + (5 / 4) * (n2 + n3)) * pdif
            - (3 * n + 3 * n2 + (21 / 8) * n3) * sin(pdif) * cos(psum)
            + (15 / 8) * (n2 + n3) * sin(2 * pdif) * cos(2 * psum)
            - (35 / 24) * n3 * sin(3 * pdif) * cos(3 * psum))

    # Calculate nu, rho, and eta^2, from [1]
    nu = a * F0 / sqrt(1 - e2 * sin(phi)**2)
    rho = a * F0 * (1 - e2) * (1 - e2 * sin(phi)**2)**(-1.5)
    eta2 = nu / rho - 1

    # Calculate equations VII, VIII, IX, X, XI, and XII, XIIA from [1]
    # But XII and XIIA taken from PDF version
    ta = tan(phi)
    ta2, ta4 = ta**2, ta**4
    nu3, nu5 = nu**3, nu**5
    sec = 1 / cos(phi)
    VII = ta / (2 * rho * nu)
    VIII = ta / (24 * rho * nu3) * (5 + 3 * ta2 + eta2 - 9 * ta2 * eta2)
    IX = ta / (720 * rho * nu5) * (61 + 90 * ta2 + 45 * ta4)
    X = sec / nu    # Typo in web version
    XI = sec / (6 * nu3) * (nu / rho + 2 * ta2)
    XII = sec / (120 * nu5) * (5 + 28 * ta2 + 24 * ta4)
    XIIA = sec / (5040 * nu**7) * (61 + 662 * ta2 + 1320 * ta4 + 720 * ta**6)

    # Calculate final phi and lambda
    dE = easting - E0
    phi = phi - VII * dE**2 + VIII * dE**4 - IX * dE**6
    lam = lambda0 + X * dE - XI * dE**3 + XII * dE**5 - XIIA * dE**7

    # Part 2, following B and B2 in the PDF, converting 3D Cartesian
    # coordinates to latitude, longitude and ellipsoid height

    # Equations B3, B4, B5, with H = 0 and replacing B2 with nu / F0 from the
    # previous section
    x1 = (nu / F0) * cos(phi) * cos(lam)
    y1 = (nu / F0) * cos(phi) * sin(lam)
    z1 = ((1 - e2) * nu / F0) * sin(phi)

    # Helmert transformation
    rx = -rxs * pio / 3600
    ry = -rys * pio / 3600
    rz = -rzs * pio / 3600
    x2 = -tx + (1 - cs) * x1 - rz * y1 + ry * z1
    y2 = -ty + rz * x1 + (1 - cs) * y1 - rx * z1
    z2 = -tz - ry * x1 + rx * y1 + (1 - cs) * z1

    # Equation B7
    p = sqrt(x2**2 + y2**2)

    # Iterative procedure with equations B2, B7 and B8
    phi = arctan2(z2, p * (1 - e2_2))
    old = 2 * pi
    while any(abs(phi - old) > tol2):
        old = phi
        nu = a_2 / sqrt(1 - e2_2 * sin(old)**2)
        phi = arctan2(z2 + e2_2 * nu * sin(old), p)

    # let mut lon = y_2.atan2(x_2);
    lam = arctan2(y2, x2)
    return lam / pio, phi / pio

