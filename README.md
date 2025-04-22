[![Ubuntu unit tests](https://github.com/MichaelClerx/bnglonglat/workflows/Ubuntu%20unit%20tests/badge.svg)](https://github.com/MichaelClerx/bnglonglat/actions?query=workflow%3A"Ubuntu+unit+tests")
[![Style](https://github.com/MichaelClerx/bnglonglat/workflows/Style/badge.svg)](https://github.com/MichaelClerx/bnglonglat/actions?query=workflow%3A"Style")

# British National Grid to Longitude/Lattitude

Pure python implementation of `convert_lonlat` from [lonlat_bng](https://github.com/urschrei/lonlat_bng).

Note: This uses the version without OSTN15 corrections, which can be several meters off.

## OSN guide

These calculations are based on [A Guide to Coordinate Systems in Great Britain](https://docs.os.uk/more-than-maps/deep-dive/a-guide-to-coordinate-systems-in-great-britain),
specifically [the equations here](https://docs.os.uk/more-than-maps/deep-dive/a-guide-to-coordinate-systems-in-great-britain/converting-between-grid-eastings-and-northings-and-ellipsoidal-latitude-and-longitude)
and constants from 
[here](https://docs.os.uk/more-than-maps/deep-dive/a-guide-to-coordinate-systems-in-great-britain/datum-ellipsoid-and-projection-information),
[here](https://docs.os.uk/more-than-maps/deep-dive/a-guide-to-coordinate-systems-in-great-britain/converting-between-3d-cartesian-and-ellipsoidal-latitude-longitude-and-height-coordinates), and
[here](https://docs.os.uk/more-than-maps/deep-dive/a-guide-to-coordinate-systems-in-great-britain/from-one-coordinate-system-to-another-geodetic-transformations/approximate-wgs84-to-osgb36-odn-transformation).

However, these pages were converted from a PDF and contain several errors and omissions.
The source code here was checked against [the PDF version](https://www.ordnancesurvey.co.uk/documents/resources/guide-coordinate-systems-great-britain.pdf).

## License

This code is available under an MIT license, See [license.txt](./license.txt).

It is based on [lonlat_bng](https://github.com/urschrei/lonlat_bng), which is Copyright (c) 2015 Stephan Hügel, and also published with an MIT license.
