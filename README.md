
dem_water DEM waterbody flattening
==================================

This tool processes a digital elevation model together with raster files
for waterbodies and modifies the DEM to flatten the waterbody areas while
keeping an overall consistent elevation model.

Warning
-------

Much of this program is highly experimental and not very robust.  It will 
crash or produce nonsense results when called with the wrong options.  Best 
to consider it just a demo for the moment.

Compiling the program
---------------------

dem_water requires the [CImg](http://cimg.sourceforge.net/) image 
processing library to be built and should be compilable on all platforms 
supported by the library.  You need to download this library separately.  
Copying CImg.h to the source directory is sufficient.

It also uses [GDAL](http://www.gdal.org/) for loading image files with 
geo-coordinates.

Program operation
-----------------

To correctly process waterbodies they need to be fully within the processing 
area of the program.  Normally this means a significantly larger area needs 
to be processed while the output can only be used for a smaller area.  There 
are two operating modes for this

* Tiled mode where a typical 1x1 degree DEM tile and its 8 neighbors are processed while the results are only for the central tile this is invoked by use of the `-is?` options.
* General mode using the `-i` option where the full area is only cropped by a single pixel and the output needs to be manually cropped further as needed.

Some of the input data is expected to cover the full processing area, some only the smaller output area.  This is mentioned in the options description.


Program options
---------------

the program offers the following command line options:

* `-i` DEM input file.  Full size (alternative to `-is?`)
* `-is1` to `-is9` DEM input 1x1 degree tiles (alternative to `-i`)
* `-m` Main waterbody mask image file, expected to be 32 bit signed integer with unique value for every waterbody values of -1 are interpreted as ocean.  Full size. (required)
* `-ms` Supersampled waterbody mask with byte valued pixels between 0 and 255.  Small size. (optional)
* `-b` Bathymetry DEM to use for ocean areas.  Small size. (optional)
* `-o` Output image file name (required)
* `-p` File for optput of POV-Ray mesh (optional)
* `-e` File for error point list in GeoJSON format (optional)
* `-ce` File for coastline error point list in GeoJSON format (optional)
* `-sr` DEM input resolution in tiled mode.  Either 1 or 3.  Default: `1`
* `-r` Radius to taper off elevation changes.  Default: 2.1
* `-z` Elevation value to use for the Ocean.  Default: 0
* `-lm` Minimum elevation value to allow for land.  Default: -32767
* `-w` Wall height to generate around water areas.  Default: 0
* `-t` Elevation difference threshold for the calculation of waterbody levels.  Default: `50`
* `-et` Elevation difference threshold for error point generation.  Default: `80`
* `-lon` tile longitude  Default: `0`
* `-lat` tile latitude  Default: `0`
* `-3` generate POV-Ray mesh in true 3d (i.e. curved earth surface).  Default: `off`
* `-oc` Generate output using CImg rather than in hgt format.  Default: `off`
* `-n` Noise to add to the output data (useful to mitigate integer meter discretization artefacts).  Default: `off`
* `-debug` Generate a large number of image files from intermediate steps in the current directory for debugging.  Default: off
* `-h` show available options

Using the program
-----------------

The program requires you to separately rasterize the waterbody data.  This can be done using something like the following:

`gdal_rasterize -ot Int32 -init -1 -burn 0 -ts ... -te ... $COASTLINE_SOURCE $WATER_MASK`

`gdal_rasterize -a "$LAKES_ID" $LAKES_SOURCE $WATER_MASK`

To generate the supersampled input mask just rasterize at higher resolution and scale down using `gdalwarp` or [ImageMagick](http://www.imagemagick.org/).

By default the program generates the output DEM in hgt format.  You can change that when using the general mode using `-oc`.

Error point output is pretty self explaining.  Attributes included are the ID from the input data and the maximum error.

If you try to use the POV-Ray mesh output you will almost certainly run into accuracy issues.

Legal stuff
-----------

This program is licensed under the GNU GPL version 3.

Copyright 2013 Christoph Hormann

