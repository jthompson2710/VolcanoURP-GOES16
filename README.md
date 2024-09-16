# VolcanoURP-GOES16
 Algorithms to download and process GOES-16 data to analysis thermal trend at volcanoes


conda create -n GEO python=3.6 spyder=4.1.5 s3fs=0.4.2 netCDF4=1.4.2 Basemap=1.2.0 matplotlib=3.2.2 gdal=2.3.2 numpy jdcal scipy

/export/home/jot52/.conda/envs/GEO/bin/python3.6 /export/home/jot52/Documents/Python/UWI/UWIAutoDownProcATI16_Setup_Log.py

0 21 * * * /export/home/jot52/.conda/envs/GEO/bin/python3.6 /export/home/jot52/Documents/Python/UWI/UWIAutoDownProcATI16_Daily_Log.py
0 3 * * * /export/home/jot52/.conda/envs/GEO/bin/python3.6 /export/home/jot52/Documents/Python/UWI/UWIROIThermalTimeseries.py
