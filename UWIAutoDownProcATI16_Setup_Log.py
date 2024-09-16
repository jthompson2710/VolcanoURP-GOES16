#----------------------------------------------------------------------------#
# Created on Mon Dec 14 13:12:33 2020
# Last edit September 14, 2024
# James Thompson - University of Texas at Austin, University of Pittsburgh
# For use by the University of West Indies
#----------------------------------------------------------------------------#
# Notes:
# This is only stable using Python 3.6 (do not use 3.7).
# Stable and tested using Spyder 4.1.5 with conda and python 3.6.12. 
# Please use the function versions listed below (some require back tracking)
#----------------------------------------------------------------------------#
# Functions
import s3fs # version 0.4.2
from netCDF4 import Dataset # version 1.4.2
import numpy as np # version 1.19.2
from numpy import zeros_like
import os # version 
import csv # version 
import shutil # version
from remap_rad import remap_rad # Import the Remap function
from remap_temp import remap_temp # Import the Remap function
from osgeo import gdal # version 2.3.2
from TES_NOSKY_ABI import *
import TES_NOSKY_ABI
import jdcal
import datetime
#----------------------------------------------------------------------------#
################################Start of Script###############################
#----------------------------------------------------------------------------#
## standards and constant variables

# image resolution (the higher the number the faster the processing is)
resolution = 2.0 # Full resolution

#Empty arrays
all_emisf = []
all_Ts = []

#cm and um band conversions
BdEq = np.array([[37.0,0.2],[56.0,0.4],[42.0,0.38],[47.0,0.5],[64.0,0.8],[66.0,1.0],[34.0,0.6]])

# Use the anonymous credentials to access public data
fs = s3fs.S3FileSystem(anon=True)

#----------------------------------------------------------------------------#
## Import logs and email data

logdir = './Logs/'
datadir = './Data/'
outdir = './Data/'

with open(logdir+'UWI_GEO_active_volcanoes.csv', newline='') as f:
    reader = csv.reader(f)
    volcano_list = list(reader)
volcano_list = volcano_list[1:]

#Setup start and end times
overall_volcano = volcano_list[0]
erupt_yr = int(overall_volcano[3])
erupt_JD = int(overall_volcano[4]) #-1 # take it back 1 day from email time stamp
erupt_hour = float(overall_volcano[5])
days_prior = float(overall_volcano[6])
days_post = float(overall_volcano[7])
erupt_JDdate = jdcal.gcal2jd(erupt_yr,1,erupt_JD)
start_JDdate = erupt_JDdate[1] - days_prior
Jdfmt = '(%Y, %m, %d, %H.%M)'
start_day = datetime.datetime.strptime(str(jdcal.jd2gcal(2400000.5,start_JDdate)), Jdfmt).timetuple().tm_yday
start_year = datetime.datetime.strptime(str(jdcal.jd2gcal(2400000.5,start_JDdate)), Jdfmt).timetuple().tm_year
end_JDdate = erupt_JDdate[1] + days_post
end_day = datetime.datetime.strptime(str(jdcal.jd2gcal(2400000.5,end_JDdate)), Jdfmt).timetuple().tm_yday
end_year = datetime.datetime.strptime(str(jdcal.jd2gcal(2400000.5,end_JDdate)), Jdfmt).timetuple().tm_year
total_days = end_JDdate - start_JDdate
total_hours = int(total_days * 24.0)
# initial extract day and time
i_dir_year=int(start_year)
i_dir_day=int(start_day)
i_dir_hour=int(erupt_hour)

#Determine reference data to check all data has been processed
#collect information on final valocano in monitoring list
final_volcano = volcano_list[-1]
final_volcano_name = final_volcano[0]

#%%
#----------------------------------------------------------------------------#
## START DOWNLAOD LOOP ## 
## download GOES ABI data ##
#----------------------------------------------------------------------------#

# Determine if in GEOS 16 data extent and then do
satellite = 'goes16'
l1_dir = datadir+'FullDisk/L1/'

## build loops
#determine the day and hour of folder to extract from
for i in range(0,total_hours):
    if (i_dir_day < 100) and (i_dir_day > 9):
        s_i_dir_day = '0'+str(i_dir_day)
    elif i_dir_day < 10:
        s_i_dir_day = '00'+str(i_dir_day)
    else:
        s_i_dir_day = str(i_dir_day)
    
    if i_dir_hour < 10:
        s_i_dir_hour = '0'+str(i_dir_hour)
    else:
        s_i_dir_hour = str(i_dir_hour)
    
    #detemines decimal day
    cJdfmt = "(%Y, %m, %d, %H.%M)"
    cJDdate = jdcal.gcal2jd(i_dir_year,1,i_dir_day)
    c_timestamp = datetime.datetime.strptime(str(jdcal.jd2gcal(cJDdate[0],cJDdate[1])), cJdfmt).timestamp()
    #datetime.datetime.fromtimestamp(c_timestamp)
    #if satellite == 'goes16':
    if (c_timestamp >= 1513663200.0): 
        f_remap_rad = remap_rad
        f_remap_temp = remap_temp
    else:
        print("no GOES 16 data available for this date ("+str(c_timestamp)+")")
        i_dir_hour = i_dir_hour + 1
        #move forward one day if required
        if (jdcal.is_leap(i_dir_year)!=False) and (i_dir_hour > 23) and (i_dir_day > 365):
            i_dir_day = 1
            i_dir_hour = 0
            i_dir_year = i_dir_year + 1
            print("Done day "+str(i_dir_day-1)+", "+str(i_dir_year-1))
        if (jdcal.is_leap(i_dir_year)==False) and (i_dir_hour > 23) and (i_dir_day > 364):
            i_dir_day = 1
            i_dir_hour = 0
            i_dir_year = i_dir_year + 1
            print("Done day "+str(i_dir_day-1)+", "+str(i_dir_year-1))
        if i_dir_hour > 23:
            i_dir_day = i_dir_day + 1
            i_dir_hour = 0
            print("Done day "+str(i_dir_day-1)+", "+str(i_dir_year))
        continue

#----------------------------------------------------------------------------#
    #determine if data has already been downloaded and L2 processed (final working_volcano_name)
    if os.path.exists(outdir+final_volcano_name+'/L2/TES/Emi/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'):
        #move time forward one hour
        i_dir_hour = i_dir_hour + 1
        #move forward one day if required
        if (jdcal.is_leap(i_dir_year)!=False) and (i_dir_hour > 23) and (i_dir_day > 365):
            i_dir_day = 1
            i_dir_hour = 0
            i_dir_year = i_dir_year + 1
            print("Done day "+str(i_dir_day-1)+", "+str(i_dir_year-1))
        if (jdcal.is_leap(i_dir_year)==False) and (i_dir_hour > 23) and (i_dir_day > 364):
            i_dir_day = 1
            i_dir_hour = 0
            i_dir_year = i_dir_year + 1
            print("Done day "+str(i_dir_day-1)+", "+str(i_dir_year-1))
        if i_dir_hour > 23:
            i_dir_day = i_dir_day + 1
            i_dir_hour = 0
            print("Done day "+str(i_dir_day-1)+", "+str(i_dir_year))
        continue
        
    #build directory
    if not os.path.exists(l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'):
        os.makedirs(l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/')
    #determine the radiance data files to extract in each folder
    rad_files = np.array(fs.ls('s3://noaa-'+satellite+'/ABI-L1b-RadF/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'))
    TIRbd1 = list(filter(lambda x: 'C10' in x, rad_files))
    TIRbd2 = list(filter(lambda x: 'C11' in x, rad_files))
    TIRbd3 = list(filter(lambda x: 'C12' in x, rad_files))
    TIRbd4 = list(filter(lambda x: 'C13' in x, rad_files))
    TIRbd5 = list(filter(lambda x: 'C14' in x, rad_files))
    TIRbd6 = list(filter(lambda x: 'C15' in x, rad_files))
    TIRbd7 = list(filter(lambda x: 'C16' in x, rad_files))

    #extract radiance data in folder 
    for i in range(len(TIRbd1)):
        if TIRbd1[i].size:
            fs.get(TIRbd1[i], l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd1[i].split('/')[-1])
    for i in range(len(TIRbd2)):
        if TIRbd2[i].size:
            fs.get(TIRbd2[i], l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd2[i].split('/')[-1])
    for i in range(len(TIRbd3)):
        if TIRbd3[i].size:
            fs.get(TIRbd3[i], l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd3[i].split('/')[-1])
    for i in range(len(TIRbd4)):
        if TIRbd4[i].size:
            fs.get(TIRbd4[i], l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd4[i].split('/')[-1])
    for i in range(len(TIRbd5)):
        if TIRbd5[i].size:
            fs.get(TIRbd5[i], l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd5[i].split('/')[-1])
    for i in range(len(TIRbd6)):
        if TIRbd6[i].size:
            fs.get(TIRbd6[i], l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd6[i].split('/')[-1])
    for i in range(len(TIRbd7)):
        if TIRbd7[i].size:
            fs.get(TIRbd7[i], l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd7[i].split('/')[-1]) 
    
#----------------------------------------------------------------------------#
## END DOWNLAOD SECTION ## 
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
## START INDIVIDUAL VOLCANO PROCESSING LOOP ## 
## extract and process all UWI volcanoes ##
#----------------------------------------------------------------------------#

    ## build loops
    #determine the day and hour of folder to extract from
    for v in range(0,len(volcano_list)):
        #individaul input variables
        working_volcano = volcano_list[v]
        working_volcano_name = working_volcano[0]  #'Klyuchevskoy' # time zone is utc+12
        working_coordinates = working_volcano[1] 
        working_newextent_range = working_volcano[2]
    
        # working_volcano_name location and extent of extraction
        lat_coord = float( working_coordinates.split(', ')[0] )
        long_coord = float( working_coordinates.split(', ')[1] )
        extent_range = [ float(working_newextent_range.split(',')[0]), float(working_newextent_range.split(',')[1]), float(working_newextent_range.split(',')[2]), float(working_newextent_range.split(',')[3]) ]
        extent = [long_coord+extent_range[1], lat_coord+extent_range[0], long_coord+extent_range[3], lat_coord+extent_range[2]] #min long, min lat, max long, max lat 
        
        # Determine if in GEOS 16 data extent and then do
        if long_coord < -5.0 and long_coord > -105.0:
            satellite = 'goes16'
            l1_dir = datadir+'FullDisk/L1/'
            L1b_dir = outdir+working_volcano_name+'/L1b/'
            l2_dir = outdir+working_volcano_name+'/L2/'
        
        # Determine if not in GEOS  data extent and then break --> continue to next iteration
        else:
            satellite = 'error'
            raise Exception('Error with geographical coverage - '+working_volcano_name+' is outside satellite coverage')
            continue
        
        #--------------------------------------------------------------------#
        #determine if data has already been downloaded and L2 processed
        if os.path.exists(l2_dir+'TES/Emi/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'):
            continue
            
        # geospacial subset and georeference data
        #determine file just downloaded and only count nc files
        rad_files = os.listdir(l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/')
        TIRbd0 = list(filter(lambda x: 'C10' in x, rad_files))
        TIRbd0=[f for f in sorted(TIRbd0) if (str(f))[-2:] == "nc"]
        TIRbd1 = list(filter(lambda x: 'C11' in x, rad_files))
        TIRbd1=[f for f in sorted(TIRbd1) if (str(f))[-2:] == "nc"]
        TIRbd2 = list(filter(lambda x: 'C12' in x, rad_files))
        TIRbd2=[f for f in sorted(TIRbd2) if (str(f))[-2:] == "nc"]
        TIRbd3 = list(filter(lambda x: 'C13' in x, rad_files))
        TIRbd3=[f for f in sorted(TIRbd3) if (str(f))[-2:] == "nc"]
        TIRbd4 = list(filter(lambda x: 'C14' in x, rad_files))
        TIRbd4=[f for f in sorted(TIRbd4) if (str(f))[-2:] == "nc"]
        TIRbd5 = list(filter(lambda x: 'C15' in x, rad_files))
        TIRbd5=[f for f in sorted(TIRbd5) if (str(f))[-2:] == "nc"]
        TIRbd6 = list(filter(lambda x: 'C16' in x, rad_files))
        TIRbd6=[f for f in sorted(TIRbd6) if (str(f))[-2:] == "nc"]
        
        #--------------------------------------------------------------------#
        # Call the reprojection funcion and reproject for working_volcano_name area
        for j in range(len(TIRbd0)):
            
            if (len(TIRbd0) == len(TIRbd1) and len(TIRbd0) == len(TIRbd2) and len(TIRbd0) == len(TIRbd3) and len(TIRbd0) == len(TIRbd4) and len(TIRbd0) == len(TIRbd5) and len(TIRbd0) == len(TIRbd6)):  
                
                #------------------------------------------------------------#
                #make output directory
                path = TIRbd0[j]
                # Search for the Scan start in the file name
                Start = (path[path.find("s")+1:path.find("_e")])
                # Search for the Scan end in the file name
                End = (path[path.find("e")+1:path.find("_c")])
                sepoutdir=L1b_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/Sep/'+Start[0:13]+'_'+End[0:13]+'/'
                
                if not os.path.exists(sepoutdir):
                    os.makedirs(sepoutdir)
                
                #------------------------------------------------------------#
                ##band 10
                grid0 = f_remap_rad(l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd0[j], extent, resolution, 'HDF5')
                #convert to wm-2sr-2um-2
                grid0cm = grid0.GetRasterBand(1).ReadAsArray()
                grid0cm = grid0cm * BdEq[0,0] / BdEq[0,1] / 1000.0
                grid0.GetRasterBand(1).WriteArray(grid0cm)              
                # Export the result to GeoTIFF
                driver = gdal.GetDriverByName('GTiff')
                driver.CreateCopy(sepoutdir+'TIR1.tif', grid0, 0)
                
                ##band 11
                grid1 = f_remap_rad(l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd1[j], extent, resolution, 'HDF5')
                #convert to wm-2sr-2um-2
                grid1cm = grid1.GetRasterBand(1).ReadAsArray()
                grid1cm = grid1cm * BdEq[1,0] / BdEq[1,1] / 1000.0
                grid1.GetRasterBand(1).WriteArray(grid1cm)   
                # Export the result to GeoTIFF
                driver = gdal.GetDriverByName('GTiff')
                driver.CreateCopy(sepoutdir+'TIR2.tif', grid1, 0)
                
                ##band 12
                grid2 = f_remap_rad(l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd2[j], extent, resolution, 'HDF5')
                #convert to wm-2sr-2um-2
                grid2cm = grid2.GetRasterBand(1).ReadAsArray()
                grid2cm = grid2cm * BdEq[2,0] / BdEq[2,1] / 1000.0
                grid2.GetRasterBand(1).WriteArray(grid2cm)   
                # Export the result to GeoTIFF
                driver = gdal.GetDriverByName('GTiff')
                driver.CreateCopy(sepoutdir+'TIR3.tif', grid2, 0)
                
                ##band 13
                grid3 = f_remap_rad(l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd3[j], extent, resolution, 'HDF5')
                #convert to wm-2sr-2um-2
                grid3cm = grid3.GetRasterBand(1).ReadAsArray()
                grid3cm = grid3cm * BdEq[3,0] / BdEq[3,1] / 1000.0
                grid3.GetRasterBand(1).WriteArray(grid3cm)   
                # Export the result to GeoTIFF
                driver = gdal.GetDriverByName('GTiff')
                driver.CreateCopy(sepoutdir+'TIR4.tif', grid3, 0)
                
                ##band 14
                grid4 = f_remap_rad(l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd4[j], extent, resolution, 'HDF5')
                #convert to wm-2sr-2um-2
                grid4cm = grid4.GetRasterBand(1).ReadAsArray()
                grid4cm = grid4cm * BdEq[4,0] / BdEq[4,1] / 1000.0
                grid4.GetRasterBand(1).WriteArray(grid4cm)   
                # Export the result to GeoTIFF
                driver = gdal.GetDriverByName('GTiff')
                driver.CreateCopy(sepoutdir+'TIR5.tif', grid4, 0)
                
                ##band 15
                grid5 = f_remap_rad(l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd5[j], extent, resolution, 'HDF5')
                #convert to wm-2sr-2um-2
                grid5cm = grid5.GetRasterBand(1).ReadAsArray()
                grid5cm = grid5cm * BdEq[5,0] / BdEq[5,1] / 1000.0
                grid5.GetRasterBand(1).WriteArray(grid5cm)   
                # Export the result to GeoTIFF
                driver = gdal.GetDriverByName('GTiff')
                driver.CreateCopy(sepoutdir+'TIR6.tif', grid5, 0)
                
                ##band 16
                grid6 = f_remap_rad(l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'+TIRbd6[j], extent, resolution, 'HDF5')
                #convert to wm-2sr-2um-2
                grid6cm = grid6.GetRasterBand(1).ReadAsArray()
                grid6cm = grid6cm * BdEq[6,0] / BdEq[6,1] / 1000.0
                grid6.GetRasterBand(1).WriteArray(grid6cm)   
                # Export the result to GeoTIFF
                driver = gdal.GetDriverByName('GTiff')
                driver.CreateCopy(sepoutdir+'TIR7.tif', grid6, 0)
                
                #------------------------------------------------------------#
                ##combine
                outputfilename = L1b_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/Rad_TIR_'+Start[0:13]+'_'+End[0:13]
                vrt_options = gdal.BuildVRTOptions(separate=True)
                gdal.BuildVRT(outputfilename+'.vrt', [sepoutdir+'TIR1.tif', sepoutdir+'TIR2.tif', sepoutdir+'TIR3.tif', sepoutdir+'TIR4.tif', sepoutdir+'TIR5.tif', sepoutdir+'TIR6.tif', sepoutdir+'TIR7.tif'], options=vrt_options)
                gdal.Translate(outputfilename+'.tif', outputfilename+'.vrt')
                gdal.Translate(outputfilename+'_template.tif', outputfilename+'.vrt')
                
                all_grids = gdal.Open(outputfilename+'.tif').ReadAsArray()
                surfrad = all_grids
                emisf_tes = zeros_like(surfrad)
                
                #------------------------------------------------------------#
                ########################T-E SEPARATION########################
                #------------------------------------------------------------#
                ecw_tir = [7.3, 8.4, 9.6, 10.3, 11.2, 12.3, 13.3]
                tes_bands=[1,3,4]
                ecw_tes = [ecw_tir[i] for i in tes_bands]
                
                #set up variables for TES algorithm: skyr_tes,skyr_nontes, surfrad_tes, and
                #surfrad_nontes
                [surfrad_tes, surfrad_nontes] = init_TES(surfrad)
                [emisf, emisf_nontes, Ts, QAmap, wave_tes, wave_nontes]=TES_for_ABI(surfrad_tes, surfrad_nontes, ecw_tes, ecw_tir)
                #combine output data 
                emisf_tes1 = np.float32(emisf_nontes[0,:,:])
                emisf_tes[0,:,:] = emisf_nontes[0,:,:]
                emisf_tes2 = np.float32(emisf[0,:,:])
                emisf_tes[1,:,:] = emisf[0,:,:]
                emisf_tes3 = np.float32(emisf_nontes[1,:,:])
                emisf_tes[2,:,:] = emisf_nontes[1,:,:]
                emisf_tes4 = np.float32(emisf[1,:,:])
                emisf_tes[3,:,:] = emisf[1,:,:]
                emisf_tes5 = np.float32(emisf[2,:,:])
                emisf_tes[4,:,:] = emisf[2,:,:]
                emisf_tes6 = np.float32(emisf_nontes[2,:,:])
                emisf_tes[5,:,:] = emisf_nontes[2,:,:]
                emisf_tes7 = np.float32(emisf_nontes[3,:,:])
                emisf_tes[6,:,:] = emisf_nontes[3,:,:]
                Ts1 = np.float32(Ts)
                
                #make Emi output directory if not already exist
                if not os.path.exists(l2_dir+'TES/Emi/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'):
                    os.makedirs(l2_dir+'TES/Emi/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/')
                #emi filename
                outputfilenameEmi = l2_dir+'TES/Emi/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/TES_Emi_'+Start[0:13]+'_'+End[0:13]
                #open tempate
                EmiTemplate = gdal.Open(outputfilename+'_template.tif')
                #get geotiff driver
                EmiDriver = gdal.GetDriverByName('GTiff')
                #EmiDriver = EmiTemplate.GetDriver()
                #Create new raster
                EmiRaster = EmiDriver.Create(outputfilenameEmi+'.tif', xsize=EmiTemplate.RasterXSize, ysize=EmiTemplate.RasterYSize, bands=7, eType=gdal.GDT_Float32)
                # Set metadata
                EmiRaster.SetGeoTransform(EmiTemplate.GetGeoTransform())
                EmiRaster.SetProjection(EmiTemplate.GetProjection())
                ##
                #read raster band as array
                EmiBand1 = EmiRaster.GetRasterBand(1).ReadAsArray()
                EmiBand1[:,:] = emisf_tes1
                EmiRaster.GetRasterBand(1).WriteArray(EmiBand1)
                EmiBand2 = EmiRaster.GetRasterBand(2).ReadAsArray()
                EmiBand2[:,:] = emisf_tes2
                EmiRaster.GetRasterBand(2).WriteArray(EmiBand2)
                EmiBand3 = EmiRaster.GetRasterBand(3).ReadAsArray()
                EmiBand3[:,:] = emisf_tes3
                EmiRaster.GetRasterBand(3).WriteArray(EmiBand3)
                EmiBand4 = EmiRaster.GetRasterBand(4).ReadAsArray()
                EmiBand4[:,:] = emisf_tes4
                EmiRaster.GetRasterBand(4).WriteArray(EmiBand4)
                EmiBand5 = EmiRaster.GetRasterBand(5).ReadAsArray()
                EmiBand5[:,:] = emisf_tes5
                EmiRaster.GetRasterBand(5).WriteArray(EmiBand5)                
                EmiBand6 = EmiRaster.GetRasterBand(6).ReadAsArray()
                EmiBand6[:,:] = emisf_tes6
                EmiRaster.GetRasterBand(6).WriteArray(EmiBand6)
                EmiBand7 = EmiRaster.GetRasterBand(7).ReadAsArray()
                EmiBand7[:,:] = emisf_tes7
                EmiRaster.GetRasterBand(7).WriteArray(EmiBand7)
                #close data file
                EmiRaster = None
                EmiTemplate = None
                
                
                #make Emi output directory if not already exist
                if not os.path.exists(l2_dir+'TES/Temp/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/'):
                    os.makedirs(l2_dir+'TES/Temp/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/')
                #temp filename
                outputfilenameTemp = l2_dir+'TES/Temp/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour+'/TES_Temp_'+Start[0:13]+'_'+End[0:13]
                #open tempate
                TempTemplate = gdal.Open(outputfilename+'_template.tif')
                #get geotiff driver
                TempDriver = gdal.GetDriverByName('GTiff')
                #Create new raster
                TempRaster = TempDriver.Create(outputfilenameTemp+'.tif', xsize=TempTemplate.RasterXSize, ysize=TempTemplate.RasterYSize, bands=1, eType=gdal.GDT_Float32)
                # Set metadata
                TempRaster.SetGeoTransform(TempTemplate.GetGeoTransform())
                TempRaster.SetProjection(TempTemplate.GetProjection())
                #read raster band as array
                TempBand = TempRaster.GetRasterBand(1).ReadAsArray()
                TempBand[:,:] = Ts1
                TempRaster.GetRasterBand(1).WriteArray(TempBand)
                # # Close datasets
                TempRaster=None
                TempTemplate=None
                #remove template file
                if os.path.isfile(outputfilename+'_template.tif'):
                    os.remove(outputfilename+'_template.tif')
                #------------------------------------------------------------#                

#----------------------------------------------------------------------------#
        # Report progress 

        working_volcano_out = [working_volcano_name,working_coordinates,working_newextent_range,str(i_dir_year),s_i_dir_day,s_i_dir_hour]

        progresslines=[]
        with open(logdir+'UWI_GEO_data_progress_log.csv', 'r') as readFile:    
            reader = csv.reader(readFile)
            for row in reader:
                progresslines.append(row)
                if row[0] == working_volcano_out[0]:
                    progresslines.remove(row)
    
        with open(logdir+'UWI_GEO_data_progress_log.csv', 'w', newline='') as writeFile:
            writer = csv.writer(writeFile)
            writer.writerows(progresslines)
            writer.writerow(working_volcano_out)
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
## END INDIVIDUAL VOLCANO PROCESSING LOOP ## 
#----------------------------------------------------------------------------#   
    #remove all L0 radaince data from previous hour
    shutil.rmtree(l1_dir+'Rad/'+str(i_dir_year)+'/'+s_i_dir_day+'/'+s_i_dir_hour)
    #--------------------------------------------------------------------#
    #move time forward one hour
    i_dir_hour = i_dir_hour + 1
    #move forward one day if required
    if (jdcal.is_leap(i_dir_year)!=False) and (i_dir_hour > 23) and (i_dir_day > 365):
        i_dir_day = 1
        i_dir_hour = 0
        i_dir_year = i_dir_year + 1
        print("Done day "+str(i_dir_day-1)+", "+str(i_dir_year-1))
    if (jdcal.is_leap(i_dir_year)==False) and (i_dir_hour > 23) and (i_dir_day > 364):
        i_dir_day = 1
        i_dir_hour = 0
        i_dir_year = i_dir_year + 1
        print("Done day "+str(i_dir_day-1)+", "+str(i_dir_year-1))
    if i_dir_hour > 23:
        i_dir_day = i_dir_day + 1
        i_dir_hour = 0
        print("Done day "+str(i_dir_day-1)+", "+str(i_dir_year))
    #--------------------------------------------------------------------#
    print('Done:  Day '+str(i_dir_day)+'/'+str(i_dir_year)+' '+str(i_dir_hour)+'hrs')
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
## END DOWNLAOD LOOP ## 
#----------------------------------------------------------------------------#