#----------------------------------------------------------------------------#
# Created on Tue Sep  7 11:36:49 2021
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
import numpy as np # version 1.19.2
import os # version 
import csv # version 
import shutil # version
from osgeo import gdal # version 2.3.2
#import pandas # version 1.1.3
import jdcal
import datetime
from getListOfFiles import getListOfFiles
import matplotlib
import matplotlib.pyplot as plt
from baseline_als_optimized import baseline_als_optimized
#from statsmodels import robust
matplotlib.use('Agg') # use a non-interactive backend
#matplotlib.use('Qt5Agg')

#----------------------------------------------------------------------------#
################################Start of Script###############################
#----------------------------------------------------------------------------#
#constants
SBcons = 5.670374419e-08
axismonth = ['01/21','04/21','07/21','10/21', \
'01/22','04/22','07/22','10/22', \
'01/23','04/23','07/23','10/23', \
'01/24','04/24','07/24','10/24', \
'01/25','04/25','07/25','10/25', \
'01/26','04/26','07/26','10/26', \
'01/27','04/27','07/27','10/27', \
'01/28','04/28','07/28','10/28', \
'01/29','04/29','07/29','10/29']

axistimestamp = [datetime.datetime.strptime(date, '%m/%y').timestamp() for date in axismonth]

#----------------------------------------------------------------------------#
## Import logs and email data
logdir = './Logs/'
indir = './Data/'

with open(logdir+'UWI_GEO_active_volcanoes.csv', newline='') as f:
    reader = csv.reader(f)
    eruption_list = list(reader)
activevolcanoes = eruption_list

#%%
for e in range(len(activevolcanoes)):
    
    #empty arrays
    data_timestamp=[]
    SummitTempMean = []
    SummitTempMean4032 = []
    SummitTempMean1008 = []
    SummitTempMean144 = []
    allmeanbaselineabovesurftemp_70 = []
    allmaxbaselineabovesurftemp_70 = []
    alltotalRadPower = []
    all_count_meanabovesurftemp = []
    allTempData = []
    allback93 = []
    SummitTempMax=[]
    
    activevolcano1 = activevolcanoes[e]
    print('staring: '+str(activevolcano1))    
    #input variables
    volcano = activevolcano1[0]
    L3summitoutdir = indir+volcano+'/L3/Summit/'
    mediadir = indir[:-5]+'Media/'+volcano+'/'
    coordinates = activevolcano1[1] 
    lat_coord = float( coordinates.split(', ')[0] )
    long_coord = float( coordinates.split(', ')[1] )
    extent_range = 0.08
    geoextent = [long_coord-extent_range, lat_coord+extent_range, long_coord+extent_range, lat_coord-extent_range] #min long, min lat, max long, max lat 
    
    allfiles = getListOfFiles(indir+volcano+'/L2/TES/Temp/')
    allfiles = [sub.replace('\\', '/') for sub in allfiles] 
    tiffiles = list(filter(lambda x: '.tif' in x, allfiles))
    L2Tempfiles = list(filter(lambda x: 'TES_Temp' in x, tiffiles))
    L2Tempfiles = sorted(L2Tempfiles)

#%%
    print('complied lists')
    nL2Tempfiles = L2Tempfiles
    for i in range(len(nL2Tempfiles)):
        filename = L2Tempfiles[i]
        filename = filename[filename.find('TES_Temp')+9:filename.find(".tif")]
        
        #determine timestamp
        Data_Jdfmt = ("%Y%j%H%M%S")
        JD_date = datetime.datetime.strptime(str(filename.split('_')[0]), Data_Jdfmt)
        file_year = filename[0:4]
        file_day = filename[4:7]
        file_hour = filename[7:9]
        
        if not os.path.exists(L3summitoutdir+'Rad/Rad_TIR_'+filename+'.tif'):        
            #rad data resample
            L1bRadfile1 = indir+volcano+'/L1b/Rad/'+file_year+'/'+file_day+'/'+file_hour+'/Rad_TIR_'+filename+'.tif'
            if os.path.exists(L1bRadfile1):
                if not os.path.exists(L3summitoutdir+'Rad/'):
                    os.makedirs(L3summitoutdir+'Rad/')
                gdal.Translate(L3summitoutdir+'Rad/Rad_TIR_'+filename+'.tif', L1bRadfile1, projWin=geoextent)
            else:
                continue
            
        if not os.path.exists(L3summitoutdir+'TES/Emi/Emi_'+filename+'.tif'):        
            #emi data resample
            L2Emifile1 = indir+volcano+'/L2/TES/Emi/'+file_year+'/'+file_day+'/'+file_hour+'/TES_Emi_'+filename+'.tif'
            if os.path.exists(L2Emifile1):
                if not os.path.exists(L3summitoutdir+'TES/Emi/'):
                    os.makedirs(L3summitoutdir+'TES/Emi/')
                gdal.Translate(L3summitoutdir+'TES/Emi/Emi_'+filename+'.tif', L2Emifile1, projWin=geoextent)
            else:
                continue
            
        if not os.path.exists(L3summitoutdir+'TES/Temp/Temp_'+filename+'.tif'):        
            #emi data resample        
            L2Tempfile1 = indir+volcano+'/L2/TES/Temp/'+file_year+'/'+file_day+'/'+file_hour+'/TES_Temp_'+filename+'.tif'
            if os.path.exists(L2Emifile1):
                if not os.path.exists(L3summitoutdir+'TES/Temp/'):
                    os.makedirs(L3summitoutdir+'TES/Temp/')
                gdal.Translate(L3summitoutdir+'TES/Temp/Temp_'+filename+'.tif', L2Tempfile1, projWin=geoextent)
            else:
                continue
            
        RadData = gdal.Open(L3summitoutdir+'Rad/Rad_TIR_'+filename+'.tif').ReadAsArray()
        EmiData = gdal.Open(L3summitoutdir+'TES/Emi/Emi_'+filename+'.tif').ReadAsArray()
        AvgEmiData = np.mean(EmiData[1:], axis=0)
        TempData = gdal.Open(L3summitoutdir+'TES/Temp/Temp_'+filename+'.tif').ReadAsArray()
        
        TempDataP69 = np.nanpercentile(TempData, 69)
        TempDataP69 = np.nan_to_num(np.asarray(TempDataP69), nan=273.15)
        SummitTempMean.append(TempDataP69)
        TempDataMax = np.nanmax(TempData)
        TempDataMax = np.nan_to_num(np.asarray(TempDataMax), nan=273.15)
        SummitTempMax.append(TempDataMax)
        
        data_timestamp.append(JD_date.timestamp())

    baseline1 = abs(baseline_als_optimized(SummitTempMean, 13104, 0.93))
   
    for i in range(len(nL2Tempfiles)):
        filename = L2Tempfiles[i]
        filename = filename[filename.find('TES_Temp')+9:filename.find(".tif")]
        file_year = filename[0:4]
        file_day = filename[4:7]
        file_hour = filename[7:9]
        
        if not os.path.exists(indir+volcano+'/L1b/Rad/'+file_year+'/'+file_day+'/'+file_hour+'/Rad_TIR_'+filename+'.tif'):
            continue
        if not os.path.exists(indir+volcano+'/L2/TES/Emi/'+file_year+'/'+file_day+'/'+file_hour+'/TES_Emi_'+filename+'.tif'):
            continue
        if not os.path.exists(indir+volcano+'/L2/TES/Temp/'+file_year+'/'+file_day+'/'+file_hour+'/TES_Temp_'+filename+'.tif'):
            continue
            
        TempData = gdal.Open(L3summitoutdir+'TES/Temp/Temp_'+filename+'.tif').ReadAsArray()
        TempData4032 = TempData - baseline1[i]
        
        TempData4032P69 = np.nanpercentile(TempData4032, 69)
        TempData4032P69 = np.nan_to_num(np.asarray(TempData4032P69), nan=0.0)
        if TempData4032P69 < 0:
            TempData4032P69 = 0.0
        SummitTempMean4032.append(TempData4032P69)

    baseline2 = abs(baseline_als_optimized(SummitTempMean4032, 1008, 0.93))
    
    for i in range(len(nL2Tempfiles)):
        filename = L2Tempfiles[i]
        filename = filename[filename.find('TES_Temp')+9:filename.find(".tif")]
        file_year = filename[0:4]
        file_day = filename[4:7]
        file_hour = filename[7:9]
        
        if not os.path.exists(indir+volcano+'/L1b/Rad/'+file_year+'/'+file_day+'/'+file_hour+'/Rad_TIR_'+filename+'.tif'):
            continue
        if not os.path.exists(indir+volcano+'/L2/TES/Emi/'+file_year+'/'+file_day+'/'+file_hour+'/TES_Emi_'+filename+'.tif'):
            continue
        if not os.path.exists(indir+volcano+'/L2/TES/Temp/'+file_year+'/'+file_day+'/'+file_hour+'/TES_Temp_'+filename+'.tif'):
            continue
        
        TempData = gdal.Open(L3summitoutdir+'TES/Temp/Temp_'+filename+'.tif').ReadAsArray()
        TempData1008 = (TempData - baseline1[i]) - baseline2[i]
        
        TempData1008P69 = np.nanpercentile(TempData1008, 69)
        TempData1008P69 = np.nan_to_num(np.asarray(TempData1008P69), nan=0.0)        
        if TempData1008P69 < 0:
            TempData1008P69 = 0.0
        SummitTempMean1008.append(TempData1008P69)

    baseline3 = abs(baseline_als_optimized(SummitTempMean1008, 36, 0.93))    
    
    for i in range(len(nL2Tempfiles)):
        filename = L2Tempfiles[i]
        filename = filename[filename.find('TES_Temp')+9:filename.find(".tif")]
        file_year = filename[0:4]
        file_day = filename[4:7]
        file_hour = filename[7:9]
        
        if not os.path.exists(indir+volcano+'/L1b/Rad/'+file_year+'/'+file_day+'/'+file_hour+'/Rad_TIR_'+filename+'.tif'):
            continue
        if not os.path.exists(indir+volcano+'/L2/TES/Emi/'+file_year+'/'+file_day+'/'+file_hour+'/TES_Emi_'+filename+'.tif'):
            continue
        if not os.path.exists(indir+volcano+'/L2/TES/Temp/'+file_year+'/'+file_day+'/'+file_hour+'/TES_Temp_'+filename+'.tif'):
            continue
        
        EmiData = gdal.Open(L3summitoutdir+'TES/Emi/Emi_'+filename+'.tif').ReadAsArray()
        AvgEmiData = np.mean(EmiData[1:], axis=0)
        TempData = gdal.Open(L3summitoutdir+'TES/Temp/Temp_'+filename+'.tif').ReadAsArray()       
        allTempData.append(TempData)

        baselineabovesurftemp = ((TempData - baseline1[i]) - baseline2[i]) - baseline3[i]
        copybaselineabovesurftemp = baselineabovesurftemp.copy()
        nancopybaselineabovesurftemp = copybaselineabovesurftemp[copybaselineabovesurftemp > 0]
        copybaselineabovesurftemp = np.nan_to_num(np.asarray(copybaselineabovesurftemp), nan=0.0)
        copybaselineabovesurftemp = copybaselineabovesurftemp[copybaselineabovesurftemp > 0]
        if copybaselineabovesurftemp.size > 0:
            back93 = abs(np.nanpercentile(nancopybaselineabovesurftemp,93))
            allback93.append(back93)
            baselineabovesurftemp_70 = baselineabovesurftemp - np.nanpercentile(nancopybaselineabovesurftemp,93)
            baselineabovesurftemp_70[baselineabovesurftemp_70 < 0] = 0.0
            baselineabovesurftemp_70 = np.nan_to_num(np.asarray(baselineabovesurftemp_70), nan=0.0)
            count_meanabovesurftemp = np.nansum(baselineabovesurftemp_70 > 0) 
            all_count_meanabovesurftemp.append(count_meanabovesurftemp)
        else:
            count_meanabovesurftemp = np.array([])
            all_count_meanabovesurftemp.append(0.0)      
            allback93.append(0.0)
            
        if count_meanabovesurftemp > 0:
            totalbackground = baseline1[i] + baseline2[i]  + back93 + baseline3[i]
            RadPower = AvgEmiData * SBcons * ((TempData**4)-(totalbackground**4)) * 4e6 * 1e-6#4e6 is the pixel size in m2 and in MW
            RadPower[RadPower < 0] = np.nan
            totalRadPower = np.nansum(RadPower)
            totalRadPower = round(totalRadPower,5)
            alltotalRadPower.append(totalRadPower)
    
            baselineabovesurftemp_70[baselineabovesurftemp_70 == 0.0] = np.nan
            meanbaselineabovesurftemp_70 = np.nanmean(baselineabovesurftemp_70)
            allmeanbaselineabovesurftemp_70.append(meanbaselineabovesurftemp_70)
            maxbaselineabovesurftemp_70 = np.nanmax(baselineabovesurftemp_70)
            allmaxbaselineabovesurftemp_70.append(maxbaselineabovesurftemp_70)
            
            baselineabovesurftemp_70[np.isnan(baselineabovesurftemp_70)] = 0
        
        else:
            alltotalRadPower.append(0.0)
            allmeanbaselineabovesurftemp_70.append(0.0)
            allmaxbaselineabovesurftemp_70.append(0.0)
   
    #temperature plot
    plt.figure(figsize=[12,6])
    plt.plot(data_timestamp,allmeanbaselineabovesurftemp_70, label='Median', color='black')
    plt.xticks(axistimestamp, axismonth, fontsize=14, rotation=0)
    plt.tick_params(axis ='both', direction='out', length=8) 
    plt.ylim(0, np.ceil(np.nanmax(allmeanbaselineabovesurftemp_70)))
    tempymax = np.ceil(np.nanmax(allmeanbaselineabovesurftemp_70))
    tempyticks = np.array([0, (tempymax/5), ((tempymax/5)*2), ((tempymax/5)*3), ((tempymax/5)*4), tempymax])
    plt.yticks(tempyticks,fontsize=14)
    plt.xlim(data_timestamp[0]-100000, data_timestamp[-1]+100000)
    plt.xlabel('Date', fontsize=16)
    plt.ylabel('Mean Temperature\nAbove Background (Kelvin)', fontsize=16)
    #build directory
    if not os.path.exists(mediadir):
        os.makedirs(mediadir)
    plt.savefig(mediadir+volcano+'MeanTemperatureSeries.png', dpi=200, bbox_inches='tight', pad_inches=0)
    plt.close()  

    #temperature plot
    plt.figure(figsize=[12,6])
    plt.plot(data_timestamp,allmaxbaselineabovesurftemp_70, label='Median', color='black')
    plt.xticks(axistimestamp, axismonth, fontsize=14, rotation=0)
    plt.tick_params(axis ='both', direction='out', length=8) 
    plt.ylim(0, np.ceil(np.nanmax(allmaxbaselineabovesurftemp_70)))
    tempymax = np.ceil(np.nanmax(allmaxbaselineabovesurftemp_70))
    tempyticks = np.array([0, (tempymax/5), ((tempymax/5)*2), ((tempymax/5)*3), ((tempymax/5)*4), tempymax])
    plt.yticks(tempyticks,fontsize=14)
    plt.xlim(data_timestamp[0]-100000, data_timestamp[-1]+100000)
    plt.xlabel('Date', fontsize=16)
    plt.ylabel('Max Temperature\nAbove Background (Kelvin)', fontsize=16)
    #build directory
    if not os.path.exists(mediadir):
        os.makedirs(mediadir)
    plt.savefig(mediadir+volcano+'MaxTemperatureSeries.png', dpi=200, bbox_inches='tight', pad_inches=0)
    plt.close() 
      
    #HF plot
    plt.figure(figsize=[12,6])
    plt.plot(data_timestamp,alltotalRadPower, color='black')
    plt.xticks(axistimestamp, axismonth, fontsize=14, rotation=0)
    plt.tick_params(axis ='both', direction='out', length=8) 
    plt.ylim(0, np.ceil(np.nanmax(alltotalRadPower)))
    HFymax = np.ceil(np.nanmax(alltotalRadPower))
    HFyticks = np.array([0, (HFymax/5), ((HFymax/5)*2), ((HFymax/5)*3), ((HFymax/5)*4), HFymax])
    plt.yticks(HFyticks,fontsize=14)
    plt.xlim(data_timestamp[0]-100000, data_timestamp[-1]+100000)
    plt.xlabel('Date', fontsize=16)
    plt.ylabel('Total Heat Flux (MW)', fontsize=16)
    #build directory
    if not os.path.exists(mediadir):
        os.makedirs(mediadir)
    plt.savefig(mediadir+volcano+'HFSeries.png', dpi=200, bbox_inches='tight', pad_inches=0)
    plt.close()   
    
    header = ['Time Stamp (datetime.timestamp)', 'Mean Above Background Temperature (Kelvin)', 'Max Above Background Temperature (Kelvin)', 'Total Radiative Heat Flux (MW)', 'SummitTempMean', 'SummitTempMax', 'baseline1', 'baseline2', 'baseline3', 'allback93']
    with open(mediadir+volcano+'_data.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for w in range(len(data_timestamp)):
            writer.writerow([data_timestamp[w],allmeanbaselineabovesurftemp_70[w],allmaxbaselineabovesurftemp_70[w],alltotalRadPower[w]])    
    
    end_datetime = datetime.datetime.fromtimestamp(data_timestamp[-1])
    s_i_dir_year = str(end_datetime.timetuple().tm_year)
    s_i_dir_day = str(end_datetime.timetuple().tm_yday)
    s_i_dir_hour = str(end_datetime.timetuple().tm_hour)
    volcano_out = [volcano,coordinates,geoextent,s_i_dir_year,s_i_dir_day,s_i_dir_hour]

    progresslines=[]
    with open(logdir+'UWI_GEO_analysis_log.csv', 'r') as readFile:    
        reader = csv.reader(readFile)
        for row in reader:
            progresslines.append(row)
            if row[0] == volcano_out[0]:
                progresslines.remove(row)

    with open(logdir+'UWI_GEO_analysis_log.csv', 'w', newline='') as writeFile:
        writer = csv.writer(writeFile)
        writer.writerows(progresslines)
        writer.writerow(volcano_out)

#%%