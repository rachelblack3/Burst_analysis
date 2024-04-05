# numerical essentials 
import numpy as np
import pandas as pd

# for plotting
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import plot_script as pls
#import old_plotting as pls
# for cdf reading
from spacepy import toolbox
#spacepy.toolbox.update(leapsecs=True)
from spacepy import pycdf
import h5py

# other numerical tools
import os
import glob
from datetime import datetime,date,timedelta,time

# for FFT
import fft_box as ft
import get_date as gd

# my own misc. functions
import funcs as funcs

# import class for holding all the FFT data 
from funcs import DataFiles, AccessSurveyAttrs

""" Defining all global variables 

Burst-related parameters
NOTE: 'burst' is the local burst file name

'Bu_sample'     - Magnetic field waveform measured in the u direction; nT
'Bv_sample'     - Magnetic field waveform measured in the v direction; nT
'Bw_sample'     - Magnetic field waveform measured in the w direction; nT
'B_cal'         - Complex magnetic spectral density callibration values (for each of 6500 frequencies, regulary spaced) ; nT/Hz
'df_cal'        - Frequency step size required for burst callibaration; Hz
'f_max'         - Maximum frequency, above which receivers 'cuttoff'; Hz
'burst_times'   - List of all burst samples on a given day; datetime elements 

"""

''' Tell me the date range - start date and end date. Function will also return the number days this range spans '''
start_date,end_date,no_days = funcs.what_dates()

""" 
The files being used for the FFTs are either the 6s burst waveform samples, or the 0.468s burst waveforms sent dwon every
15 minutes to check the onboard FFT algorithm. Choose file type.
"""

wfr_folder ='/data/spacecast/wave_database_v2/RBSP-A/L2'
burst6 = 'rbsp-a_WFR-waveform-continuous-burst_emfisis-L2_'

epoch_omni,AE= funcs.omni_dataset(start_date,end_date).omni_stats
# Plot 3/8 in summary plot - geomagnetic indicies

plot3 = {"AE": AE,
         "Epoch": epoch_omni}


for single_day in (start_date + timedelta(n) for n in range(no_days)):
    
    day_files = DataFiles(single_day)

    # String version of dates for filename  
    date_string,year,month,day =funcs.get_date_string(single_day)

    # Full filepath
    wfr_burst_path = os.path.join(wfr_folder, year, month, day, burst6 + date_string + "*_v*.cdf")

    # Getting the magnetometer data
 
    mag_file = pycdf.CDF(day_files.magnetic_data)
    mag_data = funcs.AccessL3Attrs(mag_file)
    # Getting the density data

    density_file = pycdf.CDF(day_files.l4_data)
    density_data = funcs.AccessL4Attrs(density_file)
    # Getting the LANL data

    lanl_file = h5py.File(day_files.lanl_data)
    lanl_data = funcs.AccessLANLAttrs(lanl_file)
    # Getting survey file and accessing survey frequencies, epoch and magnitude

    survey_file = pycdf.CDF(day_files.survey_data)
    survey_data = AccessSurveyAttrs(survey_file)

    survey_freq = survey_data.frequency
    survey_epoch = survey_data.epoch_convert()
    Btotal = survey_data.Bmagnitude

    # Getting gyrofrequencies and plasma frequency for the full day

    fpe = density_data.f_pe
    fpe_epoch = density_data.epoch_convert()

    fce, fce_05, fce_005 = mag_data.f_ce
    fce_epoch = mag_data.epoch_convert()

    # Getting LANL attributes
    Lstar = lanl_data.L_star
    MLT = lanl_data.MLT
    MLAT_N, MLAT_S = lanl_data.MLAT_N_S
    lanl_epoch = lanl_data.epoch

    # Plot 1/8 in summary plot - Survey PSD with trimmings

    plot1 = {"Bmagnitude": Btotal, "Frequency": survey_freq, "Epoch": survey_epoch,
             "fpe": fpe, "fpe_epoch": fpe_epoch,
             "fce": fce, "fce_05": fce_05, "fce_005": fce_005, "fce_epoch": fce_epoch}
    
    
    # Plot 2/8 in summary plot - spacecraft locations
    plot2 = {"MLT": MLT, "L*": Lstar, "MLAT North": MLAT_N, "MLAT South": MLAT_S,
             "Epoch": lanl_epoch}
    

    # add xlims to plot 3

    lower_xlim, upper_xlim = lanl_data.day_limits
    plot3["lower xlim"] =  lower_xlim
    plot3["upper xlim"] =  upper_xlim
    
    # get all burst CDFs for chosen day 
    cdfs = day_files.burst_paths
    no_cdf = len(cdfs)

    # Define empty list for storing number of records in each CDF
    no_rec = []

    if (no_cdf>0):

        # j counts every example for each day to append to plot names
        j=0

        for cdf in cdfs:
            j = j+1

            # Open burst file
            burst = pycdf.CDF(cdf)

            # Count how many records on each CDF
            no_rec = len(burst['Epoch'])
            burst_epoch = funcs.get_epoch(burst['Epoch'])
        
            B_cal = burst['BCalibrationCoef']       # B callibration values (for each of 6500 frequencies, regulary spaced)
            df_cal = burst['CalFrequencies'][0]     # Frequency step size required for burst callibaration
            f_max = burst['CalFrequencies'][-1]     # Max frequency (where the receivers 'cuttoff')

            # this is just for the callibartion coefficients
            burst_params = {
                "B_cal": B_cal,
                "df_cal": df_cal,
                "f_max": f_max}

            #time_range = [time(hour = 6,minute = 0,second=0),time(hour = 6,minute = 10,second=0)] I do want the abillity to put in sepcific events
            # Now churning through each record in a given CDF
            for i in range(no_rec):
                    
                if (i%10 == 0):
                    print(i, " bursts have now been identified")

                # date parameters for passing to functions more easily
                date_params = {
                    "year": year,
                    "month": month,
                    "day": day,
                    "single_day": single_day,
                    "burst_time": burst['Epoch'][i]}

                ''' Find gyrfofrequencies for plotting '''

                gyro1,gyrohalf,gyrolow = funcs.cross_dataset(survey_file, mag_file, date_params["burst_time"]).calc_gyro()

                fces = gyro1,gyrohalf,gyrolow

                """ 
                Bu Bv Bw waveforms from burst file 
                """
                
                Bsamples = {"Bu": burst['BuSamples'][i],
                            "Bv": burst['BvSamples'][i],
                            "Bw": burst['BwSamples'][i]}

                ''' Doing FFTs of the 0.468s samples from the burst sample
                    The ACTUAL survey is taken from the first 0.468s of the burst
                    Doing this for multiple other random sections of the burst 
                '''   
                slider = funcs.global_constants["Slider"]
                FFT_Kletzing = funcs.PerformFFT(Bsamples,burst_params,slider).process_Kletzing_windows()
                
                ''' Doing FFTs of entire 6 second data ste, with overlapping windows of size 1024
                '''
                FFT_sliding = funcs.PerformFFT(Bsamples,burst_params,slider).process_sliding_windows()  
                
                make_cdf = funcs.createPSDFiles(FFT_sliding,date_params,fces).save_FFT()
                make_cdf = funcs.createPSDFiles(FFT_Kletzing,date_params,fces).save_FFT()