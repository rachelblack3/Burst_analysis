# numerical essentials 
import numpy as np
import pandas as pd

# for plotting
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import plot_script as pls

# for cdf reading
from spacepy import toolbox
#spacepy.toolbox.update(leapsecs=True)
from spacepy import pycdf
import h5py

# other numerical tools
import os
import glob
from datetime import datetime,date,timedelta

# for FFT
import fft_box as ft
import get_date as gd

# my own misc. functions
import funcs as funcs

# Tell me the date range - start date and end date. Function will also return the number days this range spans
start_date,end_date,no_days = funcs.what_dates()

""" 
The files being used for the FFTs are either the 6s burst waveform samples, or the 0.468s burst waveforms sent dwon every
15 minutes to check the onboard FFT algorithm. Choose file type.
"""

file_type = 'Burst6'

""" 
Now choose which FFT style you would like:
recreate_survey
or
imporve_res
"""

fft_type = input("Please give FFT type; rs = recreate survey FFT, and itr = improve temporal resoloution: ")


wfr_folder ='/data/spacecast/wave_database_v2/RBSP-A/L2'
burst6 = 'rbsp-a_WFR-waveform-continuous-burst_emfisis-L2_'

epoch_omni,AE=funcs.omni_stats(start_date,end_date)


for single_day in (start_date + timedelta(n) for n in range(no_days)):
    
    # String version of dates for filename  
    date_string,year,month,day =funcs.get_date_string(single_day)

    # Full filepath
    wfr_burst_path = os.path.join(wfr_folder, year, month, day, burst6 + date_string + "*_v*.cdf")

    # Accessing the magnetometer data

    mag_file=pycdf.CDF(funcs.magnetic_data(single_day,year,month,day))   

    # Accessing survey file 

    survey = pycdf.CDF(ft.survey_data(single_day,year,month,day)) 

    # files are all CDFs for a given day
    cdfs = glob.glob(wfr_burst_path)
    no_cdf = len(cdfs)

    # Define empty list for storing number of records in each CDF
    no_rec = []

    # empty list for holding gyrofrequencies

    g1=[]
    ghalf=[]
    glow=[]
    data_set = []
    freq_set=[]
    burst_times =[]

    b_box = []
    b_mean =[]
    s = []
    b_all = []

    b_dist=[]        # burst distribution for 12 0.468s samples
    high_b = []      # high resoloution burst distribution
    color_plots = [] # saving output for each of the 12 0.468s FFTs - so should get a coloutplot of dimension (12,Number_of_frequencies)



    if (no_cdf>0):

        # j counts every example for each day to append to plot names
        j=0

        for cdf in cdfs:
            j = j+1
            if j>1:
                break
            # Open burst file
            burst = pycdf.CDF(cdf)

            # Count how many records on each CDF
            no_rec = len(burst['Epoch'])

        
            B_cal = burst['BCalibrationCoef']       # B callibration values (for each of 6500 frequencies, regulary spaced)
            df_cal = burst['CalFrequencies'][0]     # Frequency step size required for burst callibaration
            f_max = burst['CalFrequencies'][-1]     # Max frequency (where the receivers 'cuttoff')

            # Now churning through each record in a given CDF
            for i in range(no_rec):
                #if i>1:
                    #break
                burst_times.append(burst['Epoch'][i])

                if i>20:
                    """ 
                    Bu Bv Bw samples from burst file 
                    """
                    Bu_sample = burst['BuSamples'][i]
                    Bv_sample = burst['BvSamples'][i]
                    Bw_sample = burst['BwSamples'][i]


                    # This code creates two new series derived from 'wn_series'; one is a "rolling" mean, the other is the "rolling" 
                    # standard deviation. Each of these simply applies a moving window to our white noise and computes the mean and standard
                    # deviation along the moving window.
                    Bu_series = pd.Series(Bw_sample)
                    roll_mean = Bu_series.rolling(window=16384).mean()
                    roll_sd = Bu_series.rolling(window=16384).std()

                    roll_mean_short = Bu_series.rolling(window=1024).mean()
                    roll_sd_short = Bu_series.rolling(window=1024).std()

                    T = 1.0/35000.0  
                    n_s = len(roll_mean_short)
                    n_l = len(roll_mean)
                    
                    raw_times_long = np.linspace(0.468/2.,5.968-0.468/2.,n_l)
                    raw_times_short = np.linspace(0.015/2.,5.968-0.015/2.,n_s)
                    #plt.plot(raw_times_long,roll_mean,color='blue',label='Mean (0.468s windows)')
                
                    #plt.plot(raw_times_long,roll_sd,color='red',label='SD (0.468s windows)')
                    #plt.plot(raw_times_short,roll_mean_short,color='blue',label='Mean (0.015s windows)',linestyle='dashed')
                   # plt.plot(raw_times_short,roll_sd_short,color='red',label='SD (0.015s windows)',linestyle='dashed')

        

                    dt = np.linspace(T,16384*T,16384)
                    dt_short = np.linspace(T,1024*T,1024)

                    v_lines = np.linspace(0.015/2.,31*0.015-0.015/2.,31)
                    v_lines2 = np.linspace(0.015,31*0.015,31)
                    v_lines3 = np.linspace(0.015,31*0.015,31)
                    #plt.plot(dt,Bu_sample[:16384])
                    plt.plot(dt,Bu_sample[:16384],label='0.468s window')
                    plt.vlines(v_lines,min(Bu_sample[:16384]),np.max(Bu_sample[:16384]),color='black',linestyles='dashed')
                    plt.vlines(v_lines2,min(Bu_sample[:16384]),np.max(Bu_sample[:16384]),color='red',linestyles='dashed')
                    
                    plt.xlabel('Time (s)')
                    plt.ylabel('Bu Field')
                # plt.vlines(1024*T,min(Bu_sample),max(Bu_sample),color='red',linestyle='dashed')
                    #plt.vlines(16384*T,min(Bu_sample),max(Bu_sample),color='green',linestyle='dashed')
                    plt.legend()
                    plt.title('0.015s windows in 0ne 0.468s sample')
                    plt.savefig('eg_windows')
                    break