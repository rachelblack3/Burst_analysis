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
from funcs import FFTMagnitude, DataFiles, AccessSurveyAttributes

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

epoch_omni,AE=funcs.omni_stats(start_date,end_date)


for single_day in (start_date + timedelta(n) for n in range(no_days)):
    
    day_files = DataFiles(single_day)

    # String version of dates for filename  
    date_string,year,month,day =funcs.get_date_string(single_day)

    # Full filepath
    wfr_burst_path = os.path.join(wfr_folder, year, month, day, burst6 + date_string + "*_v*.cdf")

    # Getting the magnetometer data
 
    mag_file = pycdf.CDF(day_files.magnetic_data)

    # Getting the LANL data

    lanl_file = h5py.File(day_files.lanl_data)
 
    # Getting survey file and accessing survey frequencies

    #survey = pycdf.CDF(ft.survey_data(single_day,year,month,day)) 
    survey_file = pycdf.CDF(day_files.survey_data)
    survey_freq = AccessSurveyAttributes(survey_file).frequency

    # files are all CDFs for a given day
    #cdfs = glob.glob(wfr_burst_path)
    cdfs = day_files.burst_paths
    no_cdf = len(cdfs)

    # Define empty list for storing number of records in each CDF
    no_rec = []

    # empty list for holding gyrofrequencies

    g1=[]                                           # 1*gyrofreqency
    ghalf=[]                                        # 0.5*gyrofrequeuncy
    glow=[]                                         # 0.05*gyrofrequency

    # empty lists for holding the small-window FFT results 
    data_set = []                                   # B field fft data
    freq_set=[]                                     # Frequencies corresponding to FFT data 
    high_b = []                                     # high resoloution frequency-integrated burst distribution 

    burst_times =[]                                 # datetime list of burst samples for each day
    surv_set = []                                   # survey magnetic field spectral density for each day

    b_box = []                                      # frequency integrated power over first 0.468s window spectrun
    s = []                                          # frequency intergrated power over survey spectrum 
    b_all = []

    b_dist=[]                                       # burst distribution for 12 0.468s samples
    color_plots = []                                # saving output for each of the 12 0.468s FFTs - so should get a coloutplot of dimension (12,Number_of_frequencies)

    av_specs = []                                   # magnetic spectral density averaged over all small windows
    av_specs_win=[]                                 # magnetic spectral density averaged over first window

    if (no_cdf>0):

        # j counts every example for each day to append to plot names
        j=0

        for cdf in cdfs:
            j = j+1

            # Open burst file
            burst = pycdf.CDF(cdf)

            # files for passing to functions more easily
            file_access = {
                "survey_file": survey_file,
                "burst_file": burst,
                "mag_file": mag_file}

            # Count how many records on each CDF
            no_rec = len(burst['Epoch'])
            burst_epoch = funcs.epoch_convert(burst['Epoch'])
        
            B_cal = burst['BCalibrationCoef']       # B callibration values (for each of 6500 frequencies, regulary spaced)
            df_cal = burst['CalFrequencies'][0]     # Frequency step size required for burst callibaration
            f_max = burst['CalFrequencies'][-1]     # Max frequency (where the receivers 'cuttoff')

            # this is just for the callibartion coefficients
            burst_params = {
                "B_cal": B_cal,
                "df_cal": df_cal,
                "f_max": f_max}

            time_range = [time(hour = 6,minute = 0,second=0),time(hour = 6,minute = 10,second=0)]
            # Now churning through each record in a given CDF
            for i in range(no_rec):
                if (time_range[0]< burst_epoch[i].time() <time_range[1]):
                    #break
                #if 19<i<22:
                    print(i+1, " bursts have now been identified")

                # date parameters for passing to functions more easily
                    date_params = {
                        "year": year,
                        "month": month,
                        "day": day,
                        "single_day": single_day,
                        "burst_time": burst['Epoch'][i]}
                    
                    # Save the burst time to list for plotting
                    burst_times.append(burst['Epoch'][i])


                    """ 
                    Bu Bv Bw waveforms from burst file 
                    """
                    Bu_sample = burst['BuSamples'][i]
                    Bv_sample = burst['BvSamples'][i]
                    Bw_sample = burst['BwSamples'][i]
                    
                    Bsamples = {"Bu": Bu_sample,
                               "Bv": Bv_sample,
                               "Bw": Bw_sample}

                    ''' Doing FFTs of the 0.468s samples from the burst sample
                        The ACTUAL survey is taken from the first 0.468s of the burst
                        Doing this for multiple other random sections of the burst 
                    '''   

                    rebinned_468, mag_468, params_468 = funcs.process_Kletzing_windows(Bsamples,
                            burst_params,
                            date_params,
                            file_access)
                    


                    ''' Doing FFTs of entire 6 second data ste, with overlapping windows of size 1024
                    '''
                    
                    mag_030,params_030 = funcs.process_small_windows(Bu_sample,Bv_sample,Bw_sample,
                            burst_params)  
                            
                    slider = int(1024/4)
                    
                    mag_030,params_030 = funcs.process_sliding_windows(Bu_sample,Bv_sample,Bw_sample,slider,burst_params)

                    ''' Comparing the first 0.468s window to the average of the first 0.468s of 0.030s windows
                    '''              
                    
                    av_comparison= np.zeros(params_030["n_f"])

                    for m in range(params_030["n_f"]):
                        av_comparison[m] = np.mean(mag_030[:61,m])
                    
                    rebinned_comp = funcs.rebin_burst(survey_file,av_comparison,params_030["freq"],single_day,year,month,day)
                    

                    """ 
                    Find survey for equivalent time
                    """
                    Btotal = funcs.process_survey(file_access)
                    surv_set.append(Btotal)

                    ''' Find gyrfofrequencies for plotting '''

                    # First find survey time for rec

                    time_s, surv_index= funcs.surv_burst(burst['Epoch'][i],AccessSurveyAttributes(survey_file).epoch)

                    # Then add calculated fractions of the gyrofrequency to their respective lists
                    g1.append(funcs.calc_gyro(mag_file,time_s)[0])
                    ghalf.append(funcs.calc_gyro(mag_file,time_s)[1])
                    glow.append(funcs.calc_gyro(mag_file,time_s)[2])
                    

                    data_set.append(mag_030)
                    freq_set.append(params_030["freq"])

                    """" Do the survey integration """

                    survey_int = funcs.integrate_in_rebin(Btotal[surv_index,:],survey_freq)
                            
                    s.append(survey_int)    

                    """" Do the burst integration """

                    burst_int = funcs.integrate_in_rebin(rebinned_468,survey_freq)
                            
                    b_box.append(burst_int)    

                    
                    """ Do the Kletzing FFT for first 12 0.468s samples """

                    box_dist,color_samples,reb_sancheck = funcs.box_dist(Bu_sample,Bv_sample,Bw_sample,B_cal,
                                                survey_file,survey_freq,
                                                year,month,day,
                                                single_day)
                    
                    b_dist.append(box_dist)
                    color_plots.append(color_samples)
                    
                    
                    print('The burst box is',burst_int,'and the survey is',survey_int, 'and sanity check burst is',box_dist[0])


                    """" Do the high burst integration """

                    intergral_statistics = funcs.integrate_in_small_windows(mag_030,params_030)

                    """ Averaging the 0.03s windows """

                    av_rebinned, av_rebinned_win = funcs.spectrum_averaging(mag_030,params_030,file_access,date_params)
                    av_specs.append(av_rebinned)

                    av_specs_win.append(av_rebinned_win)

                    high_b.append(intergral_statistics["frequency integrated power"])
                    
                    print('The standard deviation of the burst sample distribution is',np.std(box_dist))



    epoch_survey = AccessSurveyAttributes(survey_file).epoch_convert()


    plotting_dict = {"Frequencies for sliding windows": freq_set,
                     "PSD for sliding windows":data_set,
                     "0.5 gyrofrequency": ghalf,
                     "0.05 gyrofrequency": glow,
                     "Gyrofrequency": g1,
                     "Survey frequencies": survey_freq,
                     "Survey epoch": epoch_survey,
                     "survey PSD": Btotal,
                     "Timestamps on bursts": burst_times,
                     "LANL_dict": lanl_file,
                     "AE": AE,
                     "OMNI epoch": epoch_omni,
                     "Year":year,
                     "Month": month,
                     "Day": day,
                     "Survey integral":s,
                     "First 0.468s FFT integral": b_box,
                     "0.468s window powers": b_dist,
                     "Sliding window integraed power":high_b,
                     "PSD of 0.468s windows": color_plots,
                     "Frequencies for 0.468s windows": params_468["freq"],
                     "Power averaged over small windows": av_specs,
                     "Power averaged over first 0.468s": av_specs_win,
                     "Survey PSD for each day": surv_set}
            
    """ 
    Plot results
    """
    summaries = pls.summary_plot(freq_set,data_set,
                 ghalf,glow,g1,
                 survey_freq,epoch_survey,Btotal,
                 burst_times,
                 lanl_file,AE,epoch_omni,
                 year,month,day,
                 survey_int,
                 burst_int,
                 b_dist,
                 high_b,
                 color_plots,params_468["freq"],
                 av_specs,av_specs_win,surv_set)