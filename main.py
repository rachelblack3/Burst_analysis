import numpy as np
import pandas as pd

# for plotting
import matplotlib.pyplot as plt

# for cdf reading
from spacepy import pycdf

# other numerical tools
import numpy as np
import os
import glob
from datetime import datetime,date,timedelta

# for FFT
import fft_data as ft
import get_date as gd

# Float definitions
year = 2014
month = 10
day = 3 
start_date = date(year=year, month=month, day=day)                      # Date
full_dt = datetime(2014, 10, 3, 1, 48, 35, 989845)                      # Datetime

# String definitions
date_string,year,month,day =gd.get_date_string(start_date)     

""" 
This file is the 0.5s BURST sample they send down every 15 minutes to make sure FFT analysis is working correctly
This was 0.468s long in each case
"""

short_burst=pycdf.CDF('/users/rablack75/rbsp-a_WFR-waveform_emfisis-L2_20141003_v1.3.5.cdf')
#print(random_file.keys())
print(short_burst['fftSize'][0])


""" 
B callibration values (for each of 6500 frequencies, regulary spaced)
"""
B_cal = short_burst['BCalibrationCoef']

""" 
Bu Bv Bw samples
"""
Bu_sample = short_burst['BuSamples'][4]
Bv_sample = short_burst['BvSamples'][4]
Bw_sample = short_burst['BwSamples'][4]

""" 
calling FFT routine
"""
fft_arr, freq, wms = ft.fft_dat(Bu_sample,Bv_sample,Bw_sample,B_cal)
print(np.shape(fft_arr))
""" 
take absolute values and square to get power density
"""
total_m = abs(fft_arr) * abs(fft_arr)
total_m = total_m/wms
""" 
find B field magnitude
"""
mag = np.zeros(len(freq))

for n in range(len(freq)):

    mag[n]=(total_m[0,n]+total_m[1,n]+total_m[2,n])*2/2.136                     # Take mean             
    #mag[n] = total_m[0,n]+total_m[1,n]+total_m[2,n]

""" 
Read in semi-log bins and define bin edges
"""
#bins = ft.read_bins()
survey = pycdf.CDF(ft.survey_data(start_date,year,month,day))
bins = survey['WFR_bandwidth'][0]
freq_bin, rebinned_freq = ft.bin_edges(bins,start_date,year,month,day)

""" 
Doing the rebinning
"""
# Create dataframe
rebin_dat=pd.DataFrame()

rebin_dat['Data'] = mag

# Create and save frequencies to one column
rebin_dat['Frequency']= freq

"""
pd.cut() bins all frequencies according to defined semi_log bins
groupby() groups all data in frame by these bines
then select the DATA and take the MEAN in each bin
to_numpy saves this to an array
"""

rebinned=rebin_dat.groupby(pd.cut(rebin_dat.Frequency,bins=freq_bin)).Data.mean().to_numpy()
rebinned = rebinned

""" 
Find survey for equivalent time
"""
survey = pycdf.CDF(ft.survey_data(start_date,year,month,day))
survey_freq = survey['WFR_frequencies'][0]

Bu2 = survey['BuBu']
Bv2 = survey['BvBv']
Bw2 = survey['BwBw']

# Define empty list for total mag field array 

Btotal = np.zeros(survey['BuBu'].shape)
Btot= survey['WFR_bandwidth'][:]

# Create total mag B array

for i in range(0,np.shape(Btotal)[0]):
    Btotal[i,:] =Bu2[i,:]+ Bv2[i,:] + Bw2[i,:]



""" 
Plot results
"""

fig,axs = plt.subplots(4,1)

ax1,ax2,ax3,ax4 = axs

ax1.plot(rebinned_freq,rebinned, label = r'$FFT \ B_{mag}\ (65\ central\ frequencies)$')
ax1.plot(freq,mag, label = r'$FFT\ B_{mag}\ (6500\ central\ frequencies)$')
ax1.set_ylabel(r'$Power\ Spectral\ Density\ ({nT/Hz}^2)$')
ax1.legend()

ax2.plot(rebinned_freq,rebinned, label = r'$FFT\ B_{mag}\ (65\ central\ frequencies)$')
ax2.plot(survey_freq,Btotal[635,:],linestyle='dotted', label = r'$Survey\ B_{mag}\ (65\ central\ frequencies)$')
ax2.set_ylabel(r'$Power\ Spectral\ Density\ ({nT/Hz}^2)$')
ax2.legend()

ax3.plot(rebinned_freq,rebinned, label = r'$FFT\ B_{mag}\ (65\ central\ frequencies)$')
ax3.plot(survey_freq,Btotal[635,:],linestyle='dotted', label =r'$Survey\ B_{mag}\ (65\ central\ frequencies)$')
ax3.set_ylabel(r'$Power\ Spectral\ Density\ ({nT/Hz}^2)$')
ax3.set_yscale('log')
ax3.legend()

ax4.plot(survey_freq,Btotal[635,:]/rebinned,linestyle='dotted', label = r'Survey\ B_{mag}/Burst\ B_{mag}')
ax1.legend()
ax4.set_ylabel('Ratio')
plt.gcf().set_size_inches((6, 12))
plt.tight_layout()


file_path = '/data/spacecast/wip/nmer/rablack75/simple_FFT'+year
os.makedirs(file_path, exist_ok=True)
os.makedirs(file_path+'/'+month, exist_ok=True)
plt.savefig(file_path+'/'+month+'/'+str(date)+'attempt'+'.png')

plt.savefig('plot_'+datetime.strftime(survey['Epoch'][635],'%Y-%m-%d %H:%M:%S')+'.png')



