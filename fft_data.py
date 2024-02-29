import glob
import os
import spacepy
from spacepy import pycdf

def fft_dat(udata,vdata,wdata,Bcal):
    from scipy.fft import fft, fftfreq
    import numpy as np

   
    N = 16384                                               # Number of sample points
    n_f = 5600                                              # Number of required frequencies
    T = 1.0/35000.0                                         # Time between samples
    
    
    w = np.hanning(N)                                       # hanning window
    wms = np.mean(w**2)                                     # window spectral power correction


    """ 
    Do FFTs - remember to normalise by N and multiply by hanning window
    """
    
    bu_fft = fft(udata*w*(1/N))
    bv_fft = fft(vdata*w*(1/N))
    bw_fft = fft(wdata*w*(1/N))

    """ 
    Frequencies from FFT
    """
    freq = fftfreq(N, T)[1:n_f+1]
    
 

    """ 
    Putting the calibration coeffcients into complex form
    """
   
    Bcal_c=[]
    for i in range(n_f):
        Bcal_c.append(complex(Bcal[i,0],Bcal[i,1]))
    
    """ 
    define complex array for storing FFTs
    """
    fft_arr = np.zeros((3,n_f),dtype=complex)

    """ 
    Multiply each element of fft array by correspondong complex coefficient
    """
    fft_arr[0,:] = bu_fft[1:n_f+1]*Bcal_c
    fft_arr[1,:] = bv_fft[1:n_f+1]*Bcal_c
    fft_arr[2,:] = bw_fft[1:n_f+1]*Bcal_c
    #bmag_fft = bmag_fft[1:n_f+1]*Bcal_c

    return fft_arr,freq, wms


""" survey_data: find path to survey data on that date """

def survey_data(start_date,year,month,day):
    date_string= str(start_date.strftime("%Y%m%d"))
    survey_folder ='/data/spacecast/satellite/RBSP/emfisis/data/RBSP-A/L2'
    survey_file= 'rbsp-a_WFR-spectral-matrix-diagonal_emfisis-L2_'
    survey_path = os.path.join(survey_folder, year, month, day,survey_file + date_string + "_v*.cdf")
    
    # find the latest version
    survey_path = glob.glob(survey_path)[-1]
    
    return survey_path


#def read_bins():
    # Read file with semi-logarithmic frequency bin boundaries
    #widths = []
    #f = open("widths.txt", "r")
    #y = f.read().split(',')
    #for string in y:
     #   widths.append(float(string))
   
    #widths = survey['WFR_bandwidth'].tolist()

    #return widths



def bin_edges(widths,start_date,year,month,day):
    """ 
    setting the bin edges 
    min_bin: lower edge of first bin 
    """
    survey = pycdf.CDF(survey_data(start_date,year,month,day))
    survey_freq = survey['WFR_frequencies'][0]

    min_bin = survey_freq[0]-(widths[0]/2)

    freq_bin = []
    freq_bin.append(min_bin)

    # Starting from the minimum bin, add all widths to get lower and upper bands on all semi_logarithmic bins

    for i in range(0,65):
        
        freq_bin.append(freq_bin[i]+widths[i])
        min_bin=freq_bin[i]
    
    
    return freq_bin, survey_freq
    

    

