import glob
import os
import spacepy
from spacepy import pycdf

''' FFT for full 6 second plots, with overlapping windows '''

def fft_dat(udata,vdata,wdata,Bcal,df_cal, f_max):
    from scipy.fft import fft, fftfreq,rfftfreq
    import numpy as np

   
    N = len(udata)                                                       # Number of sample points
    box_size = 1024                                                      # Number of samples in each box 
    N_box = int(N/box_size)                                              # Number of boxes
    half_window = int(box_size/2)                                        # Size of half-window
    
    n_f = 5600                                                           # Number of required frequencies
    T = 1.0/35000.0                                                      # Time between samples

    """ 
    Frequencies from FFT boxes
    """
    n_f = half_window-1
    freq = rfftfreq(box_size, T)[1:n_f+1]
    
    df = freq[1] - freq[0]

    """ 
    Putting the calibration coeffcients into complex form
    """
    cal_step = int(df/df_cal)
    Bcal_c=[]
    
    for i in range(n_f):
        f = freq[i]
        if (f < f_max):
            Bcal_c.append(complex(Bcal[i*cal_step,0],Bcal[i*cal_step,1]))
        else:
            break

    n_f = len(Bcal_c)

    # Cutting off frequencies where the reciever stopped picking them up (fmax from callibration specifications)
    freq = freq[0:n_f]
                                       
    """ 
    Do FFTs - remember to normalise by N and multiply by hanning window and complex correction
    """

    w = np.hanning(box_size)                                            # hanning window
    wms = np.mean(w**2)                                                 # hanning window correction

    n_bins = (2*N_box)-1
    fft_box = np.zeros((n_bins,n_f,3),dtype = complex)

    # Doing the first window
    # Remember to normalise the box size

    #data[0:box_size] = linear_fit(udata[0:box_size],box_size)
    #vdata[0:box_size] = linear_fit(vdata[0:box_size],box_size)
    #wdata[0:box_size] = linear_fit(wdata[0:box_size],box_size)

    fft_box[0,:,0] = (fft(w*(1/box_size)*udata[0:box_size]))[1:n_f+1]*Bcal_c    
    fft_box[0,:,1] = (fft(w*(1/box_size)*vdata[0:box_size]))[1:n_f+1]*Bcal_c 
    fft_box[0,:,2] = (fft(w*(1/box_size)*wdata[0:box_size]))[1:n_f+1]*Bcal_c 
    
    # Doing remainder of overlapping windows, and multiplying by complex callibration
    for i in range(1,n_bins):                                               
                                                    
        index_start = i*half_window-1
        index_end = (i*half_window)+box_size-1

        #udata[index_start:index_end] = linear_fit(udata[index_start:index_end],box_size)
        #vdata[index_start:index_end] = linear_fit(vdata[index_start:index_end],box_size)
        #wdata[index_start:index_end] = linear_fit(wdata[index_start:index_end],box_size)
                                                
        fft_box[i,:,0] = (fft(w*(1/box_size)*udata[index_start:index_end]))[1:n_f+1]*Bcal_c    
        fft_box[i,:,1] = (fft(w*(1/box_size)*vdata[index_start:index_end]))[1:n_f+1]*Bcal_c 
        fft_box[i,:,2] = (fft(w*(1/box_size)*wdata[index_start:index_end]))[1:n_f+1]*Bcal_c 

    fft_av = fft_box
    for i in range(1,n_bins-1):

        fft_av[i,:,0] = (fft_box[i,:,0]+fft_box[i+1,:,0])/2.
        fft_av[i,:,1] = (fft_box[i,:,1]+fft_box[i+1,:,1])/2.
        fft_av[i,:,2] = (fft_box[i,:,2]+fft_box[i+1,:,2])/2.

    T_window = box_size*T

    return fft_av,freq, wms, n_bins, df,T_window


''' FFT for 0.468s second samples, in same style as RBSP onboard FFT process '''

def fft_short(udata,vdata,wdata,Bcal):
    from scipy.fft import fft, fftfreq
    import numpy as np

    big_N = 208896
    N = 16384                                               # Number of sample points

    N_box = big_N/N

    n_f = 5600                                              # Number of required frequencies
    T = 1.0/35000.0                                         # Time between samples
    T_window = N*T
    
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
    df = freq[1]-freq[0]

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

    return fft_arr,freq, wms,df,T_window



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
    
def linear_fit(data,boxsize):
    
    import numpy as np
    dt = 1./35000.
    x = np.linspace(0,boxsize*dt,boxsize)
    coef = np.polyfit(x,data,1)
    fit = np.poly1d(coef) 

    new_data = data - fit(x)

    return new_data

