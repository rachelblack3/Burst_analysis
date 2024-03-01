import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import fft_box as ft


''' Creating waveforms to test the FFT routines - and to investigate the difference between small and large windows '''

# Just a regular sine wave 

def sin(t):

    y = np.sin(t)

    return y

# A sine wave that increases with power at somepoint during the duration

def burst_sin(t):

    y = 10*np.sin(t)

    return y

# A rising tone
# Make the frequency changes happen after one full period of the previous frequency

def rising_tone():
    
    f_s = 1/35000    
    freq = []

    f_min = 10
    n_f = 11
    f_max = n_f*f_min

    # change frequency n_f times within duration, and then create time array based on fnal duration
    
    T = 0.
    for i in range(1,n_f):

        freq_step = i*f_min                  # frequency element of rising tone
        period_step = 1./freq_step           # for the given frequency, how long to complete a full period
        T = T + period_step                  # adding up each period to get full rising tone duration

        rng = int((period_step/f_s))         # number of time elements for this frequency 
        
        freq.extend([freq_step]*(rng))       # extend the list to have this frequency for the desired number of steps
        

    N = len(freq)
    
    t = np.linspace(0.,T, N)                 # creating the time array fo plotting
    
    argument = (freq * t )                     # argument inside sin function

    data = np.sin(argument)                    # artificual rising waveform data


    return data,t,freq,N


def fft_short(udata,vdata,wdata):
    from scipy.fft import fft, fftfreq
    import numpy as np

    N = len(udata)                                               # Number of sample points

    T = 1.0/35000.0                                         # Time between samples
    T_window = N*T
    
    w = 1. #np.hanning(N)                                       # hanning window
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
    freq = fftfreq(N, T)
    df = freq[1]-freq[0]
    
   
    """ 
    define complex array for storing FFTs
    """
    fft_arr = np.zeros((3,len(freq)),dtype=complex)

    """ 
    Multiply each element of fft array by correspondong complex coefficient
    """
    fft_arr[0,:] = bu_fft
    fft_arr[1,:] = bv_fft
    fft_arr[2,:] = bw_fft
    #bmag_fft = bmag_fft[1:n_f+1]*Bcal_c

    """ 
    take absolute values and square to get power density
    """
    total = abs(fft_arr) * abs(fft_arr)

    """ 
    divide array by correction for spectral power lost by using a hanning window
    """
    total= total/wms
    
    """ 
    find B field magnitude
    """
    
    mag =(total[0,:]+total[1,:]+total[2,:])*2*T_window

    return mag,freq, bu_fft


r_tone, t, chosen_freq, N = rising_tone()

FFT_mag, FFT_freq, bu_fft = fft_short(r_tone,r_tone,r_tone)

bu_ifft = np.fft.ifft(bu_fft,len(bu_fft))

fig,axs = plt.subplots(5,1)

ax1,ax2,ax3,ax4,ax5 = axs
ax1.plot(t, r_tone, label="sin(freq * x)")
ax1.set_xlabel("Time")
ax1.set_ylabel("Amplitude")
ax1.legend()

ax2.plot(t,chosen_freq,label="Artificially chosen frequency variation")
ax2.set_xlabel("Time")
ax2.set_ylabel("Frequency")
ax2.legend()

ax3.plot(FFT_freq,bu_fft,label = "FFT spectrum")
ax3.set_xlim(-200,200)
ax3.set_xlabel("Frequency")
ax3.set_ylabel("Power spectral density")
ax3.legend()

ax4.plot(t,bu_ifft)

ax5.plot(t,bu_ifft/r_tone,label='Ratio')
ax5.legend()
ax5.set_xlabel('Time')
ax5.set_ylabel('Ratio of signal with inverse FFT')

plt.gcf().set_size_inches((12, 30))
plt.show()

    
