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
def rising_tone():
    T = 1.
    f_s = 1/35000
    N = int(T/f_s)

    t = np.linspace(0,T,N)
    freq = np.zeros_like(t)

    f_min = 10
    n_f = 20
    f_max = n_f*f_min

    # change frequency 20 times within 6s duration
    for i in range(n_f):
        for j in range(int(i*N/n_f),int((i+1)*N/n_f)):
            freq[j] = i*f_min


    x_plot = (freq * t )  # arhument inside sin function

    data = np.sin(x_plot)

    return data,t,freq


def fft_short(udata,vdata,wdata):
    from scipy.fft import fft, fftfreq
    import numpy as np

    N = len(udata)                                               # Number of sample points

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


r_tone, t, chosen_freq = rising_tone()

FFT_mag, FFT_freq, bu_fft = fft_short(r_tone,r_tone,r_tone)

bu_ifft = np.fft.ifft(bu_fft,len(bu_fft))

fig,axs = plt.subplots(4,1)

ax1,ax2,ax3,ax4 = axs
ax1.plot(t, r_tone, label="sin(freq * x)")
ax1.set_xlabel("Time")
ax1.set_ylabel("Amplitude")
ax1.legend()

ax2.plot(t,chosen_freq,label="Artificially chosen frequency variation")
ax2.set_xlabel("Time")
ax2.set_ylabel("Frequency")
ax2.legend()

ax3.plot(FFT_freq,FFT_mag,label = "FFT spectrum")
ax3.set_xlim(-200,200)
ax3.set_xlabel("Frequency")
ax3.set_ylabel("Power spectral density")
ax3.legend()

ax4.plot(t,bu_ifft)

plt.gcf().set_size_inches((12, 24))
plt.show()

    
