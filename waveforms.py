import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import fft_box as ft

# for plotting
import matplotlib.colors as mcolors
import matplotlib.dates as mdates


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

    amp = 10                                # multiply rising tone by 10 = 10x background values
    
    f_s = 1/35000    
    freq = []

    f_min = 100
    n_f = 8
    f_max = n_f*f_min

    # change frequency n_f times within duration, and then create time array based on fnal duration
    
    T = 0.
    for i in range(1,n_f):
        
        freq_step = (i/i)*2*f_min                  # frequency element of rising tone
        omega = freq_step*2*np.pi
        period_step = 1./freq_step                 # for the given frequency, how long to complete a full period
        print(freq_step,period_step)
        T = T + i*50*period_step                    # adding up each period to get full rising tone duration

        rng = int((i*50*period_step/f_s))           # number of time elements for this frequency 
        print(rng)
        freq.extend([omega]*(rng))       # extend the list to have this frequency for the desired number of steps
        f_min = freq_step   

    print("The length of rising tone data is ",len(freq))
    N = len(freq)
    
    t = np.linspace(0.,T, N)                 # creating the time array fo plotting
    
    argument = (freq * t )                   # argument inside sin function

    data = np.sin(argument)                  # artificual rising waveform data

    freq = [x/(2*np.pi) for x in freq]

    return data,t,freq,N,T


def background(f_s):
    
    duration = 6.
    n_t = int(duration/f_s)

    freq_b = 20                                                # 'background' broadband frequency

    full_t = np.linspace(0,6.,n_t)
    full_t = np.array(full_t)

    freq = np.full_like(full_t, freq_b)                        # making a frequency array filled with background value
                                                               # same length as time array
    omega = 2*np.pi*freq

    argument = (omega * full_t )
    
    data = np.sin(argument)

    return full_t,data, freq, n_t

def fft_full(udata):
    from scipy.fft import fft, fftfreq
    import numpy as np

    N = len(udata)                                               # Number of sample points

    T = 1.0/35000.0                                              # Time between samples
    T_window = N*T
    
    w = np.hanning(N)                                       # hanning window
    wms = np.mean(w**2)                                     # window spectral power correction


    """ 
    Do FFTs - remember to normalise by N and multiply by hanning window
    """
    
    bu_fft = fft(udata*w*(1/N))


    """ 
    Frequencies from FFT
    """
    freq = fftfreq(N, T)[:N//2]
    df = freq[1]-freq[0]
    
    """ 
    Multiply each element of fft array by correspondong complex coefficient
    """
    bu_fft = fft(udata*w*(1/N))[0:N//2]
    #bmag_fft = bmag_fft[1:n_f+1]*Bcal_c

    """ 
    take absolute values and square to get power density
    """
    total = abs(bu_fft) * abs(bu_fft)

    """ 
    divide array by correction for spectral power lost by using a hanning window
    """
    total= total/wms
    
    """ 
    find B field magnitude
    """
    
    mag = total*2*T_window
    

    return mag,freq, bu_fft


def fft_windows(udata,n_box):

    from scipy.fft import fft, fftfreq,rfftfreq
    import numpy as np

   
    N = len(udata)                                                       # Number of sample points
    box_size = int(n_box)                                                # Number of samples in each box 
    N_box = int(N/box_size)                                              # Number of boxes
    half_window = int(box_size/2)                                        # Size of half-window
    
    T = 1.0/35000.0                                                      # Time between samples

    """ 
    Frequencies from FFT boxes
    """
    
    freq = rfftfreq(box_size, T)[:box_size//2]
    n_f = len(freq)
    df = freq[1] - freq[0]

                                       
    """ 
    Do FFTs - remember to normalise by N and multiply by hanning window and complex correction
    """

    w = np.hanning(box_size)                                            # hanning window
    wms = np.mean(w**2)                                                 # hanning window correction

    n_bins = (2*N_box)-1
    fft_box = np.zeros((n_bins,n_f),dtype = complex)
    print('Number of bins',np.shape(fft_box))
    # Doing the first window
    # Remember to normalise the box size
    

    fft_box[0,:] = (fft(w*(1/box_size)*udata[0:box_size]))[0:box_size//2]
   
    
    # Doing remainder of overlapping windows, and multiplying by complex callibration
    for i in range(1,n_bins):                                               
                                                    
        index_start = i*half_window-1
        index_end = (i*half_window)+box_size-1

                                                
        fft_box[i,:] = (fft(w*(1/box_size)*udata[index_start:index_end]))[0:box_size//2]
        

    fft_av = fft_box
    for i in range(1,n_bins-1):

        fft_av[i,:] = (fft_box[i,:]+fft_box[i+1,:])/2.
       

    T_window = box_size*T

    print('the frequency bands for little bins are',df)


    """ 
    take absolute values and square to get power density
    """
    total = abs(fft_av) * abs(fft_av)

    """ 
    divide array by correction for spectral power lost by using a hanning window
    """
    total= total/wms
    
    """ 
    find B field magnitude
    """
    
    mag = total*2*T_window

    return mag,freq,fft_av,n_bins

def fft_over_full_interval(t,r_tone,
                           chosen_freq,
                           FFT_freq,FFT_mag):
                           #bu_ifft):
    
    """ Creating the plots for the FFT over the full interval """

    fig,axs = plt.subplots(3,1)

    ax1,ax2,ax3 = axs
    ax1.plot(t, r_tone, label="sin(freq * x)")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Amplitude")
    ax1.legend()

    ax2.plot(t,chosen_freq,label="Artificially chosen frequency variation")
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Frequency")
    ax2.legend()

    ax3.plot(FFT_freq,FFT_mag,label = "FFT spectrum")
    ax3.set_xlim(0,50)
    ax3.set_xlabel("Frequency")
    ax3.set_ylabel("Power spectral density")
    ax3.legend()

    #ax4.plot(t,bu_ifft)

    #ax5.plot(t,bu_ifft/r_tone,label='Ratio')
    #ax5.legend()
    #ax5.set_xlabel('Time')
    #ax5.set_ylabel('Ratio of signal with inverse FFT')

    plt.gcf().set_size_inches((12, 18))
    plt.show()

    return

def fft_over_windows(t,r_tone,
                    chosen_freq,
                    FFT_freq,FFT_mag,
                    window_times):
    
    """ Creating the plots for the FFT over the full interval """

    # create figure
    fig, axs = plt.subplots(3,2,gridspec_kw={"width_ratios":(1,0.05)})

    # Name plot axis
    ax1,ax2,ax3= axs[:, 0]

    # Name colourbar axis (the third/fourth one will be made invisible)
    cax1,cax2,cax3 = axs[:,1]

    
    print(np.shape(t),np.shape(r_tone))
    ax1.plot(t, r_tone, label="sin(freq * x)")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Amplitude")
    ax1.legend()
    cax1.axis('off')

    ax2.plot(t,chosen_freq,label="Artificially chosen frequency variation")
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Frequency")
    ax2.legend()
    cax2.axis('off')

    colorbar_norm = mcolors.LogNorm(vmin=10**(-11), vmax=10**(-6))

    window_plot = ax3.pcolormesh(window_times, FFT_freq, np.transpose(FFT_mag), norm=colorbar_norm, cmap="jet",label = "FFT spectrum")
        
    
    ax3.set_ylabel('Frequency (Hz)',fontsize=8)
    ax3.set_xlabel('Time',fontsize=8)
    plt.colorbar(window_plot,label='power',cax=cax3)
    #ax3.set_ylim(0,15000)
    
    ax3.legend()


    plt.gcf().set_size_inches((12, 18))
    plt.savefig('windowed.png')
    plt.show()

    return


def comparison_plot(FFT_freq_nw,FFT_mag_nw,
                    FFT_mag_w,ratio):
    
    """ Creating the plots for the FFT over the full interval """

    fig,axs = plt.subplots(2,1)

    ax1,ax2 = axs
    

    ax1.plot(FFT_freq_nw,FFT_mag_nw,label = "Full interval FFT")
    ax1.plot(FFT_freq_nw,FFT_mag_w,label = "FFT with overlapping windows")
    ax1.set_yscale('log')
    ax1.set_xlabel("Frequency")
    ax1.set_ylabel("Power spectral density")
    ax1.legend()

    
    ax2.plot(FFT_freq_nw,ratio,label = r'$PSD\_{Total\ FFT}/PSD\_{Window\ averaged}$')
    #ax2.set_yscale('log')
    ax2.set_xlabel("Frequency")
    ax2.set_ylabel("Power spectral density")
    ax2.legend()



    plt.gcf().set_size_inches((12, 30))
    plt.savefig('Comaprison_FFT.png')
    plt.show()


    return



r_tone, t, chosen_freq, N, period = rising_tone()

full_t, background_data, background_freq, n_timesteps = background(f_s = 1./35000.)
print("the lengh of the entrie 6s data is",len(background_data))
background_data[int(len(background_data)/2) - int(len(r_tone)/2) :int(len(background_data)/2) + int(len(r_tone)/2)] = r_tone
background_freq[int(len(background_data)/2) - int(len(r_tone)/2) :int(len(background_data)/2) + int(len(r_tone)/2)] = chosen_freq


rising_mag, rising_freq, rising_bu_fft = fft_full(r_tone)

comp_mag, comp_freq, comp_bu_fft = fft_full(background_data)
#comp_bu_ifft = np.fft.ifft(comp_bu_fft,len(comp_bu_fft))


#bu_ifft = np.fft.ifft(rising_bu_fft,len(rising_bu_fft))


just_rising = {"time array": t,
                    "waveform": r_tone,
                    "artifical frequencies": chosen_freq,
                    "FFT_frequencies": rising_freq,
                    "FFT_psd": rising_mag}
                    #"Inverse FFT": bu_ifft}

composite_wave = {"time array": full_t,
                    "waveform": background_data,
                    "artifical frequencies": background_freq,
                    "FFT_frequencies": comp_freq,
                    "FFT_psd": comp_mag}
                    #"Inverse FFT": comp_bu_ifft}

print("the composite background and rising tone shape is", np.shape(background_data))
print("Number of frequencies in full interval FFT is", len(composite_wave["FFT_frequencies"]))


plots_full_interval = fft_over_full_interval(just_rising["time array"],just_rising["waveform"],
                           just_rising["artifical frequencies"],
                           just_rising["FFT_frequencies"],just_rising["FFT_psd"])
                           #just_rising["Inverse FFT"])

plots_with_background = fft_over_full_interval(composite_wave["time array"],composite_wave["waveform"],
                           composite_wave["artifical frequencies"],
                           composite_wave["FFT_frequencies"],composite_wave["FFT_psd"])
                           #composite_wave["Inverse FFT"])



FFT_mag, FFT_freq, bu_fft, n_time_bins = fft_windows(r_tone,350)
comp_mag, comp_freq, comp_bu_fft, comp_time_bins = fft_windows(background_data,3500)


averaging_timesteps = np.linspace(0,period, n_time_bins)
comp_duration = 6.
comp_timesteps = np.linspace(0.,comp_duration, comp_time_bins)

windows_fft = {"time array": averaging_timesteps,
                    "rising tone waveform": r_tone,
                    "artifical frequencies": chosen_freq,
                    "FFT_frequencies": FFT_freq,
                    "FFT_psd": FFT_mag}

comp_windows_fft = {"time array": comp_timesteps,
                    "rising tone waveform": background_data,
                    "artifical frequencies": background_freq,
                    "FFT_frequencies": comp_freq,
                    "FFT_psd": comp_mag}

plots_windows = fft_over_windows(just_rising["time array"],just_rising["waveform"],
                           just_rising["artifical frequencies"],
                           windows_fft["FFT_frequencies"],windows_fft["FFT_psd"],
                           windows_fft["time array"])

plots_windows = fft_over_windows(composite_wave["time array"],composite_wave["waveform"],
                           composite_wave["artifical frequencies"],
                           comp_windows_fft["FFT_frequencies"],comp_windows_fft["FFT_psd"],
                           comp_windows_fft["time array"])



# Average the PSD over all of the windows to compare with FFT over full interval
def average_windows(mag,frequency):
    n_f = len(frequency)
    averaged_spec = np.zeros_like(frequency)  

    for m in range(n_f):
        averaged_spec[m] = np.mean(mag[:,m])
  
    return averaged_spec


def rebin(bin_edges,mag,frequency_original):
    """ 
    Doing the rebinning
    """
    import pandas as pd

    # Create dataframe
    rebin_dat=pd.DataFrame()

    rebin_dat['Data'] = mag
    
    # Create and save frequencies to one column
    rebin_dat['Frequency']= frequency_original
    
    """
    pd.cut() bins all frequencies according to defined semi_log bins
    groupby() groups all data in frame by these bines
    then select the DATA and take the MEAN in each bin
    to_numpy saves this to an array
    """
    
    rebinned=rebin_dat.groupby(pd.cut(rebin_dat.Frequency,bins=bin_edges)).Data.mean().to_numpy()
    
    return rebinned

# Average the overallping window FFT over full duration
averaged_windows = average_windows(comp_windows_fft["FFT_psd"],comp_windows_fft["FFT_frequencies"])

# Define the frequency bins for the window FFT, so that we can rebin the higher frequency resoloution long FFT into them

bin_edges = list(comp_windows_fft["FFT_frequencies"])
bin_edges.append(bin_edges[-1]+166.6)

# Rebin full FFT into window frequency bins 
rebinned_fullFFT = rebin(bin_edges,composite_wave["FFT_psd"],composite_wave["FFT_frequencies"])

ratio_of_psds = rebinned_fullFFT/averaged_windows

# Plot the rebinned full FFT next to the window FFT
do_comparison = comparison_plot(comp_windows_fft["FFT_frequencies"],rebinned_fullFFT,averaged_windows,ratio_of_psds)


""" Now, we want to check how the ratio between true frequency change with the resoloution from the FFT windows """

# In other words, we want to change the number of windows, and see how that effects the output we get out

n_steps = 4

fig,axs = plt.subplots(2,1)

ax1,ax2 = axs

for j in range(1,n_steps+1):

    # Do the windowing with different number of windows

    n_elements = len(background_data)/(10**j)

    n_windows = 10**j

    print(n_windows,'Done')

    comp_mag, comp_freq, comp_bu_fft, comp_time_bins = fft_windows(background_data,n_elements)

    comp_timesteps = np.linspace(0.,comp_duration, comp_time_bins) 


    comp_windows_fft = {"time array": comp_timesteps,
                        "rising tone waveform": background_data,
                        "artifical frequencies": background_freq,
                        "FFT_frequencies": comp_freq,
                        "FFT_psd": comp_mag}
    
    # Average the overallping window FFT over full duration
    averaged_windows = average_windows(comp_windows_fft["FFT_psd"],comp_windows_fft["FFT_frequencies"])    

    bin_edges = list(comp_windows_fft["FFT_frequencies"])
    bin_edges.append(bin_edges[-1]+166.6)

    # Rebin full FFT into window frequency bins 
    rebinned_fullFFT = rebin(bin_edges,composite_wave["FFT_psd"],composite_wave["FFT_frequencies"])

    ratio_of_psds = rebinned_fullFFT/averaged_windows


    label = "Number of windows =" + str(n_windows)

    ax1.plot(comp_windows_fft['FFT_frequencies'],rebinned_fullFFT,label = "Full FFT - "+ str(n_windows) + " windows")
    ax1.plot(comp_windows_fft["FFT_frequencies"],averaged_windows,label = " Av. FFT - "+ str(n_windows) + " windows")
    ax1.set_yscale('log')
    ax1.set_xlabel("Frequency")
    ax1.set_ylabel("Power spectral density")

    ax2.plot(comp_windows_fft["FFT_frequencies"],ratio_of_psds,label = "Ratio -"+ str(n_windows) + " windows")
    #ax2.set_yscale('log')
    ax2.set_xlabel("Frequency")
    ax2.set_ylabel("Ratio")

# Put a legend to the right of the current axis
ax1.legend(loc='center left', bbox_to_anchor=(1., 0.5),fontsize = 6)


ax2.legend(loc='center left',bbox_to_anchor=(1., 0.5),fontsize = 6)
plt.tight_layout()
plt.savefig("effect_of_nw.png")
plt.show()