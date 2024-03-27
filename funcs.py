import numpy as np
import os
import glob
from datetime import datetime,date,timedelta
import pandas as pd
from spacepy import pycdf
import fft_box as ft


class DataFiles:
    ''' class for accessing all of the data files that I need on a given day.  
    
            PARAMETERS:
            date: DateTime object 
           '''

    def __init__(self,date):
        self.date = date
        

    # finding survey data filepath for given day
    @property
    def survey_data(self):
    
        # String version of dates for filename  
        date_string,year,month,day =self.get_date_string()

        # root folder
        survey_folder ='/data/spacecast/satellite/RBSP/emfisis/data/RBSP-A/L2'
        # 'stem' name for burst CDFs
        survey_file= 'rbsp-a_WFR-spectral-matrix-diagonal_emfisis-L2_'

        survey_path = os.path.join(survey_folder, year, month, day,survey_file + date_string + "_v*.cdf")
        
        # find the latest version
        survey_path = glob.glob(survey_path)[-1]
        
        return survey_path


    # finding the magnetic data filepath for given day
    @property
    def magnetic_data(self):

        # String version of dates for filename  
        date_string,year,month,day =self.get_date_string()

        # root folder
        l3_folder ='/data/spacecast/satellite/RBSP/emfisis/data/RBSP-A/L3'
        # 'stem' name for burst CDFs
        l3_file= 'rbsp-a_magnetometer_1sec-geo_emfisis-L3_'

        l3_path = os.path.join(l3_folder, year, month, day,l3_file + date_string + "_v*.cdf")

        # find the latest version
        l3_path = glob.glob(l3_path)[-1]
        
        return l3_path

    # finding the LANL data filepath for given day 
    @property
    def lanl_data(self):
        # String version of dates for filename  
        date_string,year,month,day =self.get_date_string()

        # root folder
        lanl_folder ='/data/spacecast/satellite/RBSP/emfisis/data/RBSP-A/LANL/MagEphem'
        # 'stem' name for burst CDFs
        lanl_file= 'rbspa_def_MagEphem_TS04D_'

        lanl_path = os.path.join(lanl_folder, year, lanl_file + date_string + "_v*.h5")

        # find the latest version
        lanl_path = glob.glob(lanl_path)[-1]
        
        return lanl_path

    # Get L4 data filepath (aka density file) for given day
    @property
    def l4_data(self):

        # String version of dates for filename  
        date_string,year,month,day =self.get_date_string()

        # root folder
        l4_folder ='/data/spacecast/satellite/RBSP/emfisis/data/RBSP-A/L4'
        # 'stem' name for burst CDFs
        l4_file= 'rbsp-a_density_emfisis-L4_'

        l4_path = os.path.join(l4_folder, year, month, day,l4_file + date_string + "_v*.cdf")

        # find the latest version
        l4_path = glob.glob(l4_path)[-1]
        
        return l4_path
    
    # Get burst filepaths for all CDFs on given day
    @property
    def burst_paths(self):

        # root folder
        wfr_folder ='/data/spacecast/wave_database_v2/RBSP-A/L2'

        # 'stem' name for burst CDFs
        burst6 = 'rbsp-a_WFR-waveform-continuous-burst_emfisis-L2_'

        # String version of dates for filename  
        date_string,year,month,day =self.get_date_string()

        # Full filepath
        wfr_burst_path = os.path.join(wfr_folder, year, month, day, burst6 + date_string + "*_v*.cdf")

        # files are all CDFs for a given day
        cdf_files= glob.glob(wfr_burst_path)

        return cdf_files

    def get_date_string(self):
        ''' Method that rpovides date strings
        Outputs:
    
        date_string: string object of date
        year: string object of date year
        month: string object of date month
        day: string object of date day '''

        date_string= str(self.date.strftime("%Y%m%d"))

        if (self.date.day <10):
            day = "0"+str(self.date.day)
        else:
            day = str(self.date.day)


        if (self.date.month<10):
            month = "0"+str(self.date.month)
        else:
            month = str(self.date.month)
        

        year = str(self.date.year)
        
        return date_string,year,month,day

class AccessSurveyAttributes:
    ''' Class for finding and working on Survey data attributes
     
      PARAMETERS:
      survey_cdf: A CDF containing all survey data '''

    def __init__(self, survey_cdf):
        self.survey_cdf = survey_cdf

    @property
    def frequency(self):
        # Frequency in Hz
        frequency = self.survey_cdf['WFR_frequencies'][0]

        return frequency
    
    @property
    def epoch(self):
        # Epoch in DateTime format of: 
        epoch = self.survey_cdf['Epoch']

        return epoch
    
    def epoch_convert(self):
    # get epoch in correct format for comparisons etc.
        epoch = []

        # First - saving the epoch elements to a list as a string - acording to a particular desired format
        for i in range(len(self.survey_cdf['Epoch'])):
            epoch.append(datetime.strftime(self.survey_cdf['Epoch'],'%Y-%m-%d %H-%M-%S'))

        # Chaning back to a datetime object of same format    
        for i in range(len(epoch)):
            epoch[i] = datetime.strptime(epoch[i],'%Y-%m-%d %H-%M-%S')

        return(epoch)

# Class for the reuslts from windowed FFTs 

class FFTMagnitude:
    ''' class for storing instances of windowed FFT results, where window size and rolling window shift can be altered 
    
            PARAMETERS:
            magnitude: 2D array of FFT spectra, with shape (time_bins, frequency_bins) 
            frequency: central frequencies for each frequency band (Hz)
            n_window:  parameter deciding window size (window size = total number of samples/n_window)
            n_rolling: number of samples between rolling window '''

    def __init__(self,magnitude,frequency,n_window,n_rollling):
        self.magnitude=magnitude
        self.frequency=frequency
        self.n_window = n_window
        self.n_rolling = n_rollling

    def show_magnitude(self):
        print(f"The dimensions (time bins, frequency bins) of this FFT spectra is: {self.magnitude.shape}")

    def show_frequency(self):
        print(f"The number of frequencies is: {len(self.frequency)}")

    def __repr__(self):

        return f"magnitude shape = {self.magnitude.shape}, freq"

# Example showing how this works
test = FFTMagnitude(np.asarray([1,2]), [6],10,100)
test.show_frequency()
test.show_magnitude()
print(test.magnitude)

def what_dates():

    ''' Function that asks command prompt the required date range
        Outputs:

        start_date: Datetime object of start date
        end_date: Datetime object of end date
        no_days: Integer number of days
    '''

    # Ask for start date and make datetime object
    year,month,day =  [eval(i) for i in input("Please enter start year, month, day: ").split(",")]
    start_date = date(year=year, month=month, day=day)                                                       

    # Ask for end date and make datetime object
    year_end,month_end,day_end =  [eval(i) for i in input("Please enter end year, month, day: ").split(",")]
    end_date = date(year=year_end, month=month_end, day=day_end)

    # Find integer number of days in given date range
    no_days = end_date - start_date
    no_days = int(no_days.days)                     

    return (start_date,end_date,no_days)


def get_date_string(single_date):

    date_string= str(single_date.strftime("%Y%m%d"))

    if (single_date.day <10):
        day = "0"+str(single_date.day)
    else:
        day = str(single_date.day)


    if (single_date.month<10):
        month = "0"+str(single_date.month)
    else:
        month = str(single_date.month)
    

    year = str(single_date.year)
    
    return date_string,year,month,day


# Epoch conversions --> from 

def epoch_convert(epoch_raw):
    # get epoch in correct format for comparisons etc.
    epoch = []

    # First - saving the epoch elements to a list as a string - acording to a particular desired format
    for i in range(len(epoch_raw)):
        epoch.append(datetime.strftime(epoch_raw[i],'%Y-%m-%d %H-%M-%S'))

    # Chaning back to a datetime object of same format    
    for i in range(len(epoch)):
        epoch[i]=datetime.strptime(epoch[i],'%Y-%m-%d %H-%M-%S')

    return(epoch)

# find the closest time in a dataset from a given time

def find_closest(epoch_survey,burst_t):
    n = len(epoch_survey)

    nup = n-1
    nlow=0
    mid=np.int((nup+nlow)/2)

    while (nup - nlow)>1:
        mid= np.int((nup+nlow)/2)
        if (burst_t > epoch_survey[mid]):
            nup = nup
            nlow=mid
        else:
            nup = mid
            nlow=nlow
    
    index_burst=nup
    survey_t=epoch_survey[index_burst]

    return(survey_t,index_burst)



""" survey_data: find path to survey data on that date """

def survey_data(start_date,year,month,day):
    date_string= str(start_date.strftime("%Y%m%d"))
    survey_folder ='/data/spacecast/satellite/RBSP/emfisis/data/RBSP-A/L2'
    survey_file= 'rbsp-a_WFR-spectral-matrix-diagonal_emfisis-L2_'
    survey_path = os.path.join(survey_folder, year, month, day,survey_file + date_string + "_v*.cdf")
    
    # find the latest version
    survey_path = glob.glob(survey_path)[-1]
    
    return survey_path


# finding the magnetic data file

def magnetic_data(start_date,year,month,day):
    date_string= str(start_date.strftime("%Y%m%d"))
    survey_folder ='/data/spacecast/satellite/RBSP/emfisis/data/RBSP-A/L3'
    survey_file= 'rbsp-a_magnetometer_1sec-geo_emfisis-L3_'
    survey_path = os.path.join(survey_folder, year, month, day,survey_file + date_string + "_v*.cdf")
    # find the latest version
    survey_path = glob.glob(survey_path)[-1]
    
    return survey_path


def lanl_data(start_date,year,month,day):
    date_string= str(start_date.strftime("%Y%m%d"))
    lanl_folder ='/data/spacecast/satellite/RBSP/emfisis/data/RBSP-A/LANL/MagEphem'
    lanl_file= 'rbspa_def_MagEphem_TS04D_'

    lanl_path = os.path.join(lanl_folder, year, lanl_file + date_string + "_v*.h5")
    # find the latest version
    lanl_path = glob.glob(lanl_path)[-1]
    
    return lanl_path

# Get L4 data file (aka density file)

def l4_data(start_date,year,month,day):
    date_string= str(start_date.strftime("%Y%m%d"))
    l4_folder ='/data/spacecast/satellite/RBSP/emfisis/data/RBSP-A/L4'
    l4_file= 'rbsp-a_density_emfisis-L4_'

    l4_path = os.path.join(l4_folder, year, month, day,l4_file + date_string + "_v*.cdf")
    # find the latest version
    l4_path = glob.glob(l4_path)[-1]
    
    return l4_path


def surv_burst(time_b,surv_epoch):
    for time_s in surv_epoch:

        if (time_b < time_s <= time_b + timedelta(seconds=5.5)):
            surv_index = surv_epoch.index(time_s)
            print(surv_index,'is the index')
            break
    
    return(time_s,surv_index)


def calc_gyro(mag_file,survey_t):

    # Required constants

    q = 1.6*10**(-19)       # electron charge
    mass = 9.1*10**(-31)    # electron mass

    # get epoch in correct format for comparisons etc.

    epoch_mag = epoch_convert(mag_file['Epoch'])

    mag_field=mag_file['Magnitude']

    # Finding the closest index in this magnetometer list to the burst time object

    mag_t,mag_index = find_closest(epoch_mag,survey_t)

    # Finding the gyrofrequencies for plotting

    gyro_one= q*mag_field[mag_index]*(10**(-9))/(2*3.14*mass)
    gyro_half=0.5*gyro_one
    gyro_low=0.05*gyro_one

    return gyro_one,gyro_half,gyro_low


def omni_stats(start_date,end_date):
    import datetime
    import numpy as np

    omni_file = '/data/spacecast/wip/jinng_data/solar_wind/omni_high_res_combined_2000_to_2020.txt'
    file_format = '/data/spacecast/wip/jinng_data/solar_wind/high_res_format'
    f_format = open(f'{file_format}.txt',"r")
    line_formats=f_format.readlines()

    for line in line_formats:
        print(line)
    f_format.close()
    # From this, AE is line 11, so index line position 10
    # Initialise empty lists to store AE and omni_epoch

    AE=[]
    epoch_omni=[]
    Kp=[]
    print('the start and end are',start_date,end_date)
    start= str(start_date.strftime("%Y-%m-%d"))
    
    no_days = end_date - start_date
    no_days = int(no_days.days)

    print('but the start is',start)

    f=open(omni_file,"r")
    lines=f.readlines()
    i=0
    for line in lines:
        
        line = [x for x in line.rstrip("\n").split(" ") if x != ""]
        date = datetime.date(
            int(line[0].replace(" ", "")), 1, 1
        ) + datetime.timedelta(
            days=int(line[1]) - 1
        )
        if date == start_date:
            print(start, 'string exists in file')
            print('Line Number:', i)
            first=i
            start = start_date
            break
        i=i+1
            
    for line in lines[first:(first+(no_days*24*60))]:
    
        line = [x for x in line.rstrip("\n").split(" ") if x != ""]
        date = datetime.datetime(int(line[0].replace(" ", "")), 1, 1
        ) + datetime.timedelta(
            days=int(line[1]) - 1,
            hours=int(line[2]),
            minutes=int(line[3]),
        )
        epoch_omni.append(date)
        AE.append(float(line[10]))
    print('OMNI dataset created')
    f.close()

    return epoch_omni,AE

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

def rebin_burst(survey,mag,freq,date,
                year,month,day):

    """ 
    Read in semi-log bins and define bin edges for rebinning
    """
    
    bins = survey['WFR_bandwidth'][0]
    freq_bin, rebinned_freq = bin_edges(bins,date,year,month,day)

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
    
    return rebinned

def box_dist(Bu,Bv,Bw,B_cal,
             survey,survey_freq,
             year,month,day,
             single_day):

    N_s = len(Bu)

    ds = 16384

    n_b = int(N_s / ds)
    
    b_dist=[]

    mag_list = np.zeros((n_b,5600))
    san_check=[]
    
    for i in range(n_b):


        BuShort = Bu[i*ds:(i+1)*ds]
        BvShort = Bv[i*ds:(i+1)*ds]
        BwShort = Bw[i*ds:(i+1)*ds]
        
       # BuShort = linear_fit(BuShort,ds)
        #BvShort = linear_fit(BvShort,ds)
        #BwShort = linear_fit(BwShort,ds)

        fft_468, freq_468, wms_468,df_468,T_window = ft.fft_short(BuShort,BvShort,BwShort,B_cal)

        """ 
        take absolute values and square to get power density
        """
        total_468 = abs(fft_468) * abs(fft_468)

        """ 
        divide array by correction for spectral power lost by using a hanning window
        """
        total_468 = total_468/wms_468

        """ 
        find B field magnitude
        """
        n_f468 = len(freq_468)

        mag_468 =(total_468[0,:]+total_468[1,:]+total_468[2,:])*2*T_window  

        mag_list[i,:]=mag_468

        rebinned = rebin_burst(survey,mag_468,freq_468,single_day,year,month,day)
        san_check.append(rebinned)
        """" Do the burst integration """
        burst_int=0 
        for m in range(13,len(survey_freq)-1):
                    
            burst_int = burst_int + 0.5*(rebinned[m]+rebinned[m+1])*(survey_freq[m+1]-survey_freq[m])
                        
        b_dist.append(burst_int)
    
    return b_dist,mag_list,san_check
        

def linear_fit(data,boxsize):
    
    import numpy as np
    dt = 1./35000.
    x = np.linspace(0,boxsize*dt,boxsize)
    coef = np.polyfit(x,data,1)
    fit = np.poly1d(coef) 

    new_data = data - fit(x)

    return new_data

def process_Kletzing_windows(Bu_sample,Bv_sample,Bw_sample,
                             burst_params,
                             date_params,
                             file_access):
    
    """ Doing the FFT porcess identical to that on board the PBSP (Kletzing, 2023)
        In other words, using a 0.468s window from the beginning of the 6s burst sample, and rebinning into semu-logirtihmic bins

        Defining key parameters:

        'N_short'            - number of points in 0.468s window; integer
        'fft'                - results of the FFT; complex nT/Hz
        'freq'               - frequencies resulting from the FFT; Hz
        'wms'                - correction to the FFT power spectral densities required as a result of the Hanning window
        'Twin'               - temporal duration of window; s

        Outputs:

        'mag'               - corrected magnetic spectral density in frequency bins from FFT process; nT/Hz
        'rebinned'          - corrected magnetic spectral density rebinned in semi-logarithmic bins defined by Kletzing, 2023; nT/Hz 

    """
    
    N_short = 16384

    """ 
    calling long FFT routine

    """
    BuShort = Bu_sample[:N_short]
    BvShort = Bv_sample[:N_short]
    BwShort = Bw_sample[:N_short]
    
    fft, freq, wms,df,Twin = ft.fft_short(BuShort,BvShort,BwShort,burst_params["B_cal"])
    
    """ 
    take absolute values and square to get power density
    """
    total = abs(fft) * abs(fft)

    """ 
    divide array by correction for spectral power lost by using a hanning window
    """
    total= total/wms
    
    """ 
    find B field magnitude
    """
    
    mag =(total[0,:]+total[1,:]+total[2,:])*2*Twin

    single_day = date_params["single_day"]
    year = date_params["year"]
    month = date_params["month"]
    day = date_params["day"]

    """ 
    do rebinning in semi-logarithmic bins
    """

    rebinned = rebin_burst(file_access["survey_file"],mag,freq,single_day,year,month,day)

    params_468 = {
        "freq": freq,
        "df": df,
        "n_f": len(freq),
    } 

    return rebinned, mag, params_468


def process_survey(file_access):

    """ 
    Find survey spectral density for equivalent time

     Defining key parameters:

    'Bu2'              - power spectral density in u direction; nT/Hz 
    'Bv2'              - power spectral density in v direction; nT/Hz
    'Bw2'              - power spectral density in w direction; nT/Hz
    'Btotal'           - power spectral density magnitude; nT/Hz
    """
    survey = file_access["survey_file"]
    
    
    survey_freq = survey['WFR_frequencies'][0]

    Bu2 = survey['BuBu']
    Bv2 = survey['BvBv']
    Bw2 = survey['BwBw']

    # Define empty list for total mag field array 

    Btotal = np.zeros(survey['BuBu'].shape)
    f_band= survey['WFR_bandwidth'][:].flatten()

    # Create total mag B array

    for p in range(0,np.shape(Btotal)[0]):
        Btotal[p,:] =Bu2[p,:]+ Bv2[p,:] + Bw2[p,:]

    return Btotal


def process_small_windows(Bu_sample,Bv_sample,Bw_sample,
                            burst_params):


    ''' Doing FFTs of entire 6 second data ste, with overlapping windows of size 1024

        Defining key parameters:

        'fft'                - results of the FFT; complex nT/Hz
        'freq'               - frequencies resulting from the FFT; Hz
        'wms'                - correction to the FFT power spectral densities required as a result of the Hanning window
        'Twin'               - temporal duration of window; s
        'n_f'                - number of frequencies
        'n_bins'             - number of time bins

        Outputs:

        'mag'               - corrected magnetic spectral density in frequency bins from FFT process; nT/Hz

    '''

    """ 
    calling short FFT routine
    """
    fft_arr, freq, wms, n_bins, df,Twin = ft.fft_dat(Bu_sample,Bv_sample,Bw_sample,burst_params["B_cal"],burst_params["df_cal"], burst_params["f_max"])

    """ 
    take absolute values and square to get power density
    """
    total_m = abs(fft_arr) * abs(fft_arr)

    """ 
    divide array by correction for spectral power lost by using a hanning window
    """
    total_m = total_m/wms
    
    """ 
    find B field magnitude
    """
    n_f = len(freq)
    mag = np.zeros((n_bins,n_f))

    params_030 = {
        "freq": freq,
        "n_bins": n_bins,
        "df": df,
        "n_f": n_f,
        "n_bins": n_bins
    } 

    # There are 2*N -1 boxes
    for n in range(n_f):
        for m in range(n_bins):

            # Multiply by 2 to account for negative freqiencies from FFT 
            # Divide by frequency step size to have a PSD in units/Hz 
            mag[m,n]=(total_m[m,n,0]+total_m[m,n,1]+total_m[m,n,2])*2*Twin  

    return mag,params_030

""" FFT with sliding overlapping windows """

def process_sliding_windows(Bu_sample,Bv_sample,Bw_sample,slider,burst_params):
    from scipy.fft import fft, fftfreq,rfftfreq
    import numpy as np

    f_s = 1.0/35000.0 
    N = len(Bu_sample)                                                       # Number of sample points
    box_size = 1024 #int(N/n_windows)                                          # Number of samples in each box 
  
    n_bins = int((N - (box_size - slider))/slider)


    """ 
    Frequencies from FFT boxes
    """
    
    freq = rfftfreq(box_size, f_s)[:box_size//2]
    n_f = len(freq)

    df = freq[1] - freq[0]


    """ 
    Putting the calibration coeffcients into complex form
    """
    cal_step = int(df/burst_params["df_cal"])
    Bcal_c=[]
    
    for i in range(n_f):
        f = freq[i]
        if (f < burst_params["f_max"]):
            Bcal_c.append(complex(burst_params["B_cal"][i*cal_step,0],burst_params["B_cal"][i*cal_step,1]))
        else:
            break

    n_f = len(Bcal_c)

    # Cutting off frequencies where the reciever stopped picking them up (fmax from callibration specifications)
    freq = freq[0:n_f]                                  
    """ 
    Do FFTs - remember to normalise by N and multiply by hanning window and complex correction
    """

    fft_box = np.zeros((n_bins,n_f,3),dtype = complex)

    w = np.hanning(box_size)                                            # hanning window
    wms = np.mean(w**2)                                                 # hanning window correction

    # Remember to normalise the box size
        
    lower_edge = 0
    upper_edge = box_size
    i = 0 

    while upper_edge<N:

        fft_box[i,:,0] = (fft(w*(1/box_size)*Bu_sample[lower_edge:upper_edge]))[:n_f]*Bcal_c  
        fft_box[i,:,1] = (fft(w*(1/box_size)*Bv_sample[lower_edge:upper_edge]))[:n_f]*Bcal_c
        fft_box[i,:,2] = (fft(w*(1/box_size)*Bw_sample[lower_edge:upper_edge]))[:n_f]*Bcal_c
        
        lower_edge += slider
        upper_edge += slider

        i += 1

    fft_av = fft_box
   
    T_window = box_size*f_s

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
    PSD = np.zeros((n_bins,n_f))

    for n in range(n_f):
        for m in range(n_bins):

            # Multiply by 2 to account for negative freqiencies from FFT 
            # Divide by frequency step size to have a PSD in units/Hz 
            PSD[m,n]=(total[m,n,0]+total[m,n,1]+total[m,n,2])*2*T_window  
    
    """ 
    define the time array for plotting
    """
    duration = N*f_s
    t_array = np.linspace(0.,duration- f_s, n_bins)

    params_030 = {
        "freq": freq,
        "n_bins": n_bins,
        "df": df,
        "n_f": n_f,
        "n_bins": n_bins,
        "time_array": t_array
    } 

    return PSD,params_030


def integrate_in_rebin(B,survey_freq):

    integral=0 
    for m in range(13,len(survey_freq)-1):
                        
        integral = integral + 0.5*(B[m]+B[m+1])*(survey_freq[m+1]-survey_freq[m])

    return integral

def integrate_in_small_windows(B,params_030):

    freq = params_030["freq"]
    n_f = params_030["n_f"]
    n_bins = params_030["n_bins"]

    frequency_integral = np.zeros(n_bins)
    
    for n in range(n_bins):
        high_res_bint=0 

        # Integrate in frequency 
        for m in range(1,n_f-1):
                
            high_res_bint = high_res_bint + 0.5*(B[n,m]+B[n,m+1])*(freq[m+1]-freq[m])

        frequency_integral[n] = high_res_bint

    max_amplitude = np.max(np.sqrt(frequency_integral))

    integral_statistics = {"maximum amplitude": max_amplitude,
                           "mean power": np.mean(frequency_integral),
                           "median power": np.median(frequency_integral),
                           "power vairaince": np.std(frequency_integral)**2,
                           "frequency integrated power": frequency_integral}

    return integral_statistics



def spectrum_averaging(mag_030,params_030,file_access,date_params):

    ''' Averaging over either sectiobs or the shole of the 0.03s window results

        Defining key parameters:

        'averaged_spec'                - the average spectrum over the full 6s; nT/Hz
        'av_rebinned'                  - the average spectrum over the full 6s, rebinned in the semi_logarithmic bins; nT/Hz
        'averaged_spec_win'            - the average spectrum over first 0.468s; nT/Hz
        'av_rebinned_win'              - the average spectrum over first 0.468s, rebinned in the semi_logarithmic bins; nT/Hz

    '''

    single_day = date_params["single_day"]
    year = date_params["year"]
    month = date_params["month"]
    day = date_params["day"]

    """ over full 6s """

    averaged_spec = np.zeros(params_030["n_f"])    

    for m in range(params_030["n_f"]):
        averaged_spec[m] = np.mean(mag_030[:,m])

    av_rebinned = rebin_burst(file_access["survey_file"],averaged_spec,params_030["freq"],single_day,year,month,day)
   
    """ over first 0.468s """

    averaged_spec_win = np.zeros(params_030["n_f"])    
    for m in range(params_030["n_f"]):
        averaged_spec_win[m] = np.mean(mag_030[:32,m])

    av_rebinned_win = rebin_burst(file_access["survey_file"],averaged_spec_win,params_030["freq"],single_day,year,month,day)
    
    return av_rebinned, av_rebinned_win
    