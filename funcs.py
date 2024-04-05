import numpy as np
import os
import glob
from datetime import datetime,date,timedelta
#import datetime
import pandas as pd
from spacepy import pycdf
import fft_box as ft

global_constants = {"Duration": 4.698, 
                    "Electron q": 1.6*10**(-19),
                    "Electron m":9.1*10**(-31),
                    "Convert to nT": (10**(-9)),
                    "Pi": 3.14159,
                    "Slider": int(1024/4)}  

def get_date_string(date):
        ''' Method that rpovides date strings
        Outputs:
    
        date_string: string object of date
        year: string object of date year
        month: string object of date month
        day: string object of date day '''

        date_string= str(date.strftime("%Y%m%d"))

        if (date.day <10):
            day = "0"+str(date.day)
        else:
            day = str(date.day)


        if (date.month<10):
            month = "0"+str(date.month)
        else:
            month = str(date.month)
        

        year = str(date.year)
        
        return date_string,year,month,day

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

def h5_time_conversion(bytes_time) -> datetime:
    """
    Converts ISO 8601 datetime, which contains leap seconds, into Python
    datetime object

    :param bytes_time: time in bytes to convert
    :return: datetime in Python datetime object
    """
    date_str, time_str = bytes_time.decode("UTF-8").replace("Z", "").split("T")
    year, month, day = date_str.split("-")
    date = datetime(int(year), int(month), int(day))
    no_ms_time, dot, ms = time_str.partition(".")
    hours, minutes, seconds = no_ms_time.split(":")
    date = date + timedelta(
        hours=int(hours), minutes=int(minutes), seconds=int(seconds)
    )
    if dot:
        date = date + timedelta(milliseconds=int(ms))
    return date

def get_epoch(epoch):
    # get epoch in correct format for comparisons etc.
        epoch_new = []

        # First - saving the epoch elements to a list as a string - acording to a particular desired format
        for i in range(len(epoch)):
            epoch_new.append(datetime.strftime(epoch[i],'%Y-%m-%d %H-%M-%S'))

        # Chaning back to a datetime object of same format    
        for i in range(len(epoch_new)):
            epoch_new[i] = datetime.strptime(epoch_new[i],'%Y-%m-%d %H-%M-%S')

        return epoch_new

def find_closest(epoch_known, wanted_time):

        n = len(epoch_known)

        nup = n-1
        nlow=0
        mid=np.int((nup+nlow)/2)

        while (nup - nlow)>1:
            mid= np.int((nup+nlow)/2)
            if (wanted_time > epoch_known[mid]):
                nup = nup
                nlow=mid
            else:
                nup = mid
                nlow=nlow
        
        wanted_index = nup
        corresponding_time = epoch_known[wanted_index]

        return(corresponding_time,wanted_index)


class omni_dataset:
    ''' Class for creating omni_dataset over full duration specified '''

    def __init__(self,start,end):
        self.start_date = start
        self.end_date = end

    @property
    def omni_stats(self):
        
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
            print('the start and end are',self.start_date,self.end_date)
            start_d = str(self.start_date.strftime("%Y-%m-%d"))
            
            no_days = self.end_date - self.start_date
            no_days = int(no_days.days)


            f=open(omni_file,"r")
            lines=f.readlines()
            i=0
            for line in lines:
                
                line = [x for x in line.rstrip("\n").split(" ") if x != ""]
                Date = date(
                    int(line[0].replace(" ", "")), 1, 1
                ) + timedelta(
                    days=int(line[1]) - 1
                )
                if Date == self.start_date:
                    print(start_d, 'string exists in file')
                    print('Line Number:', i)
                    first=i
                    start_d = self.start_date
                    break
                i=i+1
                    
            for line in lines[first:(first+(no_days*24*60))]:
            
                line = [x for x in line.rstrip("\n").split(" ") if x != ""]
                Date = datetime(int(line[0].replace(" ", "")), 1, 1
                ) + timedelta(
                    days=int(line[1]) - 1,
                    hours=int(line[2]),
                    minutes=int(line[3]),
                )
                epoch_omni.append(Date)
                AE.append(float(line[10]))
            print('OMNI dataset created')
            f.close()

            return AE,epoch_omni
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



    
class AccessLANLAttrs:

    ''' Class for finding and working on Survey data attributes
     
      PARAMETERS:
      survey_cdf: A CDF containing all survey data '''

    def __init__(self, lanl_data):
        self.data = lanl_data

    @property
    def L_star(self):

        # get L*
        lstar = np.array(self.data['Lstar'][:, 0])
        
        return lstar
    
    @property
    def MLT(self):

        # get MLT
        MLT = np.array(self.data['EDMAG_MLT'])

        return MLT
    
    @property
    def MLAT_N_S(self):

        # get all MLAT
        MLAT = np.array(self.data['EDMAG_MLAT'])

        # want to plot north and south on same axis - so split up
        south_mask = MLAT<5.
        north_mask = MLAT>0.
        MLAT_north = np.where(north_mask,MLAT, np.nan)
        MLAT_south = np.where(south_mask,MLAT, np.nan)

        return MLAT_north,MLAT_south
    
    @property
    def epoch(self):

        lanl_times=np.array([h5_time_conversion(x) for x in self.data['IsoTime']])

        return lanl_times
    
    @property
    def day_limits(self):

        L_star = self.L_star
        
        left=np.array([h5_time_conversion(x) for x in self.data['IsoTime']])[L_star >= 0][0]
        right=np.array([h5_time_conversion(x) for x in self.data['IsoTime']])[L_star >= 0][-1]

        return left,right


class AccessSurveyAttrs:
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
    
    @property
    def Bmagnitude(self):

        Bu2 = self.survey_cdf['BuBu']
        Bv2 = self.survey_cdf['BvBv']
        Bw2 = self.survey_cdf['BwBw']

        # Define empty list for total mag field array 

        Btotal = np.zeros(self.survey_cdf['BuBu'].shape)

        # Create total mag B array

        for p in range(0,np.shape(Btotal)[0]):
            Btotal[p,:] =Bu2[p,:]+ Bv2[p,:] + Bw2[p,:]

        return Btotal
    
    @property
    def Emagnitude(self):

        Eu2 = self.survey_cdf['EuEu']
        Ev2 = self.survey_cdf['EvEv']
        Ew2 = self.survey_cdf['EwEw']

        # Define empty list for total mag field array 

        Etotal = np.zeros(Eu2.shape)

        # Create total mag B array

        for p in range(0,np.shape(Etotal)[0]):
            Etotal[p,:] =Eu2[p,:]+ Ev2[p,:] + Ew2[p,:]

        return Etotal
    
    def epoch_convert(self):
    # get epoch in correct format for comparisons etc.
        epoch_edit = []
        epoch = self.survey_cdf['Epoch']

        # First - saving the epoch elements to a list as a string - acording to a particular desired format
        for i in range(len(epoch)):
            epoch_edit.append(datetime.strftime(epoch[i],'%Y-%m-%d %H-%M-%S'))

        # Chaning back to a datetime object of same format    
        for i in range(len(epoch_edit)):
            epoch_edit[i] = datetime.strptime(epoch_edit[i],'%Y-%m-%d %H-%M-%S')

        return epoch_edit
    
class AccessL3Attrs:
    ''' Class for finding and working on L3 data attributes
     
      PARAMETERS:
      mag_cdf: A CDF containing all L3 data '''

    def __init__(self, mag_file):
        self.mag_cdf = mag_file

    @property
    def Bmagnitude(self):
        # Frequency in Hz
        Bmagnitude = self.mag_cdf['Magnitude']

        return Bmagnitude
    
    @property
    def epoch(self):
        # Epoch in DateTime format of: 
        epoch = self.mag_cdf['Epoch']

        return epoch
    
    def epoch_convert(self):
    # get epoch in correct format for comparisons etc.
        epoch_edit = []
        epoch = self.mag_cdf['Epoch']

        # First - saving the epoch elements to a list as a string - acording to a particular desired format
        for i in range(len(epoch)):
            epoch_edit.append(datetime.strftime(epoch[i],'%Y-%m-%d %H-%M-%S'))

        # Chaning back to a datetime object of same format    
        for i in range(len(epoch_edit)):
            epoch_edit[i] = datetime.strptime(epoch_edit[i],'%Y-%m-%d %H-%M-%S')

        return(epoch_edit)
    
    @property
    def f_ce(self):
        # Finding the gyrofrequencies for plotting
       
        gyro_1 = np.zeros(self.Bmagnitude.shape)
        gyro_05 = np.zeros(self.Bmagnitude.shape)
        gyro_005 = np.zeros(self.Bmagnitude.shape)
        
        # Clean magnetometer data
        mag_cleaned = self.clean_magnetometer(self.Bmagnitude)


        for i in range(0,len(gyro_1)):
            gyro_1[i] = global_constants["Electron q"]*mag_cleaned[i]*global_constants["Convert to nT"]/(2*global_constants["Pi"]*global_constants['Electron m'])
            gyro_05[i] = 0.5*gyro_1[i]
            gyro_005[i] = 0.05*gyro_1[i]
        
        return gyro_1, gyro_05, gyro_005

    def clean_magnetometer(self,unclean_data):

        magfill = self.mag_cdf['magFill']
        magInval = self.mag_cdf['magInvalid']
        magCal = self.mag_cdf['calState']

        # Find the indicies where we have invalid data, save to 'whereFill' 
        whereFill = []
        for i in range(len(magfill)):
            if (magfill[i] == 1 or magInval[i]==1 or magCal[i]==1):
                whereFill.append(i)
                
        # Make unclean data into list
        magnitude = np.array(unclean_data).tolist()

        # Use list of indicies to set all invlaid data points to NaNs
        for ele in sorted(whereFill, reverse = True):
            
            magnitude[ele] = np.nan
        
        return magnitude


class AccessL4Attrs:
    ''' Class for finding and working on L4 data attributes.
     
      PARAMETERS:
      density_cdf: A CDF containing all L4 data '''

    def __init__(self, density_file):
        self.density_cdf = density_file

    @property
    def density(self):
        # density in cm^(-3)
        density = self.density_cdf['density']

        return density
    
    @property 
    def f_pe(self):
        # plasma frequency in Hz
        f_pe = self.density_cdf['fpe']

        return f_pe
        
    @property
    def epoch(self):
        # Epoch in DateTime format of: 
        epoch = self.density_cdf['Epoch']

        return epoch
    
    def epoch_convert(self):
    # get epoch in correct format for comparisons etc.
        epoch_edit = []
        epoch = self.density_cdf['Epoch']

        # First - saving the epoch elements to a list as a string - acording to a particular desired format
        for i in range(len(epoch)):
            epoch_edit.append(datetime.strftime(epoch[i],'%Y-%m-%d %H-%M-%S'))

        # Chaning back to a datetime object of same format    
        for i in range(len(epoch)):
            epoch_edit[i] = datetime.strptime(epoch_edit[i],'%Y-%m-%d %H-%M-%S')

        return(epoch_edit)

class cross_dataset:
    ''' A class for performing operations across datasets '''
    def __init__(self,survey_data, l4_data, burst_time):

        self.survey_epoch = survey_data['Epoch']
        self.mag_epoch = l4_data['Epoch']
        self.Bmag = l4_data['Magnitude']
        self.burst_time = burst_time
    
    def calc_gyro(self):
        # get epoch in correct format for comparisons etc.

        epoch_mag = get_epoch(self.mag_epoch)

        mag_field=self.Bmag

        # Finding the closest index in this magnetometer list to the burst time object

        mag_t,mag_index = find_closest(epoch_mag,self.burst_time)

        # Finding the gyrofrequencies for plotting

        gyro_one= global_constants["Electron q"]*mag_field[mag_index]/(2*global_constants["Pi"]*global_constants["Electron m"]) # in T
        gyro_one = gyro_one*global_constants["Convert to nT"]                                                                   # in nT
        gyro_half=0.5*gyro_one
        gyro_low=0.05*gyro_one


        return gyro_one,gyro_half,gyro_low
    
    
    def get_epoch(epoch):
    # get epoch in correct format for comparisons etc.
        epoch_new = []

        # First - saving the epoch elements to a list as a string - acording to a particular desired format
        for i in range(len(epoch)):
            epoch_new.append(datetime.strftime(epoch[i],'%Y-%m-%d %H-%M-%S'))

        # Chaning back to a datetime object of same format    
        for i in range(len(epoch_new)):
            epoch_new[i] = datetime.strptime(epoch_new[i],'%Y-%m-%d %H-%M-%S')

        return epoch_new
    
    
    
            
class PerformFFT:
    ''' class for storing all methods related to performing FFTs'''

    def __init__(self,
               waveform_samples,
               burst_params,
               slider):

        self.waveform_samples = waveform_samples # Bu, Bv, Bw
        self.burst_params = burst_params         # B_cal, df_cal, f_max
        self.slider = slider
    
    def process_Kletzing_windows(self):
        from scipy.fft import fft,rfftfreq
        import numpy as np
    
        """ Doing the FFT porcess identical to that on board the PBSP (Kletzing, 2023)
        In other words, using a 0.468s window from the beginning of the 6s burst sample, and rebinning into semu-logirtihmic bins

        Defining key parameters:

        'N_short'            - number of points in 0.468s window; integer
        'fft'                - results of the FFT; complex nT/Hz
        'freq'               - frequencies resulting from the FFT; Hz
        'wms'                - correction to the FFT power spectral densities required as a result of the Hanning window
        'Twin'               - temporal duration of window; s

        Outputs:

        'PSD'               - corrected magnetic spectral density in frequency bins from FFT process; nT/Hz

        """

        """ 
        calling long FFT routine

        """
        Bu = self.waveform_samples["Bu"]
        Bv = self.waveform_samples["Bv"]
        Bw= self.waveform_samples["Bw"]
        
        f_s = 1.0/35000.0 
        N = len(Bu)                                                     # Number of sample points
        box_size = 16384                                                # Number of samples in each box 
        n_bins =  int(N/box_size)                                                      # Number of points in temporal space


        """ 
        Frequencies from FFT boxes
        """
        
        freq = rfftfreq(box_size, f_s)[:box_size//2]
        n_f = len(freq)

        df = freq[1] - freq[0]

        """ 
        Putting the calibration coeffcients into complex form
        """
        cal_step = int(df/self.burst_params["df_cal"])
        Bcal_c=[]
        
        for i in range(n_f):
            f = freq[i]
            if (f < self.burst_params["f_max"]):
                Bcal_c.append(complex(self.burst_params["B_cal"][i*cal_step,0],self.burst_params["B_cal"][i*cal_step,1]))
            else:
                break
        
        # Resetting number of frequencies to f_max - i.e. Cutting off frequencies where the reciever stopped picking them up (fmax from callibration specifications))
        n_f = len(Bcal_c)
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

            fft_box[i,:,0] = (fft(w*(1/box_size)*Bu[lower_edge:upper_edge]))[:n_f]*Bcal_c  
            fft_box[i,:,1] = (fft(w*(1/box_size)*Bv[lower_edge:upper_edge]))[:n_f]*Bcal_c
            fft_box[i,:,2] = (fft(w*(1/box_size)*Bw[lower_edge:upper_edge]))[:n_f]*Bcal_c
            
            lower_edge += self.slider
            upper_edge += self.slider

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
        duration = n_bins*T_window
        t_array = np.linspace(0.,duration, n_bins)

        FFTs = {"Frequencies": freq,
            "PSD": PSD,
            "Time": t_array,
            "Flag": 'Kletzing'} 

        return FFTs
    

    """ FFT with sliding overlapping windows """

    def process_sliding_windows(self):
        from scipy.fft import fft,rfftfreq
        import numpy as np

        Bu = self.waveform_samples["Bu"]
        Bv = self.waveform_samples["Bv"]
        Bw = self.waveform_samples["Bw"]

        f_s = 1.0/35000.0 
        N = len(Bu)                                                     # Number of sample points
        box_size = 1024                                                 # Number of samples in each box 
        n_bins = int((N - (box_size - self.slider))/self.slider)        # Number of points in temporal space


        """ 
        Frequencies from FFT boxes
        """
        
        freq = rfftfreq(box_size, f_s)[:box_size//2]
        n_f = len(freq)

        df = freq[1] - freq[0]


        """ 
        Putting the calibration coeffcients into complex form
        """
        cal_step = int(df/self.burst_params["df_cal"])
        Bcal_c=[]
        
        for i in range(n_f):
            f = freq[i]
            if (f < self.burst_params["f_max"]):
                Bcal_c.append(complex(self.burst_params["B_cal"][i*cal_step,0],self.burst_params["B_cal"][i*cal_step,1]))
            else:
                break
        
        # Resetting number of frequencies to f_max - i.e. Cutting off frequencies where the reciever stopped picking them up (fmax from callibration specifications))
        n_f = len(Bcal_c)
        freq = freq[0:n_f]  

        """ 
        Do FFTs - remember to normalise by N and multiply by hanning window and complex correction
        """

        fft_box = np.zeros((n_bins,n_f,3),dtype = complex)

        w = np.hanning(box_size)                                            # hanning window
        wms = np.mean(w**2)                                                 # hanning window correction

        # Do sliding window FFTs

            
        lower_edge = 0              # Lower edge of each box
        upper_edge = box_size       # Upper edge of each box
        i = 0 

        while upper_edge<N:

            fft_box[i,:,0] = (fft(w*(1/box_size)*Bu[lower_edge:upper_edge]))[:n_f]*Bcal_c  
            fft_box[i,:,1] = (fft(w*(1/box_size)*Bv[lower_edge:upper_edge]))[:n_f]*Bcal_c
            fft_box[i,:,2] = (fft(w*(1/box_size)*Bw[lower_edge:upper_edge]))[:n_f]*Bcal_c
            
            lower_edge += self.slider
            upper_edge += self.slider

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
        n_468 = int((duration/0.468)*n_bins)

        FFTs = {"Frequencies": freq,
            "PSD": PSD,
            "PSD_0468s": PSD[0:n_468,:],
            "Time": t_array,
            "Flag": 'sliding'} 

        return FFTs
    

# Class for the reuslts from windowed FFTs 

class AnalyseFFT:
    ''' class for storing all methods related to analysing FFT data
    
            PARAMETERS:
            PSD_Kletzing: 2D array of FFT spectra, with shape (n_bins, frequency_bins) correspnding to 0.468s windows, no overlap
            Frequency_Kletzing: central frequencies for each frequency band (Hz) '''
            
            

    def __init__(self,FFTs):
        self.PSD = FFTs["PSD"]
        self.Frequency = FFTs["Frequencies"]
        self.Times = FFTs["Times"]


    def show_dimensions(self):
        print(f"The dimensions (time bins, frequency bins) of FFT_Kletzing spectra is: {self.FFT.shape}")


    def show_frequency(self):
        print(f"The number of frequencies is: {len(self.frequency)}")

    def n_f(self):
        # integer number of frequency bins 

        n_f = len(self.frequency)
        return n_f
    
    def n_bins(self):
        # integer number of frequency bins 

        n_f = len(self.frequency)
        return n_f
    
    def integrate_in_frequency(self):

        n_f = len(self.frequency)
        n_bins = np.shape(self.PSD)[0]

        frequency_integral = np.zeros(n_bins)
        
        for n in range(n_bins):
            high_res_bint=0 

            # Integrate in frequency 
            for m in range(0,n_f-1):
                    
                high_res_bint = high_res_bint + 0.5*(self.PSD[n,m]+self.PSD[n,m+1])*(self.Frequency[m+1]-self.Frequency[m])

            frequency_integral[n] = high_res_bint

        # save the stats for the first 0.468s too
        if n_bins == 12: 
            # if FFT_Kletzing, then it's just the first value
            mean_468s = frequency_integral[0]
            median_468s = 0.
            Kletzing = True
        else:
            n_468 = int((global_constants["Duration"]/0.468)*n_bins)
            mean_468s = np.mean(frequency_integral[:n_468])
            median_468s = np.median(frequency_integral[:n_468])
            Kletzing = False

        max_amplitude = np.max(np.sqrt(frequency_integral))

        FFT_stats = {"Maximum amplitude": max_amplitude,
                            "Mean power": np.mean(frequency_integral),
                            "Median power": np.median(frequency_integral),
                            "Power vairaince": np.std(frequency_integral)**2,
                            "Frequency integrated power": frequency_integral,
                            "Mean power over 0.468s": mean_468s,
                            "Median power over 0.468s": median_468s,
                            "Kletzing flag": Kletzing}

        return FFT_stats
    

class rebin_data:

    def __init__(self,survey_file, PSD, PSD_Frequencies, date_params):

        # Needed from survey 
        self.bin = survey_file['WFR_bandwidth'][0]
        self.Survey_frequencies = survey_file['WFR_frequencies'][0]

        # PSD that is to be rebinned
        self.PSD = PSD
        self.frequencies = PSD_Frequencies

        # date parameters
        self.date_params = date_params



    def rebin_burst(self):

        """ 
        Read in semi-log bins and define bin edges for rebinning
        """
        freq_bin, rebinned_freq = self.bin_edges(self.bin,self.date_params["single_day"],self.date_params["year"],self.date_params["month"],self.date_params["day"])

        """ 
        Doing the rebinning
        """
        
        # Create dataframe
        rebin_dat=pd.DataFrame()

        rebin_dat['Data'] = self.PSD
        
        # Create and save frequencies to one column
        rebin_dat['Frequency']= self.frequencies
        
        """
        pd.cut() bins all frequencies according to defined semi_log bins
        groupby() groups all data in frame by these bines
        then select the DATA and take the MEAN in each bin
        to_numpy saves this to an array
        """
        
        rebinned=rebin_dat.groupby(pd.cut(rebin_dat.Frequency,bins=freq_bin)).Data.mean().to_numpy()
        
        return rebinned

    def bin_edges(self):
        """ 
        setting the bin edges 
        min_bin: lower edge of first bin 
        """

        min_bin = self.Survey_frequencies[0]-(self.bin[0]/2)

        freq_bin = []
        freq_bin.append(min_bin)

        # Starting from the minimum bin, add all widths to get lower and upper bands on all semi_logarithmic bins

        for i in range(0,65):
            
            freq_bin.append(freq_bin[i]+self.bin[i])
            min_bin=freq_bin[i]
        
        
        return freq_bin, self.Survey_frequencies

class createPSDFiles:

    ''' class for creating CDFs containg burst FFT data + statistics
    
            PARAMETERS:
             '''

    def __init__(self,FFTs,date_params,fces):

        self.PSD = FFTs["PSD"]
        self.Frequency = FFTs["Frequency"]
        self.Time = FFTs["Times"]
        self.timedate = date_params
        self.Kletzing = FFTs["Flag"]
        self.fce = fces[0]
        self.fce_05 = fces[1]
        self.fce_005 = fces[2]
    
    
    def save_FFT(self):
        
        # What is the date and time of this burst? Make a directory for that day
        file_path = '/data/emfisis_burst/wip/rablack75/rablack75/simple_FFT/PSDs/'+self.date_params['year']
        os.makedirs(file_path, exist_ok=True)
        os.makedirs(file_path+'/'+self.date_params['month'], exist_ok=True)
        os.makedirs(file_path+'/'+self.date_params['month'] + self.date_params['day'], exist_ok=True)
        # Is it Kletzing windows or normal?
        if self.Kletzing == 'Kletzing':
            cdf_name = file_path+'/'+self.date_params['month']+'/' + 'PSD_Kletzing_'+str(self.date_params['burst_time'])+'.cdf'
        else:
            cdf_name = file_path+'/'+self.date_params['month']+'/' + 'PSD_'+str(self.date_params['burst_time'])+'.cdf'

        # Create CDF for Burst Time
        cdf = pycdf.CDF(cdf_name, '')

        # Save main datasets: the PSD, the time steps, and the frequencies
        cdf['Time'] = self.Time
        cdf['PSD'] = self.PSD
        cdf['Frequencies'] = self.Frequency

        # Set units for the above
        cdf['PSD'].attrs['units'] = 'nT/Hz'
        cdf['Frequencies'].attrs['units'] = 'Hz'
        cdf['Time'].attrs['units'] = 's'

        # Set dependencies for the PSD
        cdf['PSD'].attrs['Dependency1'] = 'Time'
        cdf['PSD'].attrs['Dependency2'] = 'Frequency' 

        # important quantities
        cdf['BurstDatetime'] = self.timedate['burst_time']
        cdf['fce'] = self.fce
        cdf['fce_05'] = self.fce_05
        cdf['fce_005'] = self.fce_005

        cdf['fce'].attrs['units'] = 'Hz'
        cdf['fce_05'].attrs['units'] = 'Hz'
        cdf['fce_005'].attrs['units'] = 'Hz'
        cdf['BurstDatetime'].attrs['units'] = 'DateTime'

        # set author
        cdf.attrs['Author'] = 'Rachel Black'
        
        cdf.close()
        print("CDF saved")