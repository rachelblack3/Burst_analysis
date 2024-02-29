import numpy as np
from matplotlib import pyplot as plt, animation
from metpy.plots import add_timestamp
import os

# for plotting
import matplotlib.colors as mcolors
import matplotlib.dates as mdates

# for dealing h5py
import time_conversions

# for datetime manipulating 
from datetime import datetime

# my own misc. functions
import funcs as funcs

# Setting the upper limit for time array 
N = 208896
T = 1/35000
T_high = N*T
print('The upper time limit is', T_high)

def mlat_defs(lanl_d):

    MLAT = np.array(lanl_d['EDMAG_MLAT'])
    south_mask = MLAT<5.
    north_mask = MLAT>0.
    MLAT_north = np.where(north_mask,MLAT, np.nan)
    MLAT_south = np.where(south_mask,MLAT, np.nan)

    return MLAT_north,MLAT_south


def summary_plot(freq_set,data_set,
                 gyro_half_set,gyro_low_set,gyro_1,
                 freq_survey,epoch_survey,Btotal,
                 burst_times,
                 lanl_d,AE,epoch_omni,
                 year,month,day,
                 survey_int,
                 burst_int,
                 burst_dist,
                 highB_dist,
                 mag_lists,freq_468,
                 av_specs,av_specs_half,surv_set):
    
    
    # defining the burst timestamps
    nt_array = np.linspace(0,T_high,407)

    # get MLATS for day

    lanl_times=np.array([time_conversions.h5_time_conversion(x) for x in lanl_d['IsoTime']])
    MLAT_north, MLAT_south = mlat_defs(lanl_d)

    # create figure
    fig, axs = plt.subplots(8,2,gridspec_kw={"width_ratios":(1,0.05)})

    # Name plot axis
    axs_list = axs[:, 0]

    # Name colourbar axis (the third/fourth one will be made invisible)
    cax_list = axs[:,1]

    # 1. Survey

    def survey_plot(epoch_survey, freq_survey,
                     Btotal, burst_times,
                     plot_number,axs_list,cax_list):


        colorbar_norm2 = mcolors.LogNorm(vmin=10**(-10), vmax=10**(-5))

        survey_img = axs_list[plot_number].pcolormesh(epoch_survey, freq_survey, np.transpose(Btotal), norm=colorbar_norm2, cmap="jet")
        
        axs_list[plot_number].set_yscale('log')
        axs_list[plot_number].set_ylabel('Frequency (Hz)',fontsize=8)
        axs_list[plot_number].set_xlabel('UTC',fontsize=8)
        plt.colorbar(survey_img,label=r'$nT^2/Hz$',cax=cax_list[plot_number])

        axs_list[plot_number].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        hours = mdates.HourLocator(interval=1)
        axs_list[plot_number].xaxis.set_minor_locator(hours)

        vertical_change=axs_list[plot_number].axvline(burst_times[0],ymin=0,ymax=10**4,linestyle='dashed',linewidth=2,color='w')

        return survey_img, vertical_change
    
    def burst_small(freq_set,data_set,
                    gyro_1,gyro_low_set,
                    plot_number,axs_list,cax_list):

        axs_list[plot_number].set_xlabel('Time (s)',fontsize=8)
        axs_list[plot_number].set_ylabel('Frequency (Hz)',fontsize=8)
        colorbar_norm = mcolors.LogNorm(vmin=10**(-10), vmax=10**(-5))
        burst_img = axs_list[plot_number].pcolormesh(nt_array, freq_set[0] , data_set[0].T, norm=colorbar_norm, cmap="viridis")
        plt.colorbar(burst_img,label=r'$nT^2/Hz$',cax=cax_list[plot_number])

        gyro_line1=axs_list[plot_number].axhline(gyro_1[0]/2.,np.min(nt_array),np.max(nt_array),color='white',label=r'$0.5\ f_{ce}$')
        #gyro_line1=ax1.axhline(gyro_half_set[0],np.min(nt_array),np.max(nt_array),color='white',label=r'$0.5\ f_{ce}$')
        gyro_line2=axs_list[plot_number].axhline(gyro_low_set[0],np.min(nt_array),np.max(nt_array),color='white',label=r'$0.05\  f_{ce}$')
        axs_list[plot_number].set_xlabel('Time (s)',fontsize=8)
        axs_list[plot_number].set_ylim(0,10000)

        tex1 = r'$0.5\ f_{ce}$'
        text_1=axs_list[plot_number].text(4.5, gyro_half_set[0]+200, tex1, fontsize=10, va='bottom',color='white')
        tex2 = r'$0.05\ f_{ce}$'
        text_2=axs_list[plot_number].text(4.5, gyro_low_set[0]+200, tex2, fontsize=10, va='bottom',color='white')
        # Do the other plot 

        return burst_img,gyro_line1,gyro_line2,text_1,text_2
    
    def burst_large(freq_468, mag_lists,
                    gyro_1, gyro_low_set,
                    burst_dist,
                    plot_number,axs_list,cax_list):
        

        axs_list[plot_number].set_xlabel('Time (s)',fontsize=8)
        axs_list[plot_number].set_ylabel('Frequency (Hz)',fontsize=8)
        
        n_b = len(burst_dist[0])
        T = 1./35000.
        win = 16384
        T_win = T*win
        nt = np.linspace(T_win/2., T_high-T_win/2.,n_b)
        colorbar_norm = mcolors.LogNorm(vmin=10**(-10), vmax=10**(-5))
    

        burst_samps = axs_list[plot_number].pcolormesh(nt, freq_468 , mag_lists[0].T, norm=colorbar_norm, cmap="viridis")
        plt.colorbar(burst_samps,label=r'$nT^2/Hz$',cax= cax_list[plot_number])

        gyro_lineb=axs_list[plot_number].axhline(gyro_1[0]/2.,0,5.968,color='white',label=r'$0.5\ f_{ce}$')
        
        gyro_lineb2=axs_list[plot_number].axhline(gyro_low_set[0],0,5.968,color='white',label=r'$0.05\  f_{ce}$')

        tex1 = r'$0.5\ f_{ce}$'
        text_1=axs_list[plot_number].text(4.5, gyro_half_set[0]+200, tex1, fontsize=10, va='bottom',color='white')
        tex2 = r'$0.05\ f_{ce}$'
        text_2=axs_list[plot_number].text(4.5, gyro_low_set[0]+200, tex2, fontsize=10, va='bottom',color='white')

        return burst_samps,gyro_lineb,gyro_lineb2,text_1,text_2
    

    
    def space_locations(lanl_d,burst_times, 
                        plot_number,axs_list,cax_list):
        
        lstar = np.array(lanl_d['Lstar'][:, 0])
        MLT = np.array(lanl_d['EDMAG_MLT'])

        axs_list[plot_number].plot(np.array([time_conversions.h5_time_conversion(x) for x in lanl_d['IsoTime']])[lstar >= 0], 
            lstar[lstar >= 0])
        axs_list[plot_number].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        hours = mdates.HourLocator(interval=1)
        axs_list[plot_number].xaxis.set_minor_locator(hours)
        vertical_change2=axs_list[plot_number].axvline(burst_times[0],ymin=0,ymax=10**4,linestyle='dashed',linewidth=2,color='black')
        axs_list[plot_number].set_xlim(left=np.array([time_conversions.h5_time_conversion(x) for x in lanl_d['IsoTime']])[lstar >= 0][0],right=np.array([time_conversions.h5_time_conversion(x) for x in lanl_d['IsoTime']])[lstar >= 0][-1])
        axs_list[plot_number].set_ylabel('L*',color='blue',fontsize=8)
        axs_list[plot_number].set_xlabel('UTC',fontsize=8)
        axs_list[plot_number].tick_params(axis='y', colors='blue')


        ax3_2 = axs_list[plot_number].twinx()  # instantiate a second axes that shares the same x-axis

        color = 'green'
        ax3_2.set_ylabel('MLT', color=color,fontsize=8)  # we already handled the x-label with ax1
        ax3_2.tick_params(axis='y', colors=color)
        ax3_2.plot(np.array([time_conversions.h5_time_conversion(x) for x in lanl_d['IsoTime']]), 
            MLT,color=color)
        ax3_2.plot(lanl_times, abs(MLAT_south),color='red',ls='--',label=r'$\lambda_m\ <\ 0.$')
        ax3_2.plot(lanl_times, abs(MLAT_north),color='red',label=r'$\lambda_m\ >\ 0.$')
        ax3_2.legend(loc=(1.1, 0.35),labelcolor='red')
        ax3_2.set_ylim(0,24)
        ax3_2.set_xlim(left=np.array([time_conversions.h5_time_conversion(x) for x in lanl_d['IsoTime']])[lstar >= 0][0],right=np.array([time_conversions.h5_time_conversion(x) for x in lanl_d['IsoTime']])[lstar >= 0][-1])


        cax_list[plot_number].axis('off')

        return vertical_change2
    
    def geo_indicies(epoch_omni, AE,
                     burst_times,lanl_d,
                     plot_number,axs_list,cax_list):
        
        lstar = np.array(lanl_d['Lstar'][:, 0])

        axs_list[plot_number].plot(epoch_omni, AE,color='blue')
        axs_list[plot_number].set_ylabel('AE')
        axs_list[plot_number].set_xlabel('UTC')
        axs_list[plot_number].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        hours = mdates.HourLocator(interval=1)
        axs_list[plot_number].xaxis.set_minor_locator(hours)
        axs_list[plot_number].set_xlim(left=np.array([time_conversions.h5_time_conversion(x) for x in lanl_d['IsoTime']])[lstar >= 0][0],right=np.array([time_conversions.h5_time_conversion(x) for x in lanl_d['IsoTime']])[lstar >= 0][-1])
        vertical_change3=axs_list[plot_number].axvline(burst_times[0],ymin=0,ymax=10**4,linestyle='dashed',linewidth=2,color='black')

        cax_list[plot_number].axis('off')

        return vertical_change3
    

    def dist_large(burst_dist,
                   survey_int,burst_int,
                   plot_number,axs_list,cax_list):

        n_b = len(burst_dist[0])
        T = 1./35000.
        win = 16384
        T_win = T*win
        nt = np.linspace(T_win/2., T_high-T_win/2.,n_b)
       
        axs_list[plot_number].set_ylabel(r'$Power\ (nT^2)$',fontsize=8)
        axs_list[plot_number].set_xlabel(r'$Time\ (s)$',fontsize=8)

        b_dist = axs_list[plot_number].plot(nt, burst_dist[0], color = "r",linestyle = 'dotted', marker = '*')[0]
        dist_mean= np.mean(burst_dist[0])
        dist_sd = np.std(burst_dist[0])
        print(dist_mean,dist_sd)
        mean_12 = axs_list[plot_number].axhline(dist_mean,0,5.968,color = 'red',label='Mean of windows')
        #ax5.fill_between(x=[0.,5.968], y1=(dist_mean-dist_sd**2), y2=(dist_mean+dist_sd**2), color='red', alpha=0.2)
        survey_i = axs_list[plot_number].axhline(survey_int[0],0,5.968,color='black', label="Survey power")
        b_box=axs_list[plot_number].axhline(burst_int[0],0,5.968,color='blue',label="0.5s box --> survey")
        
        cax_list[plot_number].axis('off')
        #b_mean=ax5.hlines(burst_mean[0],np.min(nt_array),np.max(nt_array),color='green',linestyles='dotted', label="Full burst duration mean power")
        #b_max=ax5.set_ylim(0,np.max(burst_int[0])+0.02)
        #ax5.add_patch(Rectangle((box_liml,0.),0.5,0.175,
                        # edgecolor='red',
                            #facecolor='none',
                            #lw=4))
        axs_list[plot_number].legend()

        return b_dist,mean_12,survey_i,b_box,nt


    values =np.zeros(2*len(burst_dist[0]))
    value_times = np.zeros(2*len(burst_dist[0]))

    T = 1./35000.
    win = 16384
    T_win = T*win
    times = np.linspace(T_win,T_high,12)

    value_times[0] = 0 
    value_times[-1] = T_high

    for i in range(12):
        values[2*i] = burst_dist[0][i]
        values[2*i + 1] = burst_dist[0][i]

        if i<11:
            value_times[2*i + 1] = times[i]
            value_times[2*i + 2] = times[i]

    
    
    def dist_small(highB_dist,survey_int,
                values,value_times,
                plot_number,axs_list,cax_list):
        
        

        n_hb = len(highB_dist[0])

        # creating time aray 
        
        #nth = np.linspace(0.03/2., 5.968-0.03/2.,n_hb)
        axs_list[plot_number].set_ylabel(r'$Power\ (nT^2)$',fontsize=8)
        axs_list[plot_number].set_xlabel(r'$Time\ (s)$',fontsize=8)
        axs_list[plot_number].set_yscale('log')
        axs_list[plot_number].set_xlim(0,6)
        
        hb_dist = axs_list[plot_number].plot(nt_array, highB_dist[0], color = "r",linestyle = 'dotted', marker = '*')[0]
        #dist_mean= np.mean(highB_dist[0])
        #dist_sd = np.std(highB_dist[0])
        dist_median = np.median(highB_dist[0])
        #print('high_b',dist_mean,dist_sd)
        #mean_high = axs_list[plot_number].axhline(dist_mean,0,5.968,color = 'red',label='Mean of windows')
        #hb_median = axs_list[plot_number].axhline(dist_median,0,5.968,color = 'blue',label='median of windows')
        #ax5.fill_between(x=[0.,5.968], y1=(dist_mean-dist_sd**2), y2=(dist_mean+dist_sd**2), color='red', alpha=0.2)
        survey_ih = axs_list[plot_number].axhline(survey_int[0],0,5.968,color='black', label="Survey power")
        hside=axs_list[plot_number].plot(value_times,values,label='0.468s window powers')[0]
        axs_list[plot_number].legend()

        cax_list[plot_number].axis('off')

        return hb_dist,hside,survey_ih,nt_array
    

    def spectrum_comparison(av_specs,av_specs_half,
                            freq_survey,surv_set,
                            plot_number,axs_list,cax_list):
        
        av_specs_half = axs_list[plot_number].plot(freq_survey,av_specs_half[0],color='black',label='Averaged spectrum over first 0.468s')
        av_spec = axs_list[plot_number].plot(freq_survey,av_specs[0],color='black',linestyle ='dashed',label='Averaged spectrum over 6s')
       
        axs_list[plot_number].plot(freq_survey,surv_set[0][274,:],label = 'Survey spectrum')
        axs_list[plot_number].set_yscale('log')
        axs_list[plot_number].legend()
        axs_list[plot_number].set_xlabel('Frequency (Hz)',fontsize=8)
        axs_list[plot_number].set_ylabel(r'$PSD (nT^2/Hz)$',fontsize=8)
        cax_list[plot_number].axis('off')
        
        return av_spec,av_specs_half

    # set the size of the summary plot
    plt.gcf().set_size_inches((12, 24))

    # Making plot timestamps into form without milliseconds - there defo exists an easier way
    timestamp_txt=[]
    for time_stamp in burst_times:
        time_stamp=datetime.strftime(time_stamp, '%Y-%m-%d %H:%M:%S')
        timestamp_txt.append(datetime.strptime(time_stamp, '%Y-%m-%d %H:%M:%S'))

    # What order do I want the plots in:
    
    # 1. survey
    p_num = 0 #plot number

    survey_img, vertical_change_survey = survey_plot(epoch_survey, freq_survey,
                                            Btotal, burst_times, 
                                             p_num,axs_list,cax_list)
    

    

    # 2. spacecraft location
    p_num = 1
    vertical_change_loc = space_locations(lanl_d,burst_times, 
                        p_num,axs_list,cax_list)
    

    # 3. geo indicies
    p_num = 2
    vertical_change_gi = geo_indicies(epoch_omni, AE,
                                     burst_times,lanl_d,
                                     p_num,axs_list,cax_list)


    gyro_lines = np.zeros((2,2))
    gyro_text = np.zeros((2,2))
    # 4. small burst windows 
    p_num = 4

    burst_small,gyro_line1_small,gyro_line2_small,text_1_small,text_2_small = burst_small(freq_set,data_set,
                                                                gyro_1,gyro_low_set,
                                                                p_num,axs_list,cax_list)

    
    
    # 5. larger burst windows
    p_num = 3
    burst_large,gyro_line1_large,gyro_line2_large,text_1_large,text_2_large = burst_large(freq_468, mag_lists,
                                                                gyro_1, gyro_low_set,
                                                                burst_dist,
                                                                p_num,axs_list,cax_list)

   


    # 6. small burst distribution 
    dist_mean = np.zeros(2)
    

    p_num = 6

    hb_dist,hside,survey_ih,small_t_arr = dist_small(highB_dist,survey_int,
                                                     values,value_times,
                                                     p_num,axs_list,cax_list,)


    

    # 7. large burst distribution

    p_num = 5
    b_dist,mean_12,survey_i,b_box,large_t_arr = dist_large(burst_dist,
                                                        survey_int,burst_int,
                                                        p_num,axs_list,cax_list)
    

    
    

    # 8. frequency spectrum

    p_num = 7

    av_spec ,av_spec_half= spectrum_comparison(av_specs,av_specs_half,
                            freq_survey,surv_set,
                            p_num,axs_list,cax_list)
    

    
   
    
    # Animate the summary plot

    def animate(i):
        
        burst_small.set_array(data_set[i].T.flatten())

        vertical_change_survey.set_xdata(burst_times[i])
        gyro_line1_small.set_ydata(gyro_half_set[i])
        gyro_line2_small.set_ydata(gyro_low_set[i])
        vertical_change_gi.set_xdata(burst_times[i])
        vertical_change_loc.set_xdata(burst_times[i])
        text_1_small.set_position((4.5, gyro_half_set[i]+200))
        text_2_small.set_position((4.5, gyro_low_set[i]+200))
        #int_plot.set_ydata(burst_int[i])
        #b_mean.set_ydata(burst_mean[i])
        b_box.set_ydata(burst_int[i])
        survey_i.set_ydata(survey_int[i])

        gyro_line1_large.set_ydata(gyro_half_set[i])
        gyro_line2_large.set_ydata(gyro_low_set[i])
        text_1_large.set_position((4.5, gyro_half_set[i]+200))
        text_2_large.set_position((4.5, gyro_low_set[i]+200))

        values =np.zeros(2*len(burst_dist[0]))
        value_times = np.zeros(2*len(burst_dist[0]))
        
        T = 1./35000.
        win = 16384
        T_win = T*win
        times = np.linspace(T_win,T_high,12)
        

        value_times[0] = 0 
        value_times[-1] = 5.968

        for j in range(12):
            values[2*j] = burst_dist[i][j]
            values[2*j + 1] = burst_dist[i][j]

            if j<11:
                value_times[2*j + 1] = times[j]
                value_times[2*j + 2] = times[j]
        
        b_dist.set_xdata(large_t_arr)
        b_dist.set_ydata(burst_dist[i])
        dist_mean,dist_sd = np.mean(burst_dist[i]),np.std(burst_dist[i])
        
        mean_12.set_ydata(dist_mean)

        #ax5.fill_between(x=[0.,5.968], y1=dist_mean-dist_sd**2, y2=dist_mean+dist_sd**2, color='red', alpha=0.2)
        #ax5.set_ylim(np.min(burst_dist[i])-np.min(burst_dist[i])/10, np.max(burst_dist[i])+np.max(burst_dist[i])/10)

        burst_large.set_array(mag_lists[i].T.flatten())
        
        
        #b_dist = ax5.scatter(nt,burst_dist[0],color='red')
        hb_dist.set_xdata(small_t_arr)
        hb_dist.set_ydata(highB_dist[i])
        dist_mean= np.mean(highB_dist[i])
        dist_sd = np.std(highB_dist[i])
        dist_median = np.median(highB_dist[i])
        hside.set_xdata(value_times)
        hside.set_ydata(values)
        
        #ax7.set_ylim(np.min(highB_dist[i])-np.min(highB_dist[i])/10, np.max(highB_dist[i])+np.max(highB_dist[i])/10)
        
        #av_spec.set_ydata(av_specs[i])
    

        # Add timestamp
        fig.suptitle(f'Burst timestamp at: {timestamp_txt[i]} UTC', fontsize=14)
        #add_timestamp(ax1, timestamp_txt[i], y=0.9,x=0.3, time_format = '%Y-%m-%d %H:%M:%S', high_contrast=True)
        plt.savefig(str(i))

    anim = animation.FuncAnimation(fig, animate, interval=600, frames= len(burst_times)-1)


    file_path = '/data/emfisis_burst/wip/rablack75/rablack75/simple_FFT/burst/'+year
    os.makedirs(file_path, exist_ok=True)
    os.makedirs(file_path+'/'+month, exist_ok=True)
    
    anim.save(file_path+'/'+month+'/' + 'plot_'+str(day)+'.gif')
    print(str(day),'saved')
