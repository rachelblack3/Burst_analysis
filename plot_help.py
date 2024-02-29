import matplotlib.pyplot as plt
fig, axs = plt.subplots(2,2,gridspec_kw={"width_ratios":(1,0.05)})

# Name plot axis
axs_list = axs[:, 0]

# Name colourbar axis (the third/fourth one will be made invisible)
cax1,cax2=axs[:,1]

# 1. Survey

def survey_plot(axs_list,plot_number):
    import numpy as np
    

    epoch_survey = np.linspace(0,10,100)
    freq_survey = np.linspace(0,5,100)

    survey_img = axs_list[plot_number].plot(epoch_survey, freq_survey)
    
    
    axs_list[plot_number].set_ylabel('Frequency (Hz)')

    return survey_img


surv_img = survey_plot(axs_list,1)

plt.savefig('plot_helpp')