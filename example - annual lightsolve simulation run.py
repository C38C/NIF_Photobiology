"""
This code provides an example of an annual simulation for 4 views within a series of ALFA simulations as seen in:
    Evaluating the use of photobiology-driven alertness and health measures for circadian lighting design
    A Alight, JA Jakubiec
    http://jakubiec.net/papers/Alight%20and%20Jakubiec%20-%20Comparison%20of%20photobiological%20measures%20(Final).pdf
    Proceedings of BS2021

    Spectral and biological simulation methods for the design of healthy circadian lighting
    http://www.jakubiec.net/papers/Jakubiec%20and%20Alight%20-%20NIF%20Lighting%20Metrics%20(Final).pdf
    JA Jakubiec, A Alight
    Proceedings of BS2021
"""

# custom classes
import biological_model
import spectral_conversion
import helper_functions as hf

# standard Python classes
import datetime, statistics, os
import numpy as np
import math


views = [
    ["L1", 1],
    ["L2", 9],
    ["L3", 17] ] # Just a list of view names / indices from a larger ALFA simulation file with many views

    


# output directories
dir_30m = '.\\30m\\'
dir_1m = '.\\1m\\'
dir_day = '.\\day\\'

subfolder = "daylight_electric_screen"
label_add = ""

# create directories if they don't exist
if os.path.isdir(dir_30m + subfolder):
    pass
else:
    try:
        os.mkdir(dir_30m)
    except:
        pass
    os.mkdir(dir_30m + subfolder)

if os.path.isdir(dir_1m + subfolder):
    pass
else:
    try:
        os.mkdir(dir_1m)
    except:
        pass
    os.mkdir(dir_1m + subfolder)

if os.path.isdir(dir_day + subfolder):
    pass
else:
    try:
        os.mkdir(dir_day)
    except:
        pass
    os.mkdir(dir_day + subfolder)
    
dir_30m = dir_30m + subfolder + '\\'
dir_1m = dir_1m + subfolder + '\\'
dir_day = dir_day + subfolder + '\\'


fname_pre = "annual-output" # output filenames prefix

directory = '.\\example - ALFA daylight results\\' # directory containing ALFA simulation files

base_timestep = 20 # s
n_timesteps_per_day = 24 * 60 * 60 / base_timestep # timesteps per day


# ######################################
# Load in electric light sources #######
# ######################################
CCT_4100 = spectral_conversion.read_alfa_file('.\\example - ALFA electric light results\\low.csv')
CCT_6100 = spectral_conversion.read_alfa_file('.\\example - ALFA electric light results\\neutral.csv')
CCT_16000 = spectral_conversion.read_alfa_file('.\\example - ALFA electric light results\\high.csv')

# lists to keep the Eemel values per view
Eemel_4100CCT_view = []
Eemel_6100CCT_view = []
Eemel_16000CCT_view = []

for viewname,viewindex in views: # get electric lighting data / 4100 K 
    n_view = 0
    for sd in CCT_4100: # get 4100K irrad from ALFA simulations
        if (sd.type == 'view') and (n_view == viewindex): # Look for views of interest
            ill, mel_ill, m_over_p, opic_irrads, opic_elrs, opic_ders, opic_edis = spectral_conversion.calc_s026(sd.get_data())
            this_s026 = spectral_conversion.s026_metrics(ill, mel_ill, m_over_p, opic_irrads, opic_elrs, opic_ders, opic_edis)
            Eemel = opic_irrads["mel"]
        if sd.type == 'view':
            n_view += 1
    Eemel_4100CCT_view.append(Eemel) # append Ee,mel melanopic irradiance

for viewname,viewindex in views: # get electric lighting data / 6100 K 
    n_view = 0
    for sd in CCT_6100: # get 6100K irrad from ALFA simulations
        if (sd.type == 'view') and (n_view == viewindex): # Look for views of interest
            ill, mel_ill, m_over_p, opic_irrads, opic_elrs, opic_ders, opic_edis = spectral_conversion.calc_s026(sd.get_data())
            this_s026 = spectral_conversion.s026_metrics(ill, mel_ill, m_over_p, opic_irrads, opic_elrs, opic_ders, opic_edis)
            Eemel = opic_irrads["mel"]
        if sd.type == 'view':
            n_view += 1
    Eemel_6100CCT_view.append(Eemel) # append Ee,mel melanopic irradiance

for viewname,viewindex in views: # get electric lighting data / 16000 K 
    n_view = 0
    for sd in CCT_16000: # get 16000K irrad from ALFA simulations
        if (sd.type == 'view') and (n_view == viewindex): # Look for views of interest
            ill, mel_ill, m_over_p, opic_irrads, opic_elrs, opic_ders, opic_edis = spectral_conversion.calc_s026(sd.get_data())
            this_s026 = spectral_conversion.s026_metrics(ill, mel_ill, m_over_p, opic_irrads, opic_elrs, opic_ders, opic_edis)
            Eemel = opic_irrads["mel"]
        if sd.type == 'view':
            n_view += 1
    Eemel_16000CCT_view.append(Eemel) # append Ee,mel melanopic irradiance


# ##########################################
# Monitor Ee,mel sources from f.luxometer ##
# ##########################################
# https://fluxometer.com/rainbow/#!id=NEC%20PA272W%20sRGB/1900K-NEC%20PA272W%20sRGB
E_mon_6500 = 0.113 # W/m2 Ee,mel
E_mon_1900 = 0.0183 # W/m2 Ee,mel


#
# get info for simulated timesteps from lightsolve timeseries files (see helper_functions.sunposition_climate()).
#
n_days = 8
n_hours_per_day = 7

day_files = []
day_times = []
day_irrad = []

# read lightsolve timestep data
with open('.\\example - Toronto-ON-CAN_solar-times_irrad_skytype.csv', 'r') as f:
    f.readline() # skip header
    
    for n_day in range(n_days):
        this_day_files = []
        this_day_times = []
        this_day_irrad = []
        for n_step in range(n_hours_per_day):
            line = f.readline()
            
            parts = line.rstrip().split(',')
            mo_str = parts[0]
            da_str = parts[1]
            time_str = parts[2]
            hr_str = parts[3]
            mn_str = parts[4]
            alt = parts[5]
            azi = parts[6]
            irr_str = parts[7]
            
            this_fname = '%.1f_%.1f_%.2f.csv' % (int(mo_str), int(da_str), float(time_str))
            this_irrad = float(irr_str)
            this_datetime = datetime.datetime(2020, int(mo_str), int(da_str), int(hr_str), int(mn_str))
            
            this_day_files.append(this_fname)
            this_day_times.append(this_datetime)
            this_day_irrad.append(this_irrad)
        
        day_files.append(this_day_files)
        day_times.append(this_day_times)
        day_irrad.append(this_day_irrad)

def interp_ndays(I_a, I_b, n_dayrange):
    """Data interpolation function."""
    data = [[0] * len(I_a) for n in range(n_dayrange)]
    
    for h in range(len(I_a)): # loop through hour
        I_h_daygap = [I_a[h], I_b[h]]
        x_val = np.array([0,n_dayrange-1])
        x_interp = np.linspace(0, n_dayrange-1, n_dayrange) 
        
        I_h_dayrange = list( np.interp(x_interp, x_val, I_h_daygap) )
        
        # fill in data from I_h_dayrange
        for d in range(n_dayrange):
            data[d][h] = I_h_dayrange[d]
    
    return data
        
def day_index(month, day):
    """Returns 0-based day index 0-364."""
    index = 0
    
    days_in_mo = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    for m in range(month - 1):
        index += days_in_mo[m]
    
    index += day - 1
    
    return index

class median_metrics:
    """
    A custom class to keep track of interesting circadian metrics. See Postnova, et al. 2018 for more information.
    """
    
    def get_header_str(self):
        return "alpha,I,theta_L_KSS,S,F_w,V_m,H,C,D_v,p_b,p_b_max,Q_v,Q_m,W,vPVTL,KSS,vPVTRT_mean,deltaADD,deltaDSST"
    
    def get_data_str(self):
        strings = []
        # strings.append('%i' % self.day)
        # strings.append('%.6f' % self.clocktime)
        strings.append('%.6f' % self.alpha)
        strings.append('%.6f' % self.I)
        strings.append('%.6f' % self.theta_L_KSS)
        strings.append('%.6f' % self.S)
        strings.append('%.6f' % self.F_w)
        strings.append('%.6f' % self.V_m)
        strings.append('%.6f' % self.H)
        strings.append('%.6f' % self.C)
        strings.append('%.6f' % self.D_v)
        strings.append('%.6f' % self.p_b)
        strings.append('%.6f' % self.p_b_max)
        strings.append('%.6f' % self.Q_v)
        strings.append('%.6f' % self.Q_m)
        strings.append('%.6f' % self.W)
        strings.append('%.6f' % self.vPVTL)
        strings.append('%.6f' % self.KSS)
        strings.append('%.6f' % self.vPVTRT_mean)
        strings.append('%.6f' % self.deltaADD)
        strings.append('%.6f' % self.deltaDSST)
        
        return ','.join(strings)
    
    def __init__(self, list_state, list_metrics):
        self.day = statistics.mean([s.day for s in list_state])
        
        clocktime = []
        for s in list_state:
            clocktime.append( s.n / (60.0 * 60.0 / s.dt) - (s.day * 24.0))
        
        self.clocktime = statistics.mean(clocktime)
        
        self.alpha = statistics.median([s.alpha for s in list_state])
        self.vPVTL = statistics.median([s.vPVTL for s in list_metrics])
        self.KSS = statistics.median([s.KSS for s in list_metrics])
        self.vPVTRT_mean = statistics.median([s.vPVTRT_mean for s in list_metrics])
        self.deltaADD = statistics.median([s.deltaADD for s in list_metrics])
        self.deltaDSST = statistics.median([s.deltaDSST for s in list_metrics])
        self.p_b = statistics.median([s.p_b for s in list_state])
        self.p_b_max = statistics.median([s.p_b_max for s in list_state])
        self.W = statistics.median([s.W for s in list_state])
        self.S = statistics.median([s.S for s in list_state])
        self.F_w = statistics.median([s.F_w for s in list_state])
        self.V_m = statistics.median([s.V_m for s in list_state])
        self.Q_v = statistics.median([s.Q_v for s in list_state])
        self.Q_m = statistics.median([s.Q_m for s in list_state])
        self.H = statistics.median([s.H for s in list_state])
        self.C = statistics.median([s.C for s in list_state])
        self.D_v = statistics.median([s.D_v for s in list_state])
        self.I = statistics.median([s.I for s in list_state])
        self.theta_L_KSS = statistics.median([s.theta_L_KSS for s in list_state])
        
for n_view, (viewname, viewindex) in enumerate(views):
    # 
    # calculate each day
    # 
    day_eemel_24 = [] # contains 24 h of melanopic irrad for each day and the current view

    for n_day in range(n_days): # loop through days
        datetime_list = day_times[n_day]
        daylight_files = day_files[n_day]
        daylight_irrad = day_irrad[n_day]

        # calculate sunrise and sunset (although we have it in the csv file for calculations)
        simulation_timestep = datetime_list[1] - datetime_list[0] # time difference in the format datetime.timedelta
        sunrise_dt = datetime_list[0] - simulation_timestep / 2
        sunset_dt = datetime_list[-1] + simulation_timestep / 2
        
        daylight_info = [] # list of lists of tuples 
        #   # (view)[ (time)[ tuple(spectral_data, s026_metrics) ] ]

        # load spectral data from all files 
        for n_file, daylight_file in enumerate(daylight_files):
            file_spectral_data = spectral_conversion.read_alfa_file(directory + daylight_file)
            
            n_reading_view = 0
            for sd in file_spectral_data:
                
                if sd.type == 'view' and (n_reading_view == viewindex):
                    if n_file == 0 and n_day == 0:
                        print('Calculating view %s at location %.2f,%.2f,%.2f and direction %.2f,%.2f,%.2f.' % (viewname, sd.X, sd.Y, sd.Z, sd.dX, sd.dY, sd.dZ))
                    # convert to S026 metrics and class-ify
                    ill, mel_ill, m_over_p, opic_irrads, opic_elrs, opic_ders, opic_edis = spectral_conversion.calc_s026(sd.get_data())
                    this_s026 = spectral_conversion.s026_metrics(ill, mel_ill, m_over_p, opic_irrads, opic_elrs, opic_ders, opic_edis)
                    
                    daylight_info.append( (sd, this_s026) ) # append hour result to the day list
                if sd.type == 'view':
                    n_reading_view += 1

        
        E_list_sim = []
        for n_timestep in range( len(daylight_info) ):
            E_list_sim.append( daylight_info[n_timestep][1].opic_irrads["mel"] ) # messy, but this is  dailytimestep, s026_metrics, .opic_irrads["mel"] 
            
        # create interpolated irradiance for the biological model
        E_list_24h_daylight = hf.day_from_ill_datetime(E_list_sim, datetime_list, sunrise_dt, sunset_dt, timestep = base_timestep)
        
        day_eemel_24.append(E_list_24h_daylight)
    
    
    # get view coords & direction 
    # daylight_info[n_view_use].append( (sd, this_s026) )
    view_X = daylight_info[0][0].X 
    view_Y = daylight_info[0][0].Y 
    view_Z = daylight_info[0][0].Z 
    view_dX = daylight_info[0][0].dX 
    view_dY = daylight_info[0][0].dY 
    view_dZ = daylight_info[0][0].dZ 
    
    
    
    #
    # setup empty array
    #
    eemel_8760 = []
    eemel_8760_electric = []
    eemel_8760_screen = []
    for d in range(365):
        eemel_8760.append([])
        eemel_8760_electric.append([])
        eemel_8760_screen.append([])
    

    #
    # interpolate between end days of the annual loaf
    #
    total_dayrange = day_times[-1][0] - day_times[0][0]
    dayrange_loafends = 365 - total_dayrange.days 
    loafend_data = interp_ndays(day_eemel_24[-1], day_eemel_24[0], dayrange_loafends)

    last_day_index = day_index(day_times[-1][0].month, day_times[-1][0].day)
    first_day_index = day_index(day_times[0][0].month, day_times[0][0].day)
    for data_index, index in enumerate(range(last_day_index, 365, 1)):
        eemel_8760[index].extend( loafend_data[data_index][:] )

    for n_day, index in enumerate(range(0, first_day_index + 1, 1)):
        data_index = n_day + (364 - last_day_index)
        eemel_8760[index].extend( loafend_data[data_index][:] )

    #
    # interpolate between contiguous days of the year 
    #
    for d_0 in range(n_days - 1):
        d_1 = d_0 + 1
        
        # this_dayrange = day_times[d_1][0] - day_times[d_0][0]
        # this_dayrange = this_dayrange.days
        
        day_index_start = day_index(day_times[d_0][0].month, day_times[d_0][0].day)
        day_index_end = day_index(day_times[d_1][0].month, day_times[d_1][0].day)
        
        this_dayrange = day_index_end - day_index_start
        
        this_data = interp_ndays(day_eemel_24[d_0], day_eemel_24[d_1], this_dayrange)
        
        for n_day, index in enumerate(range(day_index_start, day_index_end, 1)):
            # print index, n
            eemel_8760[index].extend( this_data[n_day][:] )
    
    for n_day in range(len(eemel_8760)):
        if len(eemel_8760[n_day]) > n_timesteps_per_day:
            eemel_8760[n_day] = eemel_8760[n_day][:n_timesteps_per_day]
    
    # add in electric lighting and monitor lighting
    day_start_index = 6 * 60 * 60 / base_timestep
    day_end_index = len(eemel_8760[0])
    
    
    
    for index in range(365):
        eemel_8760_electric[index].extend( [0.0] *  (6 * 60 * 60 / base_timestep)) 
        eemel_8760_electric[index].extend( [Eemel_6100CCT_view[n_view]] *  (12 * 60 * 60 / base_timestep)) 
        eemel_8760_electric[index].extend( [Eemel_4100CCT_view[n_view]] *  (4 * 60 * 60 / base_timestep)) 
        eemel_8760_electric[index].extend( [Eemel_4100CCT_view[n_view]] *  (2 * 60 * 60 / base_timestep)) 
        
        eemel_8760_screen[index].extend( [0.0] * (18 * 60 * 60 / base_timestep)) 
        eemel_8760_screen[index].extend( [E_mon_6500] * (3 * 60 * 60 / base_timestep)) 
        eemel_8760_screen[index].extend( [E_mon_1900] * (3 * 60 * 60 / base_timestep)) 

    min30_metrics = [] # [view][time]
    min1_metrics = [] # [view][time]
    daily_metrics = []
    daily_KSS = []
    daily_vPVTRT = []

    # entrain a subject
    person = biological_model.subject(timestep=base_timestep)

    for d_sim in range(365):
        E_daylight = eemel_8760[d_sim]    
        E_electric = eemel_8760_electric[d_sim]    
        E_screen = eemel_8760_screen[d_sim]    
        
        E_total = []
        for i in range(len(E_daylight)):
            E_total.append(E_daylight[i] + E_electric[i] + E_screen[i])
        
        person.simulate_day(E_total)

    person.generate_alertness_metrics()
    
    # get daily metrics
    daily_metrics.extend(person.day_metrics)

    X = view_X
    Y = view_Y
    Z = view_Z
    dX = view_dX
    dY = view_dY
    dZ = view_dZ
        
    # get median metrics every 30 minutes
    increment = 30 * 60 / base_timestep
    for n in range(0, 365 * n_timesteps_per_day, increment):
        # if n_view == (len(daylight_info) - 1):
        #     print (n, n+increment)
        min30_metrics.append(median_metrics(person.history[n : n + increment],person.metrics[n : n + increment] ))
    
    # get median metrics every 1 minute
    increment = 1 * 60 / base_timestep
    for n in range(0, 365 *  n_timesteps_per_day, increment):
        # if n_view == (len(daylight_info) - 1):
        #     print (n, n+increment)
        min1_metrics.append(median_metrics(person.history[n : n + increment],person.metrics[n : n + increment] ))

    with open(dir_30m + fname_pre + '_%s%s.csv' % (viewname, label_add), 'w') as f:
        f.write('X,Y,Z,dX,dY,dZ,day,clocktime')
        f.write(',%s' % min30_metrics[0].get_header_str())
        f.write('\n')
        
        for n_time in range(len(min30_metrics)):
            f.write("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f," % (X,Y,Z,dX,dY,dZ))
            f.write('%i,%.3f' % (min30_metrics[n_time].day, min30_metrics[n_time].clocktime) )
            f.write(',%s' % min30_metrics[n_time].get_data_str() )
            f.write('\n')
    
    with open(dir_30m + fname_pre + '_%s%s_LIGHT.csv' % (viewname, label_add), 'w') as f:
        # get 30-minute lighting measures
        increment = 30 * 60 / base_timestep
        f.write('i,day,clocktime,I,Source\n')
        for d in range(365):
            for n in range(0, n_timesteps_per_day, increment):
                f.write('%i.%i,%i,%.2f,%.12f,%s\n' % (d, n, d, (float(n) * base_timestep) / (60.0 * 60.0), eemel_8760[d][n], "Daylight"))
                f.write('%i.%i,%i,%.2f,%.12f,%s\n' % (d, n, d, (float(n) * base_timestep) / (60.0 * 60.0), eemel_8760_electric[d][n], "Electric"))
                f.write('%i.%i,%i,%.2f,%.12f,%s\n' % (d, n, d, (float(n) * base_timestep) / (60.0 * 60.0), eemel_8760_screen[d][n], "Screen"))
    
    with open(dir_1m + fname_pre + '_%s%s_LIGHT.csv' % (viewname, label_add), 'w') as f:
        # get 1-minute lighting measures
        increment = 1 * 60 / base_timestep
        f.write('i,day,clocktime,I,Source\n')
        for d in range(365):
            for n in range(0, n_timesteps_per_day, increment):
                f.write('%i.%i,%i,%.2f,%.12f,%s\n' % (d, n, d, (float(n) * base_timestep) / (60.0 * 60.0), eemel_8760[d][n], "Daylight"))
                f.write('%i.%i,%i,%.2f,%.12f,%s\n' % (d, n, d, (float(n) * base_timestep) / (60.0 * 60.0), eemel_8760_electric[d][n], "Electric"))
                f.write('%i.%i,%i,%.2f,%.12f,%s\n' % (d, n, d, (float(n) * base_timestep) / (60.0 * 60.0), eemel_8760_screen[d][n], "Screen"))
    
    
    with open(dir_1m + fname_pre + '_%s%s.csv' % (viewname, label_add), 'w') as f:
        f.write('X,Y,Z,dX,dY,dZ,day,clocktime')
        f.write(',%s' % min1_metrics[0].get_header_str())
        f.write('\n')
        
        for n_time in range(len(min1_metrics)):
            f.write("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f," % (X,Y,Z,dX,dY,dZ))
            f.write('%i,%.3f' % (min1_metrics[n_time].day, min1_metrics[n_time].clocktime) )
            f.write(',%s' % min1_metrics[n_time].get_data_str() )
            f.write('\n')

    with open(dir_day + fname_pre + '_%s%s.csv' % (viewname,label_add), 'w') as f:
        metric_list = daily_metrics[0].header_str(day=False).split(',')
        metric_list.append('KSS_8-12')
        metric_list.append('KSS_12-18')
        metric_list.append('KSS_18-24')
        metric_list.append('KSS_day')
        metric_list.append('vPVTRT_8-12')
        metric_list.append('vPVTRT_12-18')
        metric_list.append('vPVTRT_18-24')
        metric_list.append('vPVTRT_day')
        metric_list.append('Eemel_8-12')
        metric_list.append('Eemel_12-18')
        metric_list.append('Eemel_18-24')
        metric_list.append('Eemel_day')
        
        f.write('X,Y,Z,dX,dY,dZ,')
        f.write('day,' + ','.join(metric_list) + '\n')
        
        for n_day in range(len(daily_metrics)):
            index_30min_morning = [i + n_day * 48 for i in range(16,24,1)] # 8:00 - 12:00
            index_30min_midday = [i + n_day * 48 for i in range(24,36,1)] # 12:00 - 18:00
            index_30min_evening = [i + n_day * 48 for i in range(36,48,1)] # 18:00 - 24:00
            index_30min_day = [i + n_day * 48 for i in range(16,48,1)] # 8:00 - 24:00
            
            f.write("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f," % (X,Y,Z,dX,dY,dZ))
            f.write('%i,%s' % (n_day, daily_metrics[n_day].output_str(day=False) ) )
            
            KSS_morning = []
            KSS_midday = []
            KSS_evening = []
            KSS_day = []
            vPVTRT_morning = []
            vPVTRT_midday = []
            vPVTRT_evening = []
            vPVTRT_day = []
            
            Eemel_morning = []
            Eemel_midday = []
            Eemel_evening = []
            Eemel_day = []
            
            for i in index_30min_morning:
                KSS_morning.append(min30_metrics[i].KSS)
                vPVTRT_morning.append(min30_metrics[i].vPVTRT_mean)
                Eemel_morning.append(min30_metrics[i].I)
            for i in index_30min_midday:
                KSS_midday.append(min30_metrics[i].KSS)
                vPVTRT_midday.append(min30_metrics[i].vPVTRT_mean)
                Eemel_midday.append(min30_metrics[i].I)
            for i in index_30min_evening:
                KSS_evening.append(min30_metrics[i].KSS)
                vPVTRT_evening.append(min30_metrics[i].vPVTRT_mean)
                Eemel_evening.append(min30_metrics[i].I)
            for i in index_30min_day:
                KSS_day.append(min30_metrics[i].KSS)
                vPVTRT_day.append(min30_metrics[i].vPVTRT_mean)
                Eemel_day.append(min30_metrics[i].I)
            
            f.write(',%.2f' % statistics.mean(KSS_morning))
            f.write(',%.2f' % statistics.mean(KSS_midday))
            f.write(',%.2f' % statistics.mean(KSS_evening))
            f.write(',%.2f' % statistics.mean(KSS_day))
            f.write(',%.2f' % statistics.mean(vPVTRT_morning))
            f.write(',%.2f' % statistics.mean(vPVTRT_midday))
            f.write(',%.2f' % statistics.mean(vPVTRT_evening))
            f.write(',%.2f' % statistics.mean(vPVTRT_day))
            f.write(',%.2f' % statistics.mean(Eemel_morning))
            f.write(',%.2f' % statistics.mean(Eemel_midday))
            f.write(',%.2f' % statistics.mean(Eemel_evening))
            f.write(',%.2f' % statistics.mean(Eemel_day))
            f.write('\n')