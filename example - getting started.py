"""
This code provides a getting started example of using the Postnova, et al. photobiological dynamics model as illustrated in:
    Evaluating the use of photobiology-driven alertness and health measures for circadian lighting design
    A Alight, JA Jakubiec
    http://jakubiec.net/papers/Alight%20and%20Jakubiec%20-%20Comparison%20of%20photobiological%20measures%20(Final).pdf
    Proceedings of BS2021

    Spectral and biological simulation methods for the design of healthy circadian lighting
    http://www.jakubiec.net/papers/Jakubiec%20and%20Alight%20-%20NIF%20Lighting%20Metrics%20(Final).pdf
    JA Jakubiec, A Alight
    Proceedings of BS2021
"""

import biological_model
import spectral_conversion
import statistics # pip install statistics


# You can entrain a new subject as follows...
person = biological_model.subject(timestep=20) # 20 second timestep for dynamic photobiology calculations
#       By default, entrain is set to True, and a subject is entrained to 25 days of sleep from midnight to 8 am / 0.2 W/m2 Ee,mel from 8 am to 8 pm and 0.03 W/m2 Ee,mel from 8 pm until midnight (when sleep occurs)

# Melanopic irradiance can be provided on a 24-hour daily basis and simulated per-day along with sleep / wake schedules
I = [0.0] * 6 + [1.0] * 12 + [0.2] * 6 # W/m2 Ee,mel melanopic irradince is provided as a list. 
#       The data is interpolated across a 24-hour period at the timestep interval (default 20s). 
#       In this case, 24 hours of data is provided as the input, but 30 minute, 5 minute, etc. intervals will also work.

# To simulate a day of photobiological dynamics, do as below.
person.simulate_day(I)
#       By default sleep is from midnight to 6 am and is provided as a 1 (awake) / 0 (asleep) schedule such as below,
#       sleep_schedule=[0]*6 + [1]*18
#       Different sleep schedules can be calculated using: person.simulate_day(I, S = sleep_schedule)
#       If S is set to False, sleep / wake will occur naturally based on the virtual person's H / C dynamics

# You can simulate as many days as you wish:
for i in range(6):  
    person.simulate_day(I) # now a week of a fixed sleep/wake and irradiance exposure has been calculated

# When completed with simulation, generate the alertness metrics:
person.generate_alertness_metrics()

# You can extract the 7 days of summary daily results
person.day_metrics[0].header_str() # a header string
person.day_metrics[0].output_str() # day 1's metrics
person.day_metrics[6].output_str() # day 7's metrics

# There are also some helper functions to parse ALFA-simulation export files into the CIE S 026 toolkit, including Ee,mel melanopic irradiance used by the biological_model
CCT_6100 = spectral_conversion.read_alfa_file('.\\example - ALFA electric light results\\neutral.csv')
CCT_6100[0].type # the first result is a workplane sensor
CCT_6100[13].type # the 14th result is a view sensor
print('%.2f,%.2f,%.2f\t\t%.2f,%.2f,%.2f' % (CCT_6100[13].X,CCT_6100[13].Y,CCT_6100[13].Z,CCT_6100[13].dX,CCT_6100[13].dY,CCT_6100[13].dZ)) # we can get its spatial location and direction information
ill, mel_ill, m_over_p, opic_irrads, opic_elrs, opic_ders, opic_edis = spectral_conversion.calc_s026(CCT_6100[13].get_data()) # calculate all of the CIE S 026 measures
opic_irrads["mel"] # melanopic irradiance
opic_edis["mel"] # melanopic equivalent D65 daylight illuminance

# With a little effort, the instantaneous model state data can also be extracted. A custom metrics class can facilitate this.
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

# then use the class to summarize metrics at a desired timestep, such as 30 minutes
# get median metrics every 30 minutes
min30_metrics = []
increment = int(30 * 60 / 20) # 30 minutes * 60 seconds/minute / 20 seconds/timestep
for n in range(0, int(7 * 24 * 60 * 60 / 20), increment):
    min30_metrics.append(median_metrics(person.history[n : n + increment],person.metrics[n : n + increment] ))

# write the data
with open('30-minute-data.csv', 'w') as f:
    f.write('day,clocktime')
    f.write(',%s' % min30_metrics[0].get_header_str())
    f.write('\n')
    for n_time in range(len(min30_metrics)):
        f.write('%i,%.3f' % (min30_metrics[n_time].day, min30_metrics[n_time].clocktime) )
        f.write(',%s' % min30_metrics[n_time].get_data_str() )
        f.write('\n')
# I = Ee,mel melanopic irradiance
# p_b = Melatonin in bloodstream, pmol/L
# p_bmax =  Melatonin in bloodstream without light exposure, pmol/L
# vPVTL = Reaction time, ms
# KSS = Sleepiness scale