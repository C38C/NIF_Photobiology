"""
17 Jan 2021 - Version 2 - Checking for inaccuracies in applied model... 
	- Verified that the correct v_YY equation is used! 
	- Fixed error in KSS regression function from Tekieh et al. 2020, multiplying theta_L_KSS by I in the updated equation (see _instantaneous).
	- Corrected Melatonin suppression calculation (but maybe need to implement on a day-to-day basis).
	- Implemented forced wake from 6am to 12am
	- Adjusted Forced Wake (F_w) conditions to be relative to sleep drive (D_v) > 2.4.

8 Apr 2021 - Version 3 - Updated based on email from Dr. Postnova and Clotilde's commentary
	- Corrected Sigmoid function and KSS calculation per Dr. Postnova and Clotilde
	- Theta_L_KSS normalized by S(1000.0 * 0.00080854)
	- Added linear interpolation of illuminance values rather than straight duplication. Sleep schedule is still duplication-based. 

17 May 2021 - Version 4
	- Corrected Theta_L_KSS calculation - normalized by S(1000) but with a new S_c = 1.0 / 223.5 as per Dr. Tekieh's and Dr. Postnova's emails.

19-20 June 2021 - Version 5
	- Allowed setting of day for annual calculations, useful if you want to simulate a subset of a year or start/stop at multiple locations. 
	- Calculate melatonin supression at morning / evening and entire day. 
	- Adjusted melatonin supression to be based on current entrainment level without effects of lighting, NOT the original entrained state. 
	- Phase shiftng can no longer randomly jump by 24 hours (checking for geometric solving error) 
	
"""

import math, time
from datetime import datetime
import numpy as np


class subject:

	def write_output(self, filepath, type="timestep"):
		
		if type == "timestep":
			if len(self.metrics) > 0:
				with open(filepath, 'w') as f:
					for i in range(len(self.history)):
						h = self.history[i]
						m = self.metrics[i]
						
						if i == 0:
							f.write(h.header_str() + ',' + m.header_str() + '\n')
						f.write(h.output_str() + ',' + m.output_str() + '\n')
			else:
				print("Please run subject.generate_alertness_metrics() before writing timestep output!")
				
				
				
		elif type == "daily":
			if len(self.day_metrics) > 0:
				with open(filepath, 'w') as f:
					for n,dh in enumerate(self.day_metrics):
						if n == 0:
							f.write(dh.header_str() + '\n')
						f.write(dh.output_str() + '\n')
			else:
				print("Please run subject.simulate_day(I) at least once before writing daily output!")
	
	def _interp(self, data):
		"""
		Interpolates data based on a variable number of steps, default = 72 -- 50 s
		"""
		
		steps = int( ((24.0 * 60.0 * 60.0) / len(data)) / self.timestep )
		
		interp_data = []
		for n in range( len(data) ): # 0-23 (len-1)
			this_hr = data[n]
			if n + 1 == len(data):
				next_hr = data[0]
			else:
				next_hr = data[n+1]
			for t in range(steps):
				weight_next_hr = float(t) / float(steps) #t is current time step
				weight_this_hr = 1.0 - weight_next_hr
				interp_data.append( this_hr * weight_this_hr + next_hr * weight_next_hr )
		return interp_data
	
	def t_crit(self, day_history):
		"""Based on Eq 16 / Postnova et al. / 2018		
		"""
		# This is a little weird, because -2.98 is arbitrary AND it will twice in a normally-entrained 24 hour period. Reasonably, I think, I search for when -2.98 occurs only when Y > X."
		
		phi_crit = -2.98 # when = atan(Y/X)
		dist_phi = 10000.0
		# first find which timestep phi_crit exists in using a brute force search
		
		n_crit = 100
		
		for n,s in enumerate(day_history):
			this_phi = math.atan2(s.Y, s.X)
			if abs(phi_crit - this_phi) < dist_phi:
				dist_phi = abs(phi_crit - this_phi)  
				n_crit = n
				
		t_phi_crit = float(n_crit) * self.timestep # time in seconds, start of day = 0
		
		
		
		t_0_melpeak = 0.7 * 3600.0
		t_0_cbtmin = 2.7 * 3600.0
		
		t_crit_melpeak = t_phi_crit + t_0_melpeak
		t_crit_cbtmin = t_phi_crit + t_0_cbtmin
		
		if t_crit_melpeak > 86400:
			t_crit_melpeak = t_crit_melpeak - 86400
		if t_crit_cbtmin > 86400:
			t_crit_cbtmin = t_crit_cbtmin - 86400
		
		# n_crit_melpeak = (t_crit_melpeak / self.timestep)
		
		return (t_crit_melpeak, t_crit_cbtmin)
	
	# def _set_base_melatonin(self, t_crit_melpeak):
	#	 """Function to align melatonin peak with Postnova et al. 2018."""
	#	 
	#	 n_crit_melpeak = int( t_crit_melpeak / self.timestep )
	#	 n_last_melpeak = int( self._t_crit_melpeak / self.timestep )
	#	 
	#	 # find current melatonin peak location
	#	 maxval = max(self._base_melatonin)
	#	 for n,m in enumerate(self._base_melatonin):
	#		 if m == maxval:
	#			 n_last_melpeak = n
	#			 break
	#	 
	#	 # offset melatonin in a new temporary array 
	#	 new_base_melatonin = []
	#	 max_n = len(self._base_melatonin) - 1
	#	 delta_n = n_crit_melpeak - n_last_melpeak
	#	 for n in range(len(self._base_melatonin)):  
	#		 if delta_n >= 0: # phase delay
	#			 new_base_melatonin.append( self._base_melatonin[n - delta_n] )
	#		 else: # phase advance
	#			 if n - delta_n > max_n:
	#				 new_base_melatonin.append( self._base_melatonin[n - delta_n - max_n - 1] )
	#			 else:
	#				 new_base_melatonin.append( self._base_melatonin[n - delta_n] )
	#			 
	#	 
	#	 self._base_melatonin = new_base_melatonin
	
	def generate_alertness_metrics(self):
		self.metrics = [] # empty current metrics as the delta* measures vary with total mean H and C
		
		# first calculate mean values for H and C
		mean_H = 0.0
		mean_C = 0.0   
		for s in self.history:
			mean_H += s.H
			mean_C += s.C
		mean_H = mean_H / len(self.history)
		mean_C = mean_C / len(self.history)
		
		thetaH_deltaADD = -1.35
		thetaC_deltaADD = 1.42
		
		thetaH_deltaDSST = -6.29
		thetaC_deltaDSST = 6.97
		
		for s in self.history:
			if s.S == 1: # gate metrics by sleep/wake state
				vPVTL, KSS, vPVTRT_mean, vPVTRT_median, vPVTRT_fast10, vPVTRT_slow10 = self._instantaneous(s.H, s.C, s.theta_L_KSS, s.I)
				deltaADD = thetaH_deltaADD * (s.H - mean_H) + thetaC_deltaADD * (s.C - mean_C)
				deltaDSST = thetaH_deltaDSST * (s.H - mean_H) + thetaC_deltaDSST * (s.C - mean_C)
			else:
				vPVTL, KSS, vPVTRT_mean, vPVTRT_median, vPVTRT_fast10, vPVTRT_slow10 = (0,0,0,0,0,0)
				deltaADD = 0
				deltaDSST = 0
			
			alertness = self.alertness_metrics(s.day, s.n, vPVTL, KSS, vPVTRT_mean, vPVTRT_median, vPVTRT_fast10, vPVTRT_slow10, deltaADD, deltaDSST)
			
			self.metrics.append(alertness)
		
	def _instantaneous(self, H, C, theta_L_KSS, I):
		"""Instantaneous alertness measures vPVTL, KSS, vPVTRT - Based on Postnova et al., 2018 AND Tekieh et al. 2020...
			vPVTL - visual Performance Vigilance Test - Total Lapses
			KSS - Karolinska Sleepiness Scale
			vPVTRT - visual Performance Vigilance Test - Reaction Time
				- mean, median, fastest 10%, slowest 10%"""
		
		# cPVTL - total lapses
		c_vPVTL = -137.49
		thetaH_vPVTL = 11.77
		thetaC_vPVTL = -12.47
		vPVTL = c_vPVTL + (thetaH_vPVTL ) * H + thetaC_vPVTL * C
		
		
		# KSS - sleepiness
		c_KSS = -24.34
		thetaH_KSS = 2.28
		thetaC_KSS = -1.74
		KSS = c_KSS + (thetaH_KSS + theta_L_KSS) * H + thetaC_KSS * C
		
        # KSS = -24.34 + (2.28 + theta_L_KSS) * H + -1.74 * C
		
		# vPVTRT - reaction time
		c_vPVTRT_mean = -12787.0
		thetaH_vPVTRT_mean = 1055.0
		thetaC_vPVTRT_mean = -1144.0
				
		c_vPVTRT_median = -1736.0
		thetaH_vPVTRT_median = 168.0
		thetaC_vPVTRT_median = -214.0
		
		c_vPVTRT_fast10 = -4.99
		thetaH_vPVTRT_fast10 = 17.77
		thetaC_vPVTRT_fast10 = -18.75
		
		c_vPVTRT_slow10 = -12473.0
		thetaH_vPVTRT_slow10 = 1096.0
		thetaC_vPVTRT_slow10 = -1830.0
		
		vPVTRT_mean = c_vPVTRT_mean + (thetaH_vPVTRT_mean ) * H + thetaC_vPVTRT_mean * C
		vPVTRT_median = c_vPVTRT_median + (thetaH_vPVTRT_median ) * H + thetaC_vPVTRT_median * C
		vPVTRT_fast10 = c_vPVTRT_fast10 + (thetaH_vPVTRT_fast10 ) * H + thetaC_vPVTRT_fast10 * C
		vPVTRT_slow10 = c_vPVTRT_slow10 + (thetaH_vPVTRT_slow10 ) * H + thetaC_vPVTRT_slow10 * C
		
		
		return (vPVTL, KSS, vPVTRT_mean, vPVTRT_median, vPVTRT_fast10, vPVTRT_slow10)
	
	def _internal_simulate_day(self, starting_state, I_list, S_list=False ):
		day_history = [starting_state]
		
		timestep = self.timestep
		
		total_steps = int((24 * 60 * 60) / timestep)

		for i_t in range( total_steps ): # i_t being the time index
			this_state = self.circadian_state(state = day_history[-1])
			this_state.I = I_list[i_t]
			
			this_F_w = 0.0
			
			# photic drive, P
			if S_list:
				this_S = S_list[i_t]
				this_state.S = this_S
				
				if this_S > 0.9: # and (day_history[-1].D_v > 2.4 or this_state.get_S(day_history[-1].V_m) < 0.1)
					this_F_w = 1.0
				
			else:
				this_S = day_history[-1].S			
				this_state.S = this_state.get_S(day_history[-1].V_m)
			
			this_state.F_w = this_F_w
			
			this_state.alpha = this_state.get_alpha(I_list[i_t] * this_S)
			this_state.P = day_history[-1].P + this_state.get_dP(day_history[-1].alpha, day_history[-1].P)
			this_state.D_p = this_state.get_D_p(day_history[-1].alpha, day_history[-1].P, day_history[-1].X, day_history[-1].Y)
			
			# nonphotic drive, D_n
			this_state.D_n = this_state.get_D_n(this_S, day_history[-1].X)
			
			# circadian drive, C
			this_state.X = day_history[-1].X + this_state.get_dX(day_history[-1].Y, day_history[-1].X, day_history[-1].D_p, day_history[-1].D_n)
			this_state.Y = day_history[-1].Y + this_state.get_dY(day_history[-1].Y, day_history[-1].X, day_history[-1].D_p)
			this_state.C = this_state.get_C(day_history[-1].X, day_history[-1].Y)
			
			# MA: V_m, Q_m
			this_state.V_m = day_history[-1].V_m + this_state.get_dV_m(day_history[-1].V_m, day_history[-1].Q_v, day_history[-1].W)
			this_state.Q_m = this_state.get_Q_i(day_history[-1].V_m)
			
			# Wake effort, W
			this_state.W = this_state.get_W(day_history[-1].Q_v, F_w = this_F_w)
			
			# homeostatic drive, H
			this_state.H = day_history[-1].H + this_state.get_dH(day_history[-1].Q_m, day_history[-1].H)
			
			# VLPO: V_v, Q_v 
			this_state.D_v = this_state.get_D_v(day_history[-1].H, day_history[-1].C)
			this_state.V_v = day_history[-1].V_v + this_state.get_dV_v(day_history[-1].V_v, day_history[-1].Q_m, day_history[-1].D_v)
			this_state.Q_v = this_state.get_Q_i(day_history[-1].V_v)
			
			# Melatonin dynamics
			this_state.r = this_state.get_r(I_list[i_t])
			this_state.At = day_history[-1].At + this_state.get_dAt(day_history[-1].At, day_history[-1].Y, day_history[-1].X, day_history[-1].r)
			this_state.At_max = day_history[-1].At_max + this_state.get_dAt(day_history[-1].At_max, day_history[-1].Y, day_history[-1].X, 1.0)
			this_state.p_b = day_history[-1].p_b + this_state.get_dp_b(day_history[-1].At, day_history[-1].p_b)
			this_state.p_b_max = day_history[-1].p_b_max + this_state.get_dp_b(day_history[-1].At_max, day_history[-1].p_b_max)
			
			# alerness dynamics due to light
			this_state.theta_L_KSS = day_history[-1].theta_L_KSS + this_state.get_dtheta_L(I_list[i_t], day_history[-1].theta_L_KSS)
			
			if i_t < total_steps - 1: 
				day_history.append(this_state)
			
			
		
		# max_theta = -99999
		# min_theta = 99999
		# sum_theta = 0.0
		# c = 0.0
		# for s in day_history:
		#	 if not s.S:
		#		 c += 1.0
		#		 theta = math.atan2(s.Y, s.X)
		#		 sum_theta += theta
		#		 if theta > max_theta:
		#			 max_theta = theta
		#		 if theta < min_theta:
		#			 min_theta = theta
		# 
		# print ("max", max_theta)
		# print ("min", min_theta)
		# print ("mean", sum_theta / c)
		# print "\n"
		
		return (this_state, day_history)
	
	def simulate_day(self, I, S=[0]*6 + [1]*18): # 
		if len(I) < int((24 * 60 * 60) / self.state.dt):
			# tmp = []
			
			tmp = []
			
			x_val = np.linspace(0,23, len(I))
			x_interp = np.linspace(0,23, int((24 * 60 * 60) / self.state.dt))
			
			I = np.interp(x_interp, x_val, I)
			
			# n_duplicates = int( ( 24.0 * 60.0 * 60.0 /  len(I)) / self.state.dt)
			# for this_I in I:
			# 	tmp.extend([this_I] * n_duplicates)
			# 
			# I = tmp
		if S:
			if len(S) < int((24 * 60 * 60) / self.state.dt):
				tmp = []
				
				n_duplicates = int( ( 24.0 * 60.0 * 60.0 /  len(S)) / self.state.dt)
				for this_S in S:
					tmp.extend([this_S] * n_duplicates)
				
				S = tmp
		
		
		# print len(I),len(S)
		
		self.state, day_history = self._internal_simulate_day(self.state, I,  S_list=S)
		
		t_crit_melpeak, t_crit_cbtmin = self.t_crit(day_history)

		self.history.extend(day_history)
		
		# generate day metrics
		phase_shift = t_crit_melpeak - self._t_crit_melpeak
		self._t_crit_melpeak = t_crit_melpeak # update class t_crit
		total_actual_melatonin = 0.0
		total_actual_melatonin_morning = 0.0
		total_actual_melatonin_evening = 0.0
		total_potential_melatonin = 0.0
		total_potential_melatonin_morning = 0.0
		total_potential_melatonin_evening = 0.0
		
		for n_s, s in enumerate(day_history):
			total_potential_melatonin += s.p_b_max
			if n_s < len(day_history) / 2:
				total_potential_melatonin_morning += s.p_b_max
			else:
				total_potential_melatonin_evening += s.p_b_max
		for n_s, s in enumerate(day_history):
			total_actual_melatonin += s.p_b
			if n_s < len(day_history) / 2:
				total_actual_melatonin_morning += s.p_b
			else:
				total_actual_melatonin_evening += s.p_b
		
		melatonin_supression_percent = 100.0 * ((total_potential_melatonin - total_actual_melatonin) / total_potential_melatonin)
		melatonin_supression_percent_morning = 100.0 * ((total_potential_melatonin_morning - total_actual_melatonin_morning) / total_potential_melatonin_morning)
		melatonin_supression_percent_evening = 100.0 * ((total_potential_melatonin_evening - total_actual_melatonin_evening) / total_potential_melatonin_evening)
		self.day_metrics.append( self.daily_metrics( day_history[-1].day, t_crit_melpeak, t_crit_cbtmin, phase_shift, melatonin_supression_percent, melatonin_supression_percent_morning, melatonin_supression_percent_evening) )

		self.state.day += 1 # incremement day
			
	def entrain(self, entrain_state, ndays=25, logfile=False):
		# start_time = time.time()
		# print("Beginning circadian entrainment procedure at time = %s." % datetime.now().strftime("%H:%M:%S"))
		
		entrain_history = []
		
		# forced entraining variables, 30 minute schedules ###
		# S, sleep/wake state
		S_entrain_30min = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # 00:00 - 08:00
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, # 08:00 - 20:00
			1, 1, 1, 1, 1, 1, 1, 0] # 20:00-24:00
		S_list = []
		
		# I, illuminance on eye @ 4100K CCT
		I_entrain_30min = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # 00:00 - 08:00
			250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, # 08:00 - 20:00
			40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 0] # 20:00-24:00
		I_entrain_30min = [I * 0.00080854 for I in I_entrain_30min] # convert to melanopic irradiance
		
		I_list = []
		
		# create 10-minute lists, no interpolation
		n_duplicates = int( (30 * 60) / entrain_state.dt)
		for i in range( len( S_entrain_30min) ):
			S_list.extend( [S_entrain_30min[i]] * n_duplicates )
			I_list.extend( [I_entrain_30min[i]] * n_duplicates )
		
		if logfile: # if logging is enabled at entrainment phase
			f = open(logfile, 'w')
			f.write(entrain_state.header_str() + '\n')
		
		# calculate ndays with a fixed sleep / wake schedule
		for i in range(ndays):
			# enforced sleep / wake schedule
			# forced wake at 8:00 each day
			entrain_state, entrain_day_history = self._internal_simulate_day(entrain_state, I_list, S_list = S_list)
			
			t_crit_melpeak, t_crit_cbtmin = self.t_crit(entrain_day_history)
			# self._set_base_melatonin(t_crit_melpeak)
			
			entrain_history.extend(entrain_day_history)
		
		if logfile: # if logging is enabled at entrainment phase
			for s in entrain_history:
				f.write(s.output_str() + '\n')
			f.close()
		
		# calculate AUC_baseline as per Tekieh et al. 2020 [33,36]
		AUC_melatonin_baseline = 0.0
		AUC_melatonin_baseline_morning = 0.0
		AUC_melatonin_baseline_evening = 0.0
		
		for n_s,state in enumerate(entrain_day_history):
			AUC_melatonin_baseline += state.p_b
			if n_s < len(entrain_day_history) / 2:
				AUC_melatonin_baseline_morning += state.p_b
			else:
				AUC_melatonin_baseline_evening += state.p_b
		
		return entrain_state, t_crit_melpeak, AUC_melatonin_baseline, AUC_melatonin_baseline_morning, AUC_melatonin_baseline_evening
		  
	
	class circadian_state: 
		def __init__(self, timestep=20.0, state=False, day = 0):
			"""Initializes as estimated parameters for a default, well-entrained individual as in Postnova et al. 2014, 2016, 2018. The parent class entrain() function will attempt to create a proper starting state for circadian modelling, but one could also set these model starting parameters directly."""
			# time constants, so many ###
			self.t_v = self.t_m = 50.0 # s
			self.t_H = 59.0 * 3600.0  # s
			self.t_X = self.t_Y = 24.0 * 3600.0 / (2.0 * math.pi) # s
			self.t_C = 24.2 * 3600.0 # s
			self.t_A = 1.5 * 60.0 * 60.0 # s, melatonin synthesis constant as per Tekieh et al. 2020 (NOT as per Abeysuriya et al. 2018)
			
			# defineable simulation time step - times > 30s tend to result in noise in the VLPO + MA systems
			# 10 and 20 seconds both seem to produce smooth results
			self.dt = timestep # note, if state this will be overridden by previous circadian state dt / timestep value
			
			if not state: # if no input state..
				self.n = 0
				
				self.I = 0
				
				# melatonin dynamics
				self.p_b = 200.0 # Abeysuriya et al. 2018, estimated from Fig. 2, pmol / L
				self.p_b_max = 200.0 # Abeysuriya et al. 2018, estimated from Fig. 2, pmol / L
				self.At = 0.47 / 2.0 # Abeysuriya et al. 2018, Table 1 - guess based on value for A_0
				self.At_max = 0.47 / 2.0 # Abeysuriya et al. 2018, Table 1 - guess based on value for A_0
				self.r = 1.0 # Abeysuriya et al. 2018 / suppression of melatonin synthesis 
				
				# alertness dynamics due to light 
				self.theta_L_KSS = 0.0 # starts at 0.0 as there is no irradiance present

				self.day = day
				# starting dynamic variables / differential equation-governed ###
				self.V_v = 1.5 # Philips and Robinson 2007, estimated from Fig. 5 # sleep active voltage
				self.V_m = -15.0 # Postnova et al. 2014, estimated from Fig. 2 # wake active voltage
				self.H = 13.0 # Postnova et al. 2014, estimated from Fig. 2 
				self.X = 0.04  # St Hilaire et al. 2007, estimated from p 7 / 17
				self.Y = -1.28 # estimated by solving C_adj in Postnova 2016 for X = 0, C = 1
				
				self.P = 0.0 # guess based on no photoreceptors firing during sleep
				
				# starting dynamic variables / straightforward ###
				self.Q_v = self.get_Q_i(self.V_v) # 5.55 per Postnova 2018 equation
				# 2.0 # Phillips and Robinson 2007 Table 1, NREM value # sleep active firing rate
				self.Q_m = self.get_Q_i(self.V_m) # 0.02 per Postnova 2018 equation
				# 1.1 # Phillips and Robinson 2007 Table 1 DR + LC - 1/2 SEM data, NREM values # wake 
				self.W = self.get_W(self.Q_v) # will be 0 when not during a forced wake
				# if self.W < 0: # check for wake effort < 0
				self.W = 0
				self.S = 0.0 # model starts asleep at clock time 00:00
				self.C = 1.0 # Postnova et al. 2016
				self.D_v = self.get_D_v(self.H, self.C) # 4.2 per Postnova 2018 equation, ~3 per Postnova 2014 Fig 2
				self.alpha = self.get_alpha(0.0) # light signal strength scaled from 100 - 9500 lx
				self.D_n = self.get_D_n(self.S, self.X)
				self.D_p = self.get_D_p(self.alpha, self.P, self.X, self.Y)
				
			else: # initialize from input state of class circadian_state
				self.dt = state.dt
				
				self.I = 0
				
				self.p_b = state.p_b
				self.p_b_max = state.p_b_max
				
				self.At = state.At
				self.At_max = state.At_max
				
				self.r = state.r
				
				self.theta_L_KSS = state.theta_L_KSS
				
				self.day = state.day
				self.n = state.n + 1
				self.V_v = state.V_v
				self.V_m = state.V_m
				self.H = state.H
				self.X = state.X
				self.Y = state.Y
				self.P = state.P
				self.Q_v = state.Q_v
				self.Q_m = state.Q_m
				self.W = state.W
				self.S = state.S
				self.C = state.C
				self.D_v = state.D_v
				self.alpha = state.alpha
				self.D_n = state.D_n
				self.D_p = state.D_p
		
		def header_str(self):
			return "day,clocktime,n,V_v,V_m,H,X,Y,P,Q_v,Q_m,W,S,C,D_v,alpha,D_n,D_p,I,r,At,At_max,p_b,p_b_max,theta_L_KSS"
		
		def output_str(self):			
			str_list = []
			str_list.append('%i' % self.day)
			clocktime = self.n / (60.0 * 60.0 / self.dt) - ((self.day % 7.0) * 24.0)
			str_list.append('%.3f' % clocktime)
			str_list.append('%i' % self.n)
			str_list.append('%.5f' % self.V_v)
			str_list.append('%.5f' % self.V_m)
			str_list.append('%.5f' % self.H)
			str_list.append('%.5f' % self.X)
			str_list.append('%.5f' % self.Y)
			str_list.append('%.5f' % self.P)
			str_list.append('%.5f' % self.Q_v)
			str_list.append('%.5f' % self.Q_m)
			str_list.append('%.5f' % self.W)
			str_list.append('%s' % self.S)
			str_list.append('%.5f' % self.C)
			str_list.append('%.5f' % self.D_v)
			str_list.append('%.5f' % self.alpha)
			str_list.append('%.5f' % self.D_n)
			str_list.append('%.5f' % self.D_p)
			str_list.append('%.1f' % self.I)
			str_list.append('%.5f' % self.r)
			str_list.append('%.5f' % self.At)
			str_list.append('%.5f' % self.At_max)
			str_list.append('%.5f' % self.p_b)
			str_list.append('%.5f' % self.p_b_max)
			str_list.append('%.5f' % self.theta_L_KSS)
			return ','.join(str_list)
		
		# circadian model functions ###
		# based on Postnova et al. 2018, 2016, 2014; Phillips and Robinson 2007, St. Hillaire et al. 2007
		def get_dV_v(self, V_v, Q_m, D_v):
			"""Eq 1 / Postnova et al. / 2018"""
			v_vm = -2.1 # mV		
			# V_v, mean voltage VLPO
			# Q_m, mean firing rate MA
			# D_v, total sleep drive
			
			return (v_vm * Q_m - V_v + D_v) * (self.dt / self.t_v)
		
		def get_dV_m(self, V_m, Q_v, W):
			"""Eq 2 / Postnova et al. / 2018"""
			v_mv = -1.8 # mV
			D_m = 1.3 # mV
			# V_m, mean voltage MA
			# Q_v, mean firing rate VLPO
			# W, wake effort
			
			return (v_mv * Q_v - V_m + D_m + W) * (self.dt / self.t_m)
		
		def get_dH(self, Q_m, H): 
			"""Eq 3 / Postnova et al. / 2018"""
			v_Hm = 4.57 # s
			# Q_m, mean firing rate MA
			# H, homeostatic drive
			
			return (v_Hm * Q_m - H) * (self.dt / self.t_H)
		
		def get_dX(self, Y, X, D_p, D_n):
			""""Eq 4 / Postnova et al. / 2018"""
			y = 0.13
			v_Xp = 37.0 * 60.0
			v_Xn = 0.032
			# Y, circadian variable Y
			# X, circadian variable X
			# D_p, photic drive
			# D_n, nonphotic drive
		   
			return (Y + y * (X/3.0 + (4.0/3.0) * X**3.0 - (256.0/105.0) * X**7.0) + v_Xp * D_p + v_Xn * D_n) * (self.dt / self.t_X)
		
		def get_dY(self, Y, X, D_p):
			"""Eq 5 / Postnova et al. / 2018"""
			v_YY =  (37.0 * 60.0) / 3.0 # 1.0 / (3.0 * 37.0 * 60.0) # ( # Could be 1.0 / (3.0 * 37.0 * 60.0) - paper is unclear
			v_YX = 0.55 * 37.0 * 60.0
			d = 24.0 * 3600.0 / 0.99729
			# Y, circadian variable Y
			# X, circadian variable X
			# D_p, photic drive
			
			return (D_p * (v_YY * Y - v_YX * X) - (((d / self.t_C)**2.0) * X)) * (self.dt / self.t_Y)
		
		def get_dP(self, alpha, P):
			"""Eq 6 / Postnova et al. / 2018"""
			
			beta = 0.007 / 60.0 # 1/s
			# alpha, function of photic drive
			# P, photorecptor activity
			return (alpha * (1.0 - P) - beta * P) * self.dt
		
		def get_Q_i(self, V_i):
			"""Eq 7 / Postnova et al. / 2018
			Mean population firing rate
			Used for v (VLPO) and m (MA) subscript calculations."""
			Q_max = 100.0 # 1/search_str
			theta = 10.0 # mV
			sigma_prime = 3.0 # mV
			# V_i, mean voltage of _i system
			
			return (Q_max / (1.0 + math.exp( (theta - V_i) / sigma_prime) ))
		
		def get_W(self, Q_v, F_w = 0):
			"""Eq 8 / Postnova et al. / 2018
			Wake effort"""
			V_WE = -0.07 # mV
			v_mv = -1.8 # mV
			D_m = 1.3 # mV
			# Q_v, mean firing rate VLPO
			
			return (F_w * max([0.0, V_WE - v_mv * Q_v - D_m]))
		
		def get_D_v(self, H, C):
			"""Eq 9 / Postnova et al. / 2018
			Drive to sleep active neurons (VLPO)"""
			v_vH = 1.0
			v_vC = -0.5 # mV
			A_v = -10.3 # mV
			# H, homeostatic drive
			# C, circadian drive
			
			return (v_vH * H + v_vC * C + A_v)
		
		def get_C(self, X, Y):
			"""Eq 10 / Postnova et al. / 2018
			Circadian drive"""
			# Y, circadian variable Y
			# X, circadian variable X	
			
			return 0.1 * ((1.0 + X)/2.0) + ((3.1 * X - 2.5 * Y + 4.2) / (3.7 * (X + 2.0)))**2.0
		
		def get_D_n(self, S, X):
			"""Eq 11 / Postnova et al. / 2018
			Nonphotic drive to the circadian"""
			r = 10.0
			# S, sleep / Wake state
			# X, circadian variable X	
			
			return (S - 2.0/3.0) * (1.0 - (math.tanh(r * X)))
		
		def get_D_p(self, alpha, P, X, Y):
			"""Eq 12 / Postnova et al. / 2018
			Photic drive to the circadian"""
			epsilon = 0.4
			# alpha, function of photic drive
			# P, photorecptor activity
			# Y, circadian variable Y
			# X, circadian variable X
			
			return alpha * (1.0 - P) * (1 - epsilon * X) * (1 - epsilon * Y)
		
		def get_alpha(self, I):
			"""Eq 13 / Postnova et al. / 2018"""
			alpha_0 = 0.1 / 60.0 # 1/s
			I_0 = 100 # lx
			I_1 = 9500 # lx
			
			F_4100K = 0.00080854 # F_4100K is from Tekieh et al.'s 2020 modification to the Postnova, et al. 2018 paper accounting for melanopic irradiance. Any other spectrally-weighted conversion could be used for this (i.e. Melanopix Lux, Melanopic EDI, etc.).
			
			return alpha_0 * (I / (I + I_1 * F_4100K)) * (I / (I_0 * F_4100K))**0.5
		
		# Eq 14 is just illuminance times the 1 or 0 S (sleep / wake state) variable.
		
		def get_S(self, V_m):
			"""Eq 15 / Postnova et al. / 2018
			Sleep / wake state"""
			V_th = -2.0 # mV
			# V_m, mean voltage MA

			if V_m > V_th:
				return 1.0
			else:
				return 0.0 
		
		def get_r(self, I):
			"""Eq 9 / Tekieh et al. / 2020
			Melatonin supprssion function based on melanopic _irradiance_"""
			
			# #######################################################################
			# Below is the Zeitzer equation based on instantaneous illuminance, I...
			# #######################################################################
			# """Four parameter logistic regression from Zeitzer et al, 2000 - acute supression of plasma melatonin. Returns 0 for 100% supression and 1 for no supression -- the inversion of the Zeitzer function."""
			# a = - 0.0156
			# b = 88.0
			# c = 1.0
			# d = 3.55
			# 
			# supression_ratio = abs( (a - c) / (1.0 + (I / b)**d) + c)
			# if supression_ratio < 0.0 : supression_ratio = 0.0
			# if supression_ratio > 1.0 : supression_ratio = 1.0
			# 
			# return 1.0 - supression_ratio
			
			# Tekieh et al. / 2020
			r_a = 1.0
			r_b = 0.031 # E_e,mel
			r_c = 0.82				
			
			if I <= 0.00000001:
				return 1.0
			else:
				return 1.0 - (r_a / ( 1.0 + (( I / r_b ) ** (-1.0 * r_c)) ))
		
		def get_dAt(self, At, Y, X, r):
			"""Eq 8 / Tekieh et al. / 2020
			Values from Abeysuriya 2018
			Melatonin synthesis rate"""
			theta_on = -1.44 # radians
			theta_off = 2.78 # radians
			
			theta = math.atan2(Y, X)
			
			if theta < theta_on or theta > theta_off:
				m = 1.0
			else:
				m = 0.0
			
			return ( m * r - At) * (self.dt / self.t_A)
		
		def get_dp_b(self, At, p_b):
			"""Eq 14 / Abeysuriya et al. / 2018
			Melatonin concentation in blood plasma."""
			u_star = 0.47 # pmol/(L*s), urine peak melatonin 
			p_star_b = 325.0 # pmol/L, plasma peak melatonin
			r_g = 0.9 # fraction of melatonin converted to aMT6s in the urine 
			
			return (((u_star * At) / r_g) - u_star * (p_b / (r_g * p_star_b) )) * self.dt
		
		# def get_u(self, p_b__t_u):
		#	 """Eq 15 / Abeysuriya et al. / 2018
		# 	Urine aMT6s excretion rate.
		# 	Currently not used by the code."""
		#	 p_star_b = 325.0 # pmol/L, plasma peak melatonin
		#	 t_u = 0.96 * 3600.0 # s, urine time lag
		#	 u_star = 0.47 # pmol/(L*s), urine peak aMT6s 
		#	 # p_b__t_u is the blood plasma at some time lag in the past
		#	 
		#	 return u_star * ( (p_b__t_u) / p_star_b)
		
		def get_dtheta_L(self, I, theta_L):
			"""Eq 12-14 / Tekieh et al. / 2020
			Instantaneous alerting effects of light."""
			
			# I_Eemel = I / (10.0**-2.0) # conversion from Tekieh.. confusion
			
			# Currently only for KSS estimation!
			t_L = 24.0 * 60.0 # seconds
			v_LA = -0.11 # ?? mA
 
			def S_Ee(I):
				S_b = 0.05 # W/m2, melanopic irrad
				S_c = 1.0 / 223.5 # m2/W, melanopic irrad 
				return 1.0 / (1.0 + math.exp((S_b - I) / S_c) )
			
			S_norm = ( S_Ee(I) - S_Ee(0.0) ) / ( S_Ee(1000.0) - S_Ee(0.0) )
			
			return (-1.0 * theta_L + v_LA * S_norm) * (self.dt / t_L)
			
	
	class daily_metrics: 
		"""All time measures in seconds from midnight or delta-seconds."""
		def header_str(self, day = True):
			str = "t_crit_melpeak_s,t_crit_melpeak_h,t_crit_cbtmin_s,t_crit_cbtmin_h,phase_shift_s,melatonin_supression,melatonin_supression_morning,melatonin_supression_evening"
			
			if day:
				str = "day," + str
			
			return(str)

		
		def output_str(self, day = True):
			str_list = []
			if day: 
				str_list.append('%i' % self.day)
			str_list.append('%i' % self.t_crit_melpeak_s)
			str_list.append('%.2f' % (self.t_crit_melpeak_s / 3600.0))
			str_list.append('%i' % self.t_crit_cbtmin_s)
			str_list.append('%.2f' % (self.t_crit_cbtmin_s / 3600.0))
			str_list.append('%i' % self.phase_shift_s)
			str_list.append('%.1f' % self.p_b_supression_percent)
			str_list.append('%.1f' % self.p_b_supression_percent_morning)
			str_list.append('%.1f' % self.p_b_supression_percent_evening)
			
			return ','.join(str_list)
		
		def __init__(self, day, t_crit_melpeak_s, t_crit_cbtmin_s, phase_shift_s, melatonin_supression_percent, melatonin_supression_percent_morning, melatonin_supression_percent_evening):
			self.day = day
			self.t_crit_melpeak_s = t_crit_melpeak_s
			self.t_crit_cbtmin_s = t_crit_cbtmin_s
			self.phase_shift_s = phase_shift_s
			self.p_b_supression_percent = melatonin_supression_percent
			self.p_b_supression_percent_morning = melatonin_supression_percent_morning
			self.p_b_supression_percent_evening = melatonin_supression_percent_evening
	
	class alertness_metrics:
		def header_str(self):
			return "vPVTL,KSS,vPVTRT_mean,vPVTRT_median,vPVTRT_fast10,vPVTRT_slow10,deltaADD,deltaDSST"
		
		def output_str(self):
			str_list = []
			# str_list.append('%i' % self.day)
			# str_list.append('%i' % self.n)
			str_list.append('%.5f' % self.vPVTL)
			str_list.append('%.5f' % self.KSS)
			str_list.append('%.5f' % self.vPVTRT_mean)
			str_list.append('%.5f' % self.vPVTRT_median)
			str_list.append('%.5f' % self.vPVTRT_fast10)
			str_list.append('%.5f' % self.vPVTRT_slow10)
			str_list.append('%.5f' % self.deltaADD)
			str_list.append('%.5f' % self.deltaDSST)
			return ','.join(str_list)
		
		def __init__(self, day, n, vPVTL, KSS, vPVTRT_mean, vPVTRT_median, vPVTRT_fast10, vPVTRT_slow10, deltaADD, deltaDSST):
			self.day = day
			self.n = n
			self.vPVTL = vPVTL
			self.KSS = KSS
			self.vPVTRT_mean = vPVTRT_mean
			self.vPVTRT_median = vPVTRT_median
			self.vPVTRT_fast10 = vPVTRT_fast10
			self.vPVTRT_slow10 = vPVTRT_slow10
			self.deltaADD = deltaADD
			self.deltaDSST = deltaDSST
	
	def __init__(self, timestep = 20.0, do_entrain=True, log_entrain=False, day = 0):
		"""timestep in seconds -- 20.0 is a good default.
		log_entrain should be a filepath if logging is desired."""
		
		base_melatonin_24 = [20, 30, 48, 56, 55, 60, 55, 40, 35, 20, 11, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 13]
		
		self._t_crit_melpeak = 5.0 * 60.0 * 60.0 # starting melpeak time, 5 am  
		
		self.history = [] # contains history of circadian states, always updated with simulate_day
		
		self.metrics = [] # contains list of subjective metrics equal in lenth to self.history, but is kept empty until subject.generate_alertness_metrics() is run. 
		self.day_metrics = [] # contains a list of physiological measures (phase shifting, percentage of daily melatonin supressed) that can only be assessed over a circadian period. 
		
		self.timestep = timestep
		
		self.state = self.circadian_state(timestep = self.timestep, day = 0)
		
		self._base_melatonin = self._interp(base_melatonin_24)
		self._last_melatonin = self._interp(base_melatonin_24)
		
		# print len(self._base_melatonin)
		
		if do_entrain:
			self.state, self._t_crit_melpeak, self._AUC_melatonin_baseline, self._AUC_melatonin_baseline_morning, self._AUC_melatonin_baseline_evening = self.entrain(self.state, logfile = log_entrain)
			self.state.n = 0
			self.state.day = day


# class spectral_metrics:
#	 def __init__(self, line):
#		 parts = line.rstrip().split(',')
#		 
#		 try:
#			 self.XYZ = (float(parts[6]), float(parts[7]), float(parts[8]))
#		 except:
#			 self.XYZ = (float(int(parts[6])), float(int(parts[7])), float(int(parts[8])))
#		 try:
#			 self.XYZ_dir = (float(parts[9]), float(parts[10]), float(parts[11]))
#		 except:
#			 self.XYZ = (float(int(parts[9])), float(int(parts[10])), float(int(parts[11])))
#		 
#		 self.e = float(parts[12])
#		 self.e_m = float(parts[13])
#		 self.m_over_p = float(parts[14])
#		 self.m_irrad = float(parts[19])
#		 self.m_elr = float(parts[24])
#		 self.m_der = float(parts[29])
#		 self.m_edi = float(parts[34])
# 
# def get_time(mo, da, hr, min, sky_type ="clear", viewstr="View 1", path = "c:\\Users\\ajakubiec\\Dropbox\\NSERC_USRG_2020\\Athina\\Scripts\\NgTengFong_view_output_EL.csv"):
#	 """A function to grab a certain sky condition from a cumulative ALFA->CIE C026 results file."""
#	 
#	 global spectral_metrics
#	 
#	 with open(path, 'r') as f:
#		 line = f.readline()
#		 search_str = "%s,%i,%i,%i,%i,%s" % (sky_type, mo, da, hr, min, viewstr)
#		 while not search_str in line:
#			 line = f.readline()
#	 
#	 this_metrics = spectral_metrics(line)
#	 
#	 return this_metrics
#	 
# def get_electric(type, dim_fraction = 1.0, ath = "c:\\Users\\ajakubiec\\Dropbox\\NSERC_USRG_2020\\Athina\\Scripts\\NgTengFong_view_output.csv"):
#	 """A function to grab a certain electric lighting condition from a ALFA->CIE C026 results file."""
#	 
#	 global spectral_metrics
#	 pass
#	 