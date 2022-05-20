"""
Code to produce seasonal summary files from an annual calculation.
"""

import statistics as stats

views = ["L1", "L2", "L3"]


dir_day = '.\\day\\'

subfolder = "daylight_electric_screen"

dir_day = dir_day + subfolder + '\\'

fname_pre = "annual-output"

def day_index(month, day):
	index = 0
	
	days_in_mo = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
	for m in range(month - 1):
		index += days_in_mo[m]
	
	index += day - 1
	
	return index

# grab seasonal time indices
spring_indices = range(day_index(3,1), day_index(5,31) + 1, 1)
summer_indices = range(day_index(6,1), day_index(8,31) + 1, 1)
fall_indices = range(day_index(9,1), day_index(11,30) + 1, 1)
winter_indices = range(day_index(12,1), 365 + 1, 1)
winter_indices.extend( range(1, day_index(2,28) + 1, 1) )

class metrics:
	"""Contains our annual metrics for display / summary!"""
	
	def __init__(self, name):
		self.name = name
		
		self.X = 0.0
		self.Y = 0.0
		self.Z = 0.0
		self.dX = 0.0
		self.dY = 0.0
		self.dZ = 0.0
	
		self.morning_rt_spring = []
		self.morning_rt_summer = []
		self.morning_rt_fall = []
		self.morning_rt_winter = []
		self.afternoon_rt_spring = []
		self.afternoon_rt_summer = []
		self.afternoon_rt_fall = []
		self.afternoon_rt_winter = []
		self.evening_rt_spring = []
		self.evening_rt_summer = []
		self.evening_rt_fall = []
		self.evening_rt_winter = []
		self.rt_spring = []
		self.rt_summer = []
		self.rt_fall = []
		self.rt_winter = []

		self.morning_kss_spring = []
		self.morning_kss_summer = []
		self.morning_kss_fall = []
		self.morning_kss_winter = []
		self.afternoon_kss_spring = []
		self.afternoon_kss_summer = []
		self.afternoon_kss_fall = []
		self.afternoon_kss_winter = []
		self.evening_kss_spring = []
		self.evening_kss_summer = []
		self.evening_kss_fall = []
		self.evening_kss_winter = []
		self.kss_spring = []
		self.kss_summer = []
		self.kss_fall = []
		self.kss_winter = []

		self.mel_spring = []
		self.mel_summer = []
		self.mel_fall = []
		self.mel_winter = []
		self.morning_mel_spring = []
		self.morning_mel_summer = []
		self.morning_mel_fall = []
		self.morning_mel_winter = []
		self.evening_mel_spring = []
		self.evening_mel_summer = []
		self.evening_mel_fall = []
		self.evening_mel_winter = []
		
		self.I_spring = []
		self.I_summer = []
		self.I_fall = []
		self.I_winter = []
		self.morning_I_spring = []
		self.morning_I_summer = []
		self.morning_I_fall = []
		self.morning_I_winter = []
		self.afternoon_I_spring = []
		self.afternoon_I_summer = []
		self.afternoon_I_fall = []
		self.afternoon_I_winter = []
		self.evening_I_spring = []
		self.evening_I_summer = []
		self.evening_I_fall = []
		self.evening_I_winter = []

		self.ps_spring = []
		self.ps_summer = []
		self.ps_fall = []
		self.ps_winter = []
		
	def seasonal_average(self):
		self.morning_rt_spring = stats.mean(self.morning_rt_spring)
		self.morning_rt_summer = stats.mean(self.morning_rt_summer)
		self.morning_rt_fall = stats.mean(self.morning_rt_fall)
		self.morning_rt_winter = stats.mean(self.morning_rt_winter)
		self.afternoon_rt_spring = stats.mean(self.afternoon_rt_spring)
		self.afternoon_rt_summer = stats.mean(self.afternoon_rt_summer)
		self.afternoon_rt_fall = stats.mean(self.afternoon_rt_fall)
		self.afternoon_rt_winter = stats.mean(self.afternoon_rt_winter)
		self.evening_rt_spring = stats.mean(self.evening_rt_spring)
		self.evening_rt_summer = stats.mean(self.evening_rt_summer)
		self.evening_rt_fall = stats.mean(self.evening_rt_fall)
		self.evening_rt_winter = stats.mean(self.evening_rt_winter)
		self.rt_spring = stats.mean(self.rt_spring)
		self.rt_summer = stats.mean(self.rt_summer)
		self.rt_fall = stats.mean(self.rt_fall)
		self.rt_winter = stats.mean(self.rt_winter)
		
		self.morning_kss_spring = stats.mean(self.morning_kss_spring)
		self.morning_kss_summer = stats.mean(self.morning_kss_summer)
		self.morning_kss_fall = stats.mean(self.morning_kss_fall)
		self.morning_kss_winter = stats.mean(self.morning_kss_winter)
		self.afternoon_kss_spring = stats.mean(self.afternoon_kss_spring)
		self.afternoon_kss_summer = stats.mean(self.afternoon_kss_summer)
		self.afternoon_kss_fall = stats.mean(self.afternoon_kss_fall)
		self.afternoon_kss_winter = stats.mean(self.afternoon_kss_winter)
		self.evening_kss_spring = stats.mean(self.evening_kss_spring)
		self.evening_kss_summer = stats.mean(self.evening_kss_summer)
		self.evening_kss_fall = stats.mean(self.evening_kss_fall)
		self.evening_kss_winter = stats.mean(self.evening_kss_winter)
		self.kss_spring = stats.mean(self.kss_spring)
		self.kss_summer = stats.mean(self.kss_summer)
		self.kss_fall = stats.mean(self.kss_fall)
		self.kss_winter = stats.mean(self.kss_winter)
		
		self.mel_spring = stats.mean(self.mel_spring)
		self.mel_summer = stats.mean(self.mel_summer)
		self.mel_fall = stats.mean(self.mel_fall)
		self.mel_winter = stats.mean(self.mel_winter)
		self.morning_mel_spring = stats.mean(self.morning_mel_spring)
		self.morning_mel_summer = stats.mean(self.morning_mel_summer)
		self.morning_mel_fall = stats.mean(self.morning_mel_fall)
		self.morning_mel_winter = stats.mean(self.morning_mel_winter)
		self.evening_mel_spring = stats.mean(self.evening_mel_spring)
		self.evening_mel_summer = stats.mean(self.evening_mel_summer)
		self.evening_mel_fall = stats.mean(self.evening_mel_fall)
		self.evening_mel_winter = stats.mean(self.evening_mel_winter)
		
		self.ps_spring = sum(self.ps_spring) / 60.0
		self.ps_summer = sum(self.ps_summer) / 60.0
		self.ps_fall = sum(self.ps_fall) / 60.0 
		self.ps_winter = sum(self.ps_winter) / 60.0 
		
		self.morning_I_spring = stats.mean(self.morning_I_spring)
		self.morning_I_summer = stats.mean(self.morning_I_summer)
		self.morning_I_fall = stats.mean(self.morning_I_fall)
		self.morning_I_winter = stats.mean(self.morning_I_winter)
		self.afternoon_I_spring = stats.mean(self.afternoon_I_spring)
		self.afternoon_I_summer = stats.mean(self.afternoon_I_summer)
		self.afternoon_I_fall = stats.mean(self.afternoon_I_fall)
		self.afternoon_I_winter = stats.mean(self.afternoon_I_winter)
		self.evening_I_spring = stats.mean(self.evening_I_spring)
		self.evening_I_summer = stats.mean(self.evening_I_summer)
		self.evening_I_fall = stats.mean(self.evening_I_fall)
		self.evening_I_winter = stats.mean(self.evening_I_winter)
		self.I_spring = stats.mean(self.I_spring)
		self.I_summer = stats.mean(self.I_summer)
		self.I_fall = stats.mean(self.I_fall)
		self.I_winter = stats.mean(self.I_winter)


with open(fname_pre + '_seasonal-summary.csv', 'w') as f_w:
	f_w.write('view,season,x,y,z,dx,dy,dz,phaseshift,melsuppress_day,melsuppress_morning,melsuppress_evening,kss_day,kss_morning,kss_afternoon,kss_evening,vpvtrt_day,vpvtrt_morning,vpvtrt_afternoon,vpvtrt_evening,eemel_day,eemel_morning,eemel_afternoon,eemel_evening\n')

for view in views:
	vm = metrics(view) # metrics!
	
	fname = dir_day + fname_pre + "_%s.csv" % view
	
	with open(fname, 'r') as f:
		n_day = 0
		for n, line in enumerate(f):
			if n > 0:
				# X,Y,Z,dX,dY,dZ,day,t_crit_melpeak_s,t_crit_melpeak_h,t_crit_cbtmin_s,t_crit_cbtmin_h,phase_shift_s,melatonin_supression,KSS_8-12,KSS_12-18,KSS_18-24,KSS_day,vPVTRT_8-12,vPVTRT_12-18,vPVTRT_18-24,vPVTRT_day																							12
				parts = line.rstrip().split(',')
				
				if n == 1:
					X = float(parts[0])
					Y = float(parts[1])
					Z = float(parts[2])
					dX = float(parts[3])
					dY = float(parts[4])
					dZ = float(parts[5])
					
					metrics.X = X
					metrics.Y = Y
					metrics.Z = Z
					metrics.dX = dX
					metrics.dY = dY
					metrics.dZ = dZ
				
				day = int(parts[6])
				
				ps = float(parts[11])
				if (ps < -80000):
					ps = (ps + 24.0 * 60.0 * 60.0) * -1.0
				elif (ps > 80000):
					ps = (ps - 24.0 * 60.0 * 60.0) * -1.0
				elif (ps > 10000):
					ps = 0
				elif (ps < -10000):
					ps = 0
				
				mel = float(parts[12])
				mel_morning = float(parts[13])
				mel_evening = float(parts[14])
				
				KSS_morning = float(parts[15])
				KSS_afternoon = float(parts[16])
				KSS_evening = float(parts[17])
				KSS_day = float(parts[18])
				
				RT_morning = float(parts[19])
				RT_afternoon = float(parts[20])
				RT_evening = float(parts[21])
				RT_day = float(parts[22])
				
				I_morning = float(parts[23])
				I_afternoon = float(parts[24])
				I_evening = float(parts[25])
				I_day = float(parts[26])
				
				if day in spring_indices:
					vm.mel_spring.append(mel)
					vm.morning_mel_spring.append(mel_morning)
					vm.evening_mel_spring.append(mel_evening)
					
					vm.morning_kss_spring.append(KSS_morning)
					vm.afternoon_kss_spring.append(KSS_afternoon)
					vm.evening_kss_spring.append(KSS_evening)
					vm.kss_spring.append(KSS_day)
					
					vm.morning_rt_spring.append(RT_morning)
					vm.afternoon_rt_spring.append(RT_afternoon)
					vm.evening_rt_spring.append(RT_evening)
					vm.rt_spring.append(RT_day)
					
					vm.morning_I_spring.append(I_morning)
					vm.afternoon_I_spring.append(I_afternoon)
					vm.evening_I_spring.append(I_evening)
					vm.I_spring.append(I_day)
					
					vm.ps_spring.append(ps)
				
				if day in summer_indices:
					vm.mel_summer.append(mel)
					vm.morning_mel_summer.append(mel_morning)
					vm.evening_mel_summer.append(mel_evening)
					
					vm.morning_kss_summer.append(KSS_morning)
					vm.afternoon_kss_summer.append(KSS_afternoon)
					vm.evening_kss_summer.append(KSS_evening)
					vm.kss_summer.append(KSS_day)
					
					vm.morning_rt_summer.append(RT_morning)
					vm.afternoon_rt_summer.append(RT_afternoon)
					vm.evening_rt_summer.append(RT_evening)
					vm.rt_summer.append(RT_day)
					
					vm.morning_I_summer.append(I_morning)
					vm.afternoon_I_summer.append(I_afternoon)
					vm.evening_I_summer.append(I_evening)
					vm.I_summer.append(I_day)
					
					vm.ps_summer.append(ps)
				
				if day in fall_indices:
					vm.mel_fall.append(mel)
					vm.morning_mel_fall.append(mel_morning)
					vm.evening_mel_fall.append(mel_evening)
					
					vm.morning_kss_fall.append(KSS_morning)
					vm.afternoon_kss_fall.append(KSS_afternoon)
					vm.evening_kss_fall.append(KSS_evening)
					vm.kss_fall.append(KSS_day)
					
					vm.morning_rt_fall.append(RT_morning)
					vm.afternoon_rt_fall.append(RT_afternoon)
					vm.evening_rt_fall.append(RT_evening)
					vm.rt_fall.append(RT_day)
					
					vm.morning_I_fall.append(I_morning)
					vm.afternoon_I_fall.append(I_afternoon)
					vm.evening_I_fall.append(I_evening)
					vm.I_fall.append(I_day)
					
					vm.ps_fall.append(ps)
				
				if day in winter_indices:
					vm.mel_winter.append(mel)
					vm.morning_mel_winter.append(mel_morning)
					vm.evening_mel_winter.append(mel_evening)
					
					vm.morning_kss_winter.append(KSS_morning)
					vm.afternoon_kss_winter.append(KSS_afternoon)
					vm.evening_kss_winter.append(KSS_evening)
					vm.kss_winter.append(KSS_day)
					
					vm.morning_rt_winter.append(RT_morning)
					vm.afternoon_rt_winter.append(RT_afternoon)
					vm.evening_rt_winter.append(RT_evening)
					vm.rt_winter.append(RT_day)
					
					vm.morning_I_winter.append(I_morning)
					vm.afternoon_I_winter.append(I_afternoon)
					vm.evening_I_winter.append(I_evening)
					vm.I_winter.append(I_day)
					
					vm.ps_winter.append(ps)
	
	vm.seasonal_average()
	
	with open(fname_pre + '_seasonal-summary.csv', 'a') as f_w:
		# view,season,x,y,z,dx,dy,dz,phaseshift,melsuppress,kss_morning,kss_afternoon,kss_evening,vpvtrt_morning,vpvtrt_afternoon,vpvtrt_evening
		
		# spring
		f_w.write('%s,spring,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,' % (view, X, Y, Z, dX, dY, dZ))
		f_w.write('%.1f,%.2f,%.2f,%.2f,' % (vm.ps_spring, vm.mel_spring, vm.morning_mel_spring, vm.evening_mel_spring))
		f_w.write('%.1f,%.1f,%.1f,%.1f,' % (vm.kss_spring, vm.morning_kss_spring, vm.afternoon_kss_spring, vm.evening_kss_spring))
		f_w.write('%.1f,%.1f,%.1f,%.1f,' % (vm.rt_spring, vm.morning_rt_spring, vm.afternoon_rt_spring, vm.evening_rt_spring))
		f_w.write('%.5f,%.5f,%.5f,%.5f\n' % (vm.I_spring, vm.morning_I_spring, vm.afternoon_I_spring, vm.evening_I_spring))
		
		# summer
		f_w.write('%s,summer,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,' % (view, X, Y, Z, dX, dY, dZ))
		f_w.write('%.1f,%.2f,%.2f,%.2f,' % (vm.ps_summer, vm.mel_summer, vm.morning_mel_summer, vm.evening_mel_summer))
		f_w.write('%.1f,%.1f,%.1f,%.1f,' % (vm.kss_summer, vm.morning_kss_summer, vm.afternoon_kss_summer, vm.evening_kss_summer))
		f_w.write('%.1f,%.1f,%.1f,%.1f,' % (vm.rt_summer, vm.morning_rt_summer, vm.afternoon_rt_summer, vm.evening_rt_summer))
		f_w.write('%.5f,%.5f,%.5f,%.5f\n' % (vm.I_summer, vm.morning_I_summer, vm.afternoon_I_summer, vm.evening_I_summer))
		
		# fall
		f_w.write('%s,fall,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,' % (view, X, Y, Z, dX, dY, dZ))
		f_w.write('%.1f,%.2f,%.2f,%.2f,' % (vm.ps_fall, vm.mel_fall, vm.morning_mel_fall, vm.evening_mel_fall))
		f_w.write('%.1f,%.1f,%.1f,%.1f,' % (vm.kss_fall, vm.morning_kss_fall, vm.afternoon_kss_fall, vm.evening_kss_fall))
		f_w.write('%.1f,%.1f,%.1f,%.1f,' % (vm.rt_fall, vm.morning_rt_fall, vm.afternoon_rt_fall, vm.evening_rt_fall))
		f_w.write('%.5f,%.5f,%.5f,%.5f\n' % (vm.I_fall, vm.morning_I_fall, vm.afternoon_I_fall, vm.evening_I_fall))
		
		# winter
		f_w.write('%s,winter,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,' % (view, X, Y, Z, dX, dY, dZ))
		f_w.write('%.1f,%.2f,%.2f,%.2f,' % (vm.ps_winter, vm.mel_winter, vm.morning_mel_winter, vm.evening_mel_winter))
		f_w.write('%.1f,%.1f,%.1f,%.1f,' % (vm.kss_winter, vm.morning_kss_winter, vm.afternoon_kss_winter, vm.evening_kss_winter))
		f_w.write('%.1f,%.1f,%.1f,%.1f,' % (vm.rt_winter, vm.morning_rt_winter, vm.afternoon_rt_winter, vm.evening_rt_winter))	
		f_w.write('%.5f,%.5f,%.5f,%.5f\n' % (vm.I_winter, vm.morning_I_winter, vm.afternoon_I_winter, vm.evening_I_winter))	
		
		
		
		
		
		
		