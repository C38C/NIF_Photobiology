import datetime


def day_from_ill_datetime(E_list, datetime_list, sunrise_dt, sunset_dt, timestep=20):
	"""To deal with interpolation into an even timestep for use in biological_model code implementation of Postnova et al. / 2018. 
	
	timestep input is in seconds!"""
	
	E_interp = [0.0] * int((24 * 60 * 60) / timestep) # create blank array!
	
	current_time = datetime.datetime(year=2020, month=sunrise_dt.month, day=sunrise_dt.day, hour = 0, minute = 0)
	
	timestep_td = datetime.timedelta(seconds = timestep)
	
	for n in range( len(E_interp) ): # loop through all timesteps...
		if (current_time > sunrise_dt) and (current_time < sunset_dt):
			
			# find first datetime in datetime_list that you are larger than
			if current_time < datetime_list[0]:
				dt_last = sunrise_dt
				E_last = 0.0
				dt_next = datetime_list[0]
				E_next = E_list[0]
			elif (current_time >= datetime_list[0]) and (current_time <= datetime_list[-1]):
				# find location within!
				for index in range( len(datetime_list) ):
					if (current_time > datetime_list[index]) and (current_time < datetime_list[index + 1]):
						dt_last = datetime_list[index]
						E_last = E_list[index]
						dt_next = datetime_list[index + 1]
						E_next = E_list[index + 1]
						break
					if (current_time == datetime_list[index]):
						dt_last = datetime_list[index]
						dt_next = datetime_list[index]
						E_last = E_list[index]
						E_next = E_list[index]
						break
				
			else:
				dt_last = datetime_list[-1]
				E_last = E_list[-1]
				dt_next = sunset_dt
				E_next = 0.0
			
			# calculate linear scale values	
			if abs(E_last - E_next) < 0.001:
				this_E = E_last
			else:
				dist_last_td = current_time - dt_last
				dist_last_s = dist_last_td.seconds
				
				dist_next_ts = dt_next - current_time
				dist_next_s = dist_next_ts.seconds
				
				timestep_dist_td = dt_next - dt_last
				timestep_dist_s = timestep_dist_td.seconds
				
				this_E = (float(dist_next_s) / float(timestep_dist_s)) * E_last + (float(dist_last_s) / float(timestep_dist_s)) * E_next
			
			E_interp[n] = this_E
			
		current_time += timestep_td # increase current_time by timestep
	
	return E_interp
		
