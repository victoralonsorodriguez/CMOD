

# This file creates the required constraints template for Galfit
def create_constraints(file,n_range,re_range,mag_range,xy_range,q_range,pa_range):
	
	file.write('# Component   parameter   constraint  (Optional and written with #)Commentary\n')
	file.write(f'# operation values \n')
	file.write(f'1 n {n_range[0]} to {n_range[1]} #polis email correction \n')
	file.write(f'1 re {re_range[0]} to {re_range[1]}  \n')
	file.write(f'1 mag {mag_range[0]} {mag_range[1]}  \n')
	file.write(f'1 x {xy_range[0]} {xy_range[1]}  \n')
	file.write(f'1 y {xy_range[0]} {xy_range[1]}  \n')
	file.write(f'1 q {q_range[0]} to {q_range[1]} \n')
	file.write(f'1 pa {pa_range[0]} to {pa_range[1]}')
