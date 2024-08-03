

# This file creates the required constraints template for Galfit
def create_initiaL_params(file,galaxy,version,eff_rad,ser_ind,ax_rat,pos_ang):
	
	file.write(f'# {galaxy} initial parameters introduced manually for {version}\n')
	file.write(f'Effective radius: {eff_rad} [pix]\n')
	file.write(f'Sérsic Index: {ser_ind}  \n')
	file.write(f'Axis ratio: {ax_rat}  \n')
	file.write(f'Position angle: {pos_ang} [degrees, up=0, left=90]  \n')
