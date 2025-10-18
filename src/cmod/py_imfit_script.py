
# This file creates the required template for Galfit
def create_conf_imfit(file,funct,galaxy,img_center_x,img_center_y,
                     pa_ser=[None,None,None],ell_ser=[None,None,None],n=[None,None,None],I_e=[None,None,None],r_e=[None,None,None], # Sersic
                     pa_disk=[None,None,None],ell_disk=[None,None,None],I_0=[None,None,None],h=[None,None,None],                      # Exponential
                     pa_bar=[None,None,None],ell_bar=[None,None,None],n_bar=[None,None,None],a_bar=[None,None,None],c0=[None,None,None]): # Ferrer
	
    if len(funct) != 0:
        
        file.write(f'# Configuration file for the galaxy {galaxy}\n')
        file.write(f'X0    {img_center_x[0]}    {img_center_x[1]},{img_center_x[2]}\n')
        file.write(f'Y0    {img_center_y[0]}    {img_center_y[1]},{img_center_y[2]}\n')

        for f in funct:
            if f == 'Sersic':
                file.write(f'FUNCTION {f}\n')
                file.write(f'PA    {pa_ser[0]}    {pa_ser[1]},{pa_ser[2]}\n')
                file.write(f'ell    {ell_ser[0]}    {ell_ser[1]},{ell_ser[2]}\n')
                file.write(f'n    {n[0]}    {n[1]},{n[2]}\n')
                file.write(f'I_e    {I_e[0]}    {I_e[1]},{I_e[2]}\n')
                file.write(f'r_e    {r_e[0]}    {r_e[1]},{r_e[2]}\n')

            elif f == 'Exponential':
                file.write(f'FUNCTION {f}\n')
                file.write(f'PA    {pa_disk[0]}    {pa_disk[1]},{pa_disk[2]}\n')
                file.write(f'ell    {ell_disk[0]}    {ell_disk[1]},{ell_disk[2]}\n')
                file.write(f'I_0    {I_0[0]}    {I_0[1]},{I_0[2]}\n')
                file.write(f'h    {h[0]}    {h[1]},{h[2]}\n')
                
            elif f == 'FerrersBar2D':
                file.write(f'FUNCTION {f}\n')
                file.write(f'PA    {pa_bar[0]}    {pa_bar[1]},{pa_bar[2]}\n')
                file.write(f'ell    {ell_bar[0]}    {ell_bar[1]},{ell_bar[2]}\n')
                file.write(f'I_0    {I_0[0]}    {I_0[1]},{I_0[2]}\n')
                file.write(f'n    {n_bar[0]}    {n_bar[1]},{n_bar[2]}\n')
                file.write(f'a_bar    {a_bar[0]}    {a_bar[1]},{a_bar[2]}\n')
                file.write(f'c0    {c0[0]}    {c0[1]},{c0[2]}\n')

    else:
        print('No functions for the fitting were selected')
        return None
	
	
