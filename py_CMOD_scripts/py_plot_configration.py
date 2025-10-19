

def plot_prop_conf(mag_pairs_plot):
    
    mag_plot_prop_dict = {}
    
    # Plots configuration
    for mag_pair in mag_pairs_plot:
        for mag in mag_pair:
            
            mag_plot_style = 'sct'
            mag_legend_label = ' '
            mag_alpha = 0.75
            mag_zorder = 3
            mag_Y_ax_inversion = False
            
            if mag == 'Ax_rat':
                mag_label = 'Axis Ratio'
                mag_color = 'limegreen'
                mag_marker = '<'
                mag_marker_size = 150
                
            elif mag == 'Eff_rad_arcsec':
                mag_label = 'Effective Radius [arcsec]'
                mag_color = 'dodgerblue'
                mag_marker = '*'
                mag_marker_size = 200
                
            elif mag == 'Eff_rad_kpc':
                mag_label = 'Effective Radius [kpc]'
                mag_color = 'dodgerblue'
                mag_marker = '*'
                mag_marker_size = 200
            
            elif mag == 'Ind':
                mag_label = 'Sersic Index'
                mag_color = 'red'
                mag_marker = '>'
                mag_marker_size = 150
                
            elif mag == 'Mag':
                mag_label = 'Total Magnitude'
                mag_color = 'darkorange'
                mag_marker = 'o'
                mag_marker_size = 100
                mag_Y_ax_inversion = True
            
            elif mag == 'I_e':
                mag_label = 'Effective Intensity [counts]'
                mag_color = 'darkorange'
                mag_marker = 'o'
                mag_marker_size = 100    
            
            elif mag == 'Pos_ang':
                mag_label = 'Position Angle [°]'
                mag_color = 'gold'
                mag_marker = 'P'
                mag_marker_size = 150
                
            elif mag == 'X_cent':
                mag_label = 'Center X Position'
                mag_color = 'darkorchid'
                mag_marker = 'D'
                mag_marker_size = 100
                
            elif mag == 'Y_cent':
                mag_label = 'Center Y Position'
                mag_color = 'mediumturquoise'
                mag_marker = 'X'
                mag_marker_size = 150
            
            


            mag_plot_prop_dict[mag] = {}
            mag_plot_prop_dict[mag]['plot_style'] = mag_plot_style
            mag_plot_prop_dict[mag]['label'] = mag_label
            mag_plot_prop_dict[mag]['leg_lab'] = mag_legend_label
            mag_plot_prop_dict[mag]['color'] = mag_color
            mag_plot_prop_dict[mag]['marker'] = mag_marker
            mag_plot_prop_dict[mag]['marker_size'] = mag_marker_size
            mag_plot_prop_dict[mag]['alpha'] = mag_alpha
            mag_plot_prop_dict[mag]['zorder'] = mag_zorder
            mag_plot_prop_dict[mag]['Y_ax_inversion'] = mag_Y_ax_inversion
            
    return mag_plot_prop_dict