
import numpy as np
import sys

from scipy.optimize import curve_fit

from cmod.utils import round_number
from py_plot_general import plot_gen

def sersic_profile(r, n, R_e, I_e):

    b_n = 2 * n - 1 / 3 + 0.009876 / n
    I_r = I_e * np.exp(-b_n * ((r / R_e) ** (1 / n) - 1))
    
    return I_r

def sersic_profile_log(r, n, R_e, a):
    if 0.6<n<6 and 1<R_e<100 and 1<a<6.9:
        b_n = 2 * n - 1 / 3 + 0.009876 / n
        I_r_log = a + (-b_n * ((r / R_e) ** (1 / n) - 1))
        
    elif 0.6<n<6 and 1<R_e<100:
        b_n = 2 * n - 1 / 3 + 0.009876 / n
        I_r_log = a + (-b_n * ((r / R_e) ** (1 / n) - 1))

    else:
        I_r_log = 10e10

    return I_r_log

def exponential_disk(r, I_0, h):
    I_r = I_0 * np.exp(-r / h)
    return I_r

def exponential_disk_linear(x, m, b):
    y = m * x + b
    return y

def two_line_model(x, a1, b1, a2, b2, x_break):
    
    return np.where(x < x_break, a1 * x + b1, a2 * x + b2)

def calculate_chi_squared(x, y_obs, y_err, model, params):

    y_model = model(x, *params)  
    residuals = (y_obs - y_model) / y_err  
    chi_squared = np.sum(residuals**2)  
    return chi_squared


def initial_conditions(img_name,df, x_col, y_col,y_err_col,const,output_path):
    
    print('Computing the initial conditions')
    
    # Loading original data
    x_data = df[x_col].values
    y_data = df[y_col].values
    
    # FIRST STEP: Fitting two straight lines
    # The break point is computed automatically
    print(f'\nFittig two straight lines')
    
    # Converting original data into ln to work with
    # straight lines for a disk like part
    x_data_two = x_data
    y_data_two = np.log(df[y_col].values)
    y_err_two = np.log(df[y_err_col].values)
    
    # Initial estimation for parameters: a1, b1, a2, b2, x_break
    initial_guess_two = [1, 0, -1, 0, max(x_data_two) * (1/3)]
    
    # Fitting the two straight lines model
    popt_two, pcov_two = curve_fit(two_line_model, x_data_two, y_data_two, 
                                   p0=initial_guess_two,sigma=y_err_two)
    
    a1, b1, a2, b2, x_break_two = popt_two
    
    # Computing the ordinary parameter form
    I_0_two = np.exp(b2)
    h_two = -1 / a2
    
    print(f'I_0: {I_0_two:.3f}')
    print(f'h: {h_two:.3f}')
    print(f'X break point: {int(x_break_two)}')

    chi_squared_two = calculate_chi_squared(x_data_two, y_data_two, y_err_two, two_line_model, popt_two)
    print(f"Chi^2: {chi_squared_two:.3f}")
    
    # Obtainting the break position between the two lines
    break_pos = list(x_data_two).index(min(x_data_two[x_data_two >= x_break_two]))
    print(f'X break point pos: {break_pos}')

    # Representing the lines
    y_fit_line_1 = a1 * x_data_two[x_data_two < x_break_two] + b1
    y_fit_line_2 = a2 * x_data_two[x_data_two >= x_break_two] + b2
    
    y_fit_line_1 = list(y_fit_line_1) + [np.nan] * (len(x_data_two) - len(y_fit_line_1))
    y_fit_line_2 = [np.nan] * (len(x_data_two) - len(y_fit_line_2)) + list(y_fit_line_2)
    
    x_data_fig = [x_data,x_data,x_data]
    y_data_fig = [y_data_two,y_fit_line_1,y_fit_line_2]
    
    # Plotting the lines
    fig_name = f'{img_name}_IC_00'
    label_list = ['Intensity Profile','Line 1', 'Line 2']
    
    plot_gen(x_data=x_data_fig,
             y_data=y_data_fig,
             label_list=label_list,
             fig_name=fig_name,
             fig_save_path=output_path,
             line_color = ['gold','dodgerblue','black'],
             x_axis_label = ['$X\\, dimension\\, [\\mathrm{{pix}}]$'],
             x_ticks_dec = 0,
             y_axis_label = ['$\\ln(Intensity/\\mathrm{{counts}})$'],
             y_ticks_dec = 2,
             width_style = [0.75,2,2],
             zorder = [1,2,2],
             style_plot = ['sct','line','line'],
             plot_show= False)
    
    
    # SECOND STEP: disk fitting
    print(f'\nDisk-like fitting')
    
    # The disk-like fitting will be carry out with the
    # external part from the break position
    # Transforming data to ln
    x_data_disk = df[x_col].values[break_pos:]
    y_data_disk = np.log(df[y_col].values[break_pos:])
    y_err_disk = np.log(df[y_err_col].values[break_pos:])

    initial_guess_disk = [-1, np.mean(x_data_disk)]
    
    # Fitting the model to an exponential converted to
    # a straight line
    popt_disk, pcov_disk = curve_fit(exponential_disk_linear, x_data_disk, y_data_disk, 
                                     p0=initial_guess_disk,sigma=y_err_disk,nan_policy='omit')
    m_disk, b_disk = popt_disk
    
    # Changing from the straight line to
    # the exponential parameters
    I_0_disk = np.exp(b_disk)
    h_disk = -1/m_disk

    print(f'I_0: {I_0_disk:.3f}')
    print(f'h: {h_disk:.3f}')
        
    chi_squared_disk = calculate_chi_squared(x_data_disk, y_data_disk, y_err_disk, exponential_disk_linear, popt_disk)
    print(f"Chi^2 disk: {chi_squared_disk:.3f}")
    
    # Creating the exponential 
    y_fit_disk = exponential_disk_linear(x_data_disk, m_disk, b_disk)
    y_fit_disk = exponential_disk(x_data_disk,I_0_disk,h_disk)
    
    
    # THIRD STEP: bulge fitting
    print(f'\nBulge-like fitting')

    # For this fit an exponential is applied
    # affecting only to the intensity
    break_pos_bul = int(break_pos)
    x_data_bul = df[x_col].values[:break_pos_bul]
    y_data_bul = np.log(df[y_col].values[:break_pos_bul])
    y_err_bul = np.log(df[y_err_col].values[:break_pos_bul])
    
    if '-hr' in sys.argv:
        
        x_data_bul = [val for val in x_data_bul if list(x_data_bul).index(val)%2 == 0]
        y_data_bul = [val for val in y_data_bul if list(y_data_bul).index(val)%2 == 0]
        y_err_bul = [val for val in y_err_bul if list(y_err_bul).index(val)%2 == 0]

    initial_guess_bul = [3,40,3]
    
    # Fitting the model
    popt_bul, pcov_bul = curve_fit(sersic_profile_log, x_data_bul, y_data_bul, 
                                   p0=initial_guess_bul,sigma=y_err_bul)
    n_bul, r_e_bul, a = popt_bul
    
    # Computing the ordinary parameter form
    I_e_bul = np.exp(a)
    
    print(f'n: {n_bul:.3f}')
    print(f'R_e: {r_e_bul:.3f}')
    print(f'I_e: {I_e_bul:.3f}')
    chi_squared_bul = calculate_chi_squared(x_data_bul, y_data_bul, y_err_bul, sersic_profile_log, popt_bul)
    print(f'Chi^2: {chi_squared_bul:.3f}')
    
    # Crear los ajustes
    y_fit_bul = sersic_profile(x_data_bul, n_bul, r_e_bul, I_e_bul)
    
    # Ploting the fitting
    
    y_fit_bul = np.array(list(y_fit_bul) + [np.nan] * (len(x_data_two) - len(y_fit_bul)))
    y_fit_disk = np.array([np.nan] * (len(x_data_two) - len(y_fit_disk)) + list(y_fit_disk))
    
    x_data_fig = [x_data,x_data,x_data]
    y_data_fig = [y_data,y_fit_bul,y_fit_disk]
    
    label_list = ['Intensity Profile','Sersic Fitting - Bulge-like', 'Exponential Fitting - Disk-like']
    fig_name = f'{img_name}_IC_01'    
    
    plot_gen(x_data=x_data_fig,
             y_data=y_data_fig,
             label_list=label_list,
             fig_name=fig_name,
             fig_save_path=output_path,
             line_color = ['gold','dodgerblue','black'],
             x_axis_label =['$X\\, dimension\\, [\\mathrm{{pix}}]$'],
             x_ticks_dec = 0,
             y_axis_label = ['$Intensisty\\, [\mathrm{{counts}}]$'],
             y_ticks_dec = 2,
             width_style = [0.75,2,2],
             zorder = [1,2,2],
             style_plot = ['sct','line','line'],
             plot_show=False)
    
    return (break_pos, break_pos_bul, round_number(I_0_disk,3), round_number(h_disk,3), 
            round_number(n_bul,3), round_number(r_e_bul,3), round_number(I_e_bul,3))