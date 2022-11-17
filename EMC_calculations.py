import flow_calculations as fc
import pandas as pd
import numpy as np
import matplotlib as plt

# HACK - Hard-coded project parameters - Need to read from template
surface_area_sf = 1262
dim_width_ft = 20
dim_length_ft = 10
side_slope = .33
radius_in = 6
length_ft = 27.583
nMannings = np.arange(0.009, 0.033, 0.002)
slope = 0.002
gravity_SI = 9.81

# Hard-coded program parameters, including conversion factors
sf_to_m2 = 0.0929
in_to_m = 0.0254
ft_to_m = 0.3048
min_to_sec = 60.0
cfs_to_cms = 0.0283
Q_guess_factor = 2
Q_guess_interval = 1000


def read_template():
    global y1_depth, y2_depth, y_time, Qex_flow, Qex_time


    # Use pandas to read data template
    data = pd.read_excel(r'C:/Users/edwardt/SCCWRP/SMC BMP Projects - Regional monitoring network/Flow Calculations from Depth - M1 Backwater Curve/Flow Calculation from Depth.xlsx', sheet_name='For Python')
    df_depth = pd.DataFrame(data, columns=['Elapsed Time Depth (min)', 'y1 (in)', 'y2 (in)'])
    df_flow = pd.DataFrame(data, columns=['Elapsed Time Flow (min)', 'Exfiltration (cfs)'])
    df_result = pd.DataFrame(data, columns=['SFSH Q (cfs)', 'GVPF Q (cfs)'])

    # Strip out NA's from dataframes, store as 1D-list variables
    y1_depth = df_depth['y1 (in)'][df_depth['y1 (in)'].notna()]
    y2_depth = df_depth['y2 (in)'][df_depth['y2 (in)'].notna()]
    y_time = df_depth['Elapsed Time Depth (min)'][df_depth['Elapsed Time Depth (min)'].notna()]
    Qex_time = df_flow['Elapsed Time Flow (min)'][df_flow['Elapsed Time Flow (min)'].notna()]
    Qex_flow = df_flow['Exfiltration (cfs)'][df_flow['Exfiltration (cfs)'].notna()]

    return


def settings_SI():
    global surface_area_m2, radius_m, length_m, dim_width_m, dim_length_m, zz_drop_m, area_full_m2, wp_full_m, rh_full_m

    # Convert units to SI
    surface_area_m2 = surface_area_sf * sf_to_m2
    radius_m = radius_in * in_to_m
    length_m = length_ft * ft_to_m
    dim_width_m = dim_width_ft * ft_to_m
    dim_length_m = dim_length_ft * ft_to_m

    zz_drop_m = slope * length_m
    area_full_m2 = np.pi * radius_m**2
    wp_full_m = np.pi * radius_m * 2
    rh_full_m = area_full_m2/wp_full_m

    # Convert min -> sec, in -> m, cfs -> cms
    for ii in range(len(y_time)):
        y_time[ii] = y_time[ii] * min_to_sec
        y1_depth[ii] = y1_depth[ii] * cfs_to_cms
        y2_depth[ii] = y2_depth[ii] * cfs_to_cms

        # If instrument-recorded depths are negative, set them to zero
        if y1_depth[ii] < 0:
            y1_depth[ii] = 0
        if y2_depth[ii] < 0:
            y2_depth[ii] = 0

    for ii in range(len(Qex_time)):
        Qex_time[ii] = float(Qex_time[ii] * min_to_sec)
        Qex_flow[ii] = Qex_flow[ii] * cfs_to_cms

    return 


def main():
    read_template()
    settings_SI()

    Qin_SFSH, Q_exfil = fc.SandFilter_CV(y2_depth, y_time, Qex_flow, Qex_time)

    print(Qin_SFSH, Q_exfil)

    return

if __name__ == "__main__":
    main()



