# import flow_backwater as fbw
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from bisect import bisect_left

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
comp_vol = 1000

# Hard-coded program parameters, including conversion factors
sf_to_m2 = 0.0929
in_to_m = 0.0254
ft_to_m = 0.3048
min_to_sec = 60.0
cfs_to_cms = 0.0283
Q_guess_factor = 2
Q_guess_interval = 1000

# HACK - flag for whether backwater flow calculation is needed
fbw_needed = False


def read_template():
    global flow_time, sample_time, flow_value


    # Use pandas to read data template
    data = pd.read_excel(r'C:/Users/edwardt/SCCWRP/SMC BMP Projects - Regional monitoring network/Flow Calculations from Depth - M1 Backwater Curve/EMC_calculator_repo/Data Product QA/EMC_Template.xlsx', 
        sheet_name=r'Auckland Wetland Outlet')
    df = pd.DataFrame(data, columns=['Flow Elapsed Time (min)', 'Flow (m3/s)', 'Sample Elapsed Time (min)'])

    # Strip out NA's from dataframes, store as 1D-list variables
    flow_time = df['Flow Elapsed Time (min)'][df['Flow Elapsed Time (min)'].notna()]
    sample_time = df['Sample Elapsed Time (min)'][df['Sample Elapsed Time (min)'].notna()]
    flow_value = df['Flow (m3/s)'][df['Flow (m3/s)'].notna()]

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
    for ii in range(len(flow_time)):
        flow_time[ii] = int(flow_time[ii] * min_to_sec)
        flow_value[ii] = flow_value[ii] # * cfs_to_cms
    
    for ii in range(len(sample_time)):
        sample_time[ii] = int(sample_time[ii] * min_to_sec)

    return 


def volume_calc(y_series, x_series):
    """
    This subroutine calculates the trapezoidal approximation of the area under the curve
    """

    volume = np.trapz(y_series, x_series)
    return volume


def match_sample_to_inflow(sample_time, flow_time, flow_value):
    """
    This subroutine is used to map hydrograph data onto the sample times - this is for plotting purposes
    """
    sample_value = []
    sample_time = list(map(int, sample_time))
    flow_time = list(map(int, flow_time))
    flow_value = list(flow_value)


    # When the sample timeseries elapsed time is aligned with the flow timeseries elapsed time, include value 
    for ii in range(len(flow_time)):
        if flow_time[ii] in sample_time:
            sample_value.append(flow_value[ii])

    return sample_value


def volume_alloc(sample_time, flow_time, flow_value, comp_vol):
    """
    This subroutine is the main EMC calculator

    Loop through the samples
        assign hydrograph segments to each sample
        calculate the relative volume associated with that segment
    """

    # Recast as lists for easier indexing
    sample_time = list(sample_time)
    sample_volumes = []
    flow_time = list(flow_time)
    flow_value = list(flow_value)
    sample_time = [0] + sample_time + [flow_time[-1]]

    # # Hard-code for matching Auckland Wetland data better
    # sample_leftedge = 8*60
    # time_index_left = 4
    # sample_rightedge = np.average([sample_time[1], sample_time[2]])
    # time_index_right = flow_time.index(take_closest(flow_time, sample_rightedge))

    # Generalized to clean data

    # Interpolate between sample points to draw time-boundaries on hydrograph
    sample_leftedge = 0
    time_index_left = 0
    sample_rightedge = np.average([sample_time[1], sample_time[2]])
    time_index_right = flow_time.index(take_closest(flow_time, sample_rightedge))

    print(flow_time[time_index_left]/60, flow_time[time_index_right]/60)

    # Save the sample-assigned hydrograph time/values to variables
    flow_value_snip = flow_value[time_index_left:time_index_right+1]
    flow_time_snip = flow_time[time_index_left:time_index_right+1]

    # Calculate the total volume in that hydrograph snip
    volume_snip = volume_calc(flow_value_snip, flow_time_snip)
    sample_volumes.append(volume_snip)

    # Left and right-most boundaries are not included in loop because they must snap to the ends of the timeseries
    for ii in range(2, len(sample_time)-2):

        sample_leftedge = np.average([sample_time[ii-1], sample_time[ii]])
        time_index_left = flow_time.index(take_closest(flow_time, sample_leftedge))

        sample_rightedge = np.average([sample_time[ii], sample_time[ii+1]])
        time_index_right = flow_time.index(take_closest(flow_time, sample_rightedge))

        print(flow_time[time_index_left]/60, flow_time[time_index_right]/60)
        
        flow_value_snip = flow_value[time_index_left:time_index_right+1]
        flow_time_snip = flow_time[time_index_left:time_index_right+1]
        volume_snip = volume_calc(flow_value_snip, flow_time_snip)
        sample_volumes.append(volume_snip)


    # # Hard-code for matching Auckland Wetland better
    # sample_leftedge = sample_rightedge
    # time_index_left = flow_time.index(take_closest(flow_time, sample_leftedge))

    # sample_rightedge = 1502.0*60
    # time_index_right = flow_time.index(sample_rightedge)

    # print(flow_time[time_index_left]/60, flow_time[time_index_right]/60)

    # flow_value_snip = flow_value[time_index_left:time_index_right+1]
    # flow_time_snip = flow_time[time_index_left:time_index_right+1]
    # volume_snip = volume_calc(flow_value_snip, flow_time_snip)
    # sample_volumes.append(volume_snip)

    sample_leftedge = np.average([sample_time[-2], sample_time[-1]])
    time_index_left = flow_time.index(take_closest(flow_time, sample_leftedge))

    print(flow_time[time_index_left]/60, flow_time[time_index_right]/60)

    flow_value_snip = flow_value[time_index_left:]
    flow_time_snip = flow_time[time_index_left:]
    volume_snip = volume_calc(flow_value_snip, flow_time_snip)
    sample_volumes.append(volume_snip)

    vol_total = sum(sample_volumes)
    volume_allocated = np.empty(len(sample_volumes))
    for ss in range(len(sample_volumes)):
        volume_allocated[ss] = comp_vol * sample_volumes[ss] / vol_total
    return volume_allocated


def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before



def main():
    """
    This subroutine drives the EMC calculator
    """

    # Read and convert data to SI units
    read_template()
    settings_SI()

    # Map sample times to hydrograph
    sample_value = match_sample_to_inflow(sample_time, flow_time, flow_value)
    last_sample = list(sample_time)[-1]

    # Calculate the EMC composite volumes in the units given by comp_vol
    volume = volume_alloc(sample_time, flow_time, flow_value, comp_vol)

    print(len(sample_time))
    print(len(volume))
    round_volume = [round(num) for num in list(volume)]
    print("EMC Compositing Volumes (mL)")
    print(round_volume)
    print("")
    print(sum(volume))

    # Create hydrograph plot with samples overlaid
    fig, ax = plt.subplots()
    fig.suptitle('Wetland Outlet Hydrograph - Python')
    ax.plot(flow_time, flow_value, color='Blue', label='Flow')
    ax.plot(sample_time, sample_value, 'o', color='Green', label='Sample')
    ax.set(xlabel="Time since start (s)", ylabel="Flowrate (cms)")
    ax.legend()
    plt.ylim([0, 0.012])
    plt.xlim([0, last_sample])
    plt.show()

    return

if __name__ == "__main__":
    main()



