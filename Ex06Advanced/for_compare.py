from bluecoefficients import bluecoeff

def blue4pressure(series, n_blocks, P1 = 0.80, P2 = 0.5704, dur_ratio = None, return_max = True, return_loc_scale = False):
    '''
    From a time series, blue4pressure estimates
    extremes of positive and negative values based on Lieblein's BLUE 
    (Best Linear Unbiased Estimate) method applied to n_blocks epochs. 
    
    Extremes are estimated for the duration of the record and for a ratio of it for probabilities of non-exceedance 
    P1 and P2 of the Gumbel distribution fitted to the epochal peaks.
    n = integer, dur need NOT be an integer.
    Written by Dat Duthinh 8_25_2015, 2_2_2016, 2_6_2017

    Reference: 
    1) Julius Lieblein "Efficient Methods of Extreme-Value
       Methodology" NBSIR 74-602 OCT 1974 for n = 4:16
    2) Nicholas John Cook "The designer's guide to wind loading of
       building structures" part 1, British Research Establishment 1985 Table C3
       pp. 321-323 for n = 17:24. Extension to n=100 by Adam Pintar Feb 12 2016.
    3) INTERNATIONAL STANDARD, ISO 4354 (2009-06-01), 2nd edition, “Wind 
       actions on structures,” Annex D (informative) “Aerodynamic pressure and 
       force coefficients,” Geneva, Switzerland, p. 22
    
    INPUT 
    series = vector of time history of pressure coefficients
    n = number of epochs (integer)of series data, 4 <= n <= 100
    dur = number of epochs for estimation of extremes. Default dur = n
          dur need not be an integer
    NOTE:
    replaced dur by dur_ratio to have the same as maxminest
    P1, P2 = probabilities of non-exceedance of extremes in EV1 (Gumbel).  
    P1 defaults to 0.80 (ISO) and P2 to 0.5704 (mean).

    OUTPUT 
    suffix max for + peaks, min for - peaks of pressure coeff.

    NOTE:
    changed: default returning extremes for the duration of the record.
             if dur_ratio is NOT given this is the return value.
             if dur_ratio is given the extremes for the duration of dur = dur_ratio * len(series) is returned 
             e.g. series is 2h long with dur_ratio = 0.5 the extreme within a period of 1h can be calculated. 
             => like this it is anlogous to NIST maxminest/qnt

    these are the returned values:

    p1_rmax (p1_rmin)= extreme value of positive (negative) peaks with probability of non-exceedance P1 for duration of series
    p2_rmax (p2_rmin)= extreme value of positive (negative) peaks with probability of non-exceedance P2 for for dur duration of series
        
    Computed but not returned are:

    p1_max (p1_min)= extreme value of positive (negative) peaks with probability of non-exceedance P1 for duration of 1 epoch
    p2_max (p2_min)= extreme value of positive (negative) peaks withprobability of exceedance P2 for duration of 1 epoch

    series_max (series_min)= vector of n positive (negative) epochal peaks
    u_max, b_max (u_min, b_min) = location and scale parameters of EV1
    (Gumbel) for positive (negative) peaks
    '''
    
    import numpy as np
    import math
    # Size of series array
    t = len(series)

    # Initialize variables
    series_max = np.zeros([n_blocks,1])
    series_min = np.zeros([n_blocks,1])

    # Find the peaks for each of the n user-defined epochs
    # and store in seriesmax and seriesmin arrays

    # Separate cases if n evenly divides t or not
    r = np.fmod(t,n_blocks)
    if r == 0:
        for i in np.arange(0,n_blocks):
            a = series[int(i*t/n_blocks):int((i+1)*t/n_blocks)]
            series_max[i] = a.max()
            series_min[i] = a.min()

    elif r > n_blocks/2:
        q = int(np.fix(t/n_blocks)+1)
        for i in np.arange(0,n_blocks-1):
            a = series[i*q:(i+1)*q]
            series_max[i] = a.max()
            series_min[i] = a.min()

    else:
        q = int(np.fix(t/n_blocks))
        for i in np.arange(0,n_blocks-1):
            a = series[i*q:(i+1)*q]
            series_max[i] = a.max()
            series_min[i] = a.min()
        
    # Coefficients for all n
    [ai,bi]= bluecoeff(n_blocks)

    # Organize values in ascending or descending order
    x_max = np.sort(series_max,axis=0)
    x_min = np.sort(series_min,axis=0)
    x_min = x_min[::-1]

    # defaults    
    if not dur_ratio:
        dur = n_blocks
    else:
        dur = dur_ratio * n_blocks
    # ************************** MAX CASE PEAK ***************************
    u = 0 # location parameter
    b = 0 # scale parameter

    # Calculate parameters of location and scale
    # Lieblein eq. 4
    for j in np.arange(0,n_blocks):
        u = u + ai[j]*x_max[j]
        b = b + bi[j]*x_max[j]
    
    p1_max = u - b*np.log(-np.log(P1)) # for 1 epoch
    p1_rmax = p1_max + b*np.log(dur) # for longer duration
    p2_max = u - b*np.log(-np.log(P2)) # for 1 epoch
    p2_rmax = p2_max + b*np.log(dur) # for longer duration
    u_max = u
    b_max = b
    # ************************** MIN CASE PEAK ***************************
    u = 0
    b = 0

    # Calculate parameters of location and scale
    for j in np.arange(0,n_blocks):
        u = u + ai[j]*x_min[j]
        b = b + bi[j]*x_min[j]

    p1_min = u - b*np.log(-np.log(P1))  # for 1 epoch
    p1_rmin = p1_min + b*np.log(dur)  # for longer duration
    p2_min = u - b*np.log(-np.log(P2))  # for 1 epoch
    p2_rmin = p2_min + b*np.log(dur)  # for longer duration
    u_min = u
    b_min = abs(b)
    
    #print(p1_max, p2_max, p1_rmax, p2_rmax, u_max, b_max, series_max, p1_min, p2_min, p1_rmin, p2_rmin, u_min, b_min, series_min)
    if return_max:
        return p1_rmax, p1_rmin, p2_rmax, p2_rmin
    if return_loc_scale:
        return (p1_rmax, p1_rmin, max_loc, max_scale, min_loc, min_scale, p2_rmax, p2_rmin)