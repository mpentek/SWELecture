from bluecoeff import bluecoeff

def blue4pressure(cp, n, P1, P2, dur):
    # From a time series of pressure coefficients, blue4pressure estimates
    # extremes of positive and negative pressures based on Lieblein's BLUE 
    # (Best Linear Unbiased Estimate) method applied to n epochs. Extremes 
    # are estimated for 1 and dur epochs for probabilities of non-exceedance 
    # P1 and P2 of the Gumbel distribution fitted to the epochal peaks.
    # n = integer, dur need not be an integer.
    # Written by Dat Duthinh 8_25_2015, 2_2_2016, 2_6_2017
    # Reference: 1) Julius Lieblein "Efficient Methods of Extreme-Value
    # Methodology" NBSIR 74-602 OCT 1974 for n = 4:16
    # 2) Nicholas John Cook "The designer's guide to wind loading of
    # building structures" part 1, British Research Establishment 1985 Table C3
    # pp. 321-323 for n = 17:24. Extension to n=100 by Adam Pintar Feb 12 2016.
    # 3) INTERNATIONAL STANDARD, ISO 4354 (2009-06-01), 2nd edition, “Wind 
    # actions on structures,” Annex D (informative) “Aerodynamic pressure and 
    # force coefficients,” Geneva, Switzerland, p. 22

    # INPUT 
    # cp = vector of time history of pressure coefficients
    # n = number of epochs (integer)of cp data, 4 <= n <= 100
    # dur = number of epochs for estimation of extremes. Default dur = n
    # dur need not be an integer
    # P1, P2 = probabilities of non-exceedance of extremes in EV1 (Gumbel).  
    # P1 defaults to 0.80 (ISO)and P2 to 0.5704 (mean).
    # OUTPUT 
    # suffix max for + peaks, min for - peaks of pressure coeff.
    # p1_max (p1_min)= extreme value of positive (negative) peaks with
    # probability of non-exceedance P1 for 1 epoch
    # p2_max (p2_min)= extreme value of positive (negative) peaks with
    # probability of exceedance P2 for 1 epoch
    # p1_rmax (p1_rmin)= extreme value of positive (negative) peaks with
    # probability of non-exceedance P1 for dur epochs
    # p2_rmax (p2_rmin)= extreme value of positive (negative) peaks with
    # probability of non-exceedance P2 for for dur epochs
    # cp_max (cp_min)= vector of n positive (negative) epochal peaks
    # u_max, b_max (u_min, b_min) = location and scale parameters of EV1
    # (Gumbel) for positive (negative) peaks

    import numpy as np
    import math

    # Size of cp array
    t = len(cp)

    # Initialize variables
    cp_max = np.zeros([n,1])
    cp_min = np.zeros([n,1])

    # Find the peaks for each of the n user-defined epochs
    # and store in cpmax and cpmin arrays

    # Separate cases if n evenly divides t or not
    r = np.fmod(t,n)
    if r == 0:
        for i in np.arange(0,n):
            a = cp[int(i*t/n):int((i+1)*t/n)]
            cp_max[i] = a.max()
            cp_min[i] = a.min()
    elif r > n/2:
        q = np.fix(t/n)+1
        for i in np.arange(0,n-1):
            a = cp[int(i*q):int((i+1)*q)]
            cp_max[i] = a.max()
            cp_min[i] = a.min()

        a = cp[n*q:t]
        cp_max[n] = max[a]
        cp_min[n] = min[a]
    else:
        q = np.fix(t/n)
        for i in np.arange(0,n-1):
            a = cp[i*q:(i+1)*q]
            cp_max[i] = max[a]
            cp_min[i] = min[a]
        
        # a = cp(1+(n-1)*q:t)
        cp_max[n] = max[a]
        cp_min[n] = min[a]

    # Coefficients for all n
    [ai,bi]= bluecoeff(n)

    # Organize values in ascending or descending order
    x_max = np.sort(cp_max,axis=0)
    x_min = np.sort(cp_min,axis=0)
    x_min = x_min[::-1]

    if (P1==0):
        P1 = 0.80

    if (P2==0):
        P2 = 0.5704
    
    if (dur==0):
        dur = n
    
    # ************************** MAX CASE PEAK ***************************
    u = 0; # location parameter
    b = 0; # scale parameter

    # Calculate parameters of location and scale
    for j in np.arange(0,n):
        u = u + ai[j]*x_max[j]
        b = b + bi[j]*x_max[j]
    
    p1_max = u - b*math.log(-math.log(P1)); # for 1 epoch
    p1_rmax = p1_max + b*math.log(dur); # for longer duration
    p2_max = u - b*math.log(-math.log(P2)); # for 1 epoch
    p2_rmax = p2_max + b*math.log(dur); # for longer duration
    u_max = u
    b_max = b
    # ************************** MIN CASE PEAK ***************************
    u = 0
    b = 0

    # Calculate parameters of location and scale
    for j in np.arange(0,n):
        u = u + ai[j]*x_min[j]
        b = b + bi[j]*x_min[j]

    p1_min = u - b*math.log(-math.log(P1));  # for 1 epoch
    p1_rmin = p1_min + b*math.log(dur);  # for longer duration
    p2_min = u - b*math.log(-math.log(P2));  # for 1 epoch
    p2_rmin = p2_min + b*math.log(dur);  # for longer duration
    u_min = u
    b_min = abs(b)
    
    #print(p1_max, p2_max, p1_rmax, p2_rmax, u_max, b_max, cp_max, p1_min, p2_min, p1_rmin, p2_rmin, u_min, b_min, cp_min)
    return(p1_max, p2_max, p1_rmax, p2_rmax, u_max, b_max, cp_max, p1_min, p2_min, p1_rmin, p2_rmin, u_min, b_min, cp_min)