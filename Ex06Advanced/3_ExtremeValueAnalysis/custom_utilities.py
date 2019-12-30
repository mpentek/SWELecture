#===============================================================================
'''
Project: Lecture - Structural Wind Engineering WS19-20 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        Various custom utilities for statistical evaluations

Author: mate.pentek@tum.de, anoop.kodakkal@tum.de
         
Description: This is a script containing necessary function definitions for examples 
        and tasks.

Note: ...

Created on:  30.11.2015
Last update: 27.09.2019
'''
#===============================================================================
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import gaussian_kde
import timeit

#===============================================================================

def get_pdf_kde(given_series):
    '''
    The function get_pdf_kde evaluates the probability distribution function (pdf)
    of the samples by using a non-parametric estimation technique called Kernal Desnity 
    Estimation (KDE). More details can be found at 
    https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.gaussian_kde.html.
    '''

    series_max = np.max(given_series)
    series_min = np.min(given_series)
    kde = gaussian_kde(given_series)
    pdf_x = np.linspace(series_min, series_max, 1000)
    pdf_y = kde(pdf_x)
    return pdf_x, pdf_y

def get_pdf_normal(given_series):
    '''
    The function get_pdf_normal estimates the normal pdf of the signal from the mean 
    and standard deviation of the samples. Recall the fact that a Normal distribution 
    can be entirely defined by two parameters, namely the mean and standard deviation. 
    More details about the function mlab.normpdf can be found at 
    https://matplotlib.org/api/mlab_api.html. 
    '''    

    series_max = np.max(given_series)
    series_std = np.std(given_series)
    series_m = np.mean(given_series)
    series_min = np.min(given_series)
    series_step = (series_max - series_min)/ 1000
    series_pdf = mlab.normpdf(np.arange(series_min, 
                                            series_max + series_step, 
                                            series_step), 
                                  series_m, 
                                  series_std)
    pdf_x = np.arange(series_min, series_max + series_step, series_step)
    pdf_y = series_pdf
    return pdf_x, pdf_y


def get_pdf_const(given_series):
    '''
    The function get_pdf_const mimcs the pdf of a constant signal which is a vertical
    line extending to infinity and the cdf is a vertical line to unity. 
    As 'KDE' and 'Normal' are unable to deliver this behavious use this function for 
    constant series   
    '''    

    series_std = np.std(given_series)
    series_m = np.mean(given_series)
    series_step = (2 * series_m - 0)/ 1000
    if series_std ==0.0:
        pdf_x = np.arange(0, 2 * series_m , series_step)
        pdf_y = np.zeros(len(pdf_x))
        pdf_y[int(len(pdf_x)/2)] = len(given_series)
    else : 
        raise Exception("The given series is not a Constant signal, use 'KDE' or 'Normal'")
    return pdf_x, pdf_y
    

def get_pdf(given_series, case='KDE'):
    if case == 'KDE':
        return get_pdf_kde(given_series)
    elif case == 'Normal':
        return get_pdf_normal(given_series)
    elif case == 'Constant':
        return get_pdf_const(given_series)
    else:
        raise NotImplementedError("PDF type not implemented, choose either KDE, Normal or Constant")
    
def get_fft(given_series, sampling_freq):
    '''
    The function get_fft estimates the Fast Fourier transform of the given signal 
    '''

    signal_length=len(given_series)

    freq_half =  np.arange(0, 
                           sampling_freq/2 - sampling_freq/signal_length + sampling_freq/signal_length, 
                           sampling_freq/signal_length)

    # single sided fourier
    series_fft = np.fft.fft(given_series)
    series_fft = np.abs(series_fft[0:int(np.floor(signal_length/2))])/np.floor(signal_length/2)  
    
    max_length = len(freq_half)
    if max_length < len(series_fft):
        max_length = len(series_fft)
    
    freq_half = freq_half[:max_length-1]
    series_fft = series_fft[:max_length-1]
    
    return freq_half, series_fft
    
def get_ecdf(series_pdf_x, series_pdf_y):
    '''
    The function get_ecdf computes the emperital CDF of the given PDF.
    '''
    
    # set up data
    dx = series_pdf_x[1] - series_pdf_x[0]    
    Y = series_pdf_y
    # normalize data
    Y /= (dx * series_pdf_y).sum()
    # compute ecdf
    CY = np.cumsum(Y * dx)
    
    return CY
    
def get_pot(given_series, threshold_value):
    '''
    The function get_pot computes the Peak Over Threshold for a given threshold value.
    '''
    pot_index = []
    pot_extreme = []

    for i in range(len(given_series)): 
         if  threshold_value < np.abs(given_series[i]): 
             pot_index.append(i)
             pot_extreme.append(np.abs(given_series[i]))
   
    return pot_index, pot_extreme
    
def get_bm(given_series, block_size):
    '''
    The function get_bm computes the Block Maxima of the signal for a given block size
    '''

    block_max = np.abs(given_series[0])
    
    bm_index = []
    bm_extreme = [] 

    for i in range(1,len(given_series)):
         if  block_max < np.abs(given_series[i]):
             block_max = np.abs(given_series[i])
             block_max_idx = i

         if  (i+1) % block_size == 0: 
             bm_index.append(block_max_idx)
             bm_extreme.append(block_max)
            
             block_max = 0
             
    return bm_index, bm_extreme
    
def get_pot_runtime(given_series, threshold_param):
    '''
    The function get_pot_runtime does a runtime evaluation of the Peak Over Threshold 
    for a given threshold parameter
    '''
    tic = timeit.default_timer()    
    
    # renaming
    values = given_series
    all_values = len(given_series)
    
    #setting up necessary variables
    meannew = values[0]
    meanold = meannew
    rmsnew = np.abs(values[0])
    rmsold = rmsnew
    part1 = values[0]*values[0]
    part2 = values[0]

    arraysorted = np.zeros(all_values) #preallocate with zeros
    arraysorted[0] = values[0] # create sorted array with new values
    
    pot_index = []
    pot_extreme = []
    
    res_m = np.zeros(all_values)
    res_m[0] = meannew
    res_rms = np.zeros(all_values)
    res_rms[0] = rmsnew

    res_std = np.zeros(all_values)
    
    res_med = np.zeros(all_values)
    res_med[0] = values[0]    
    
    res_skew = np.zeros(all_values)

    res_thres = np.zeros(all_values)
    res_thres[0] = meannew + threshold_param * 0
    
    
    # starting calculation loop
    for i in range(1,all_values): 
        # calculate mean
        meannew = (meanold * (i-1) + values[i])/i
        meanold = meannew
        
        # calculate rms
        rmsnew = np.sqrt(((i-1)*rmsold*rmsold + values[i]*values[i])/i)
        rmsold = rmsnew
        
        # calculate standard deviation, initally summing up parts
        part1 = part1 + values[i] * values[i]
        part2 = part2 + values[i]
        standarddev = np.sqrt((part1 - 2* meannew * part2 + meannew*meannew* i)/(i-1))
        
        
        # calculate median 
        if values[i]>= arraysorted[i-1]:
          arraysorted[i] = values[i]
        elif values[i]<= arraysorted[0]:
          arraysorted[1:i+1] = arraysorted[0:i]
          arraysorted[0] = values[i]
        else:
          j = i-1
          push_here = j
          while values[i] < arraysorted[j]:
            j = j - 1
            push_here = j
    
          arraysorted[push_here+1:i+1] = arraysorted[push_here:i]
          arraysorted[push_here+1] = values[i] 
    
        if (i % 2 == 1): # check if it matches matlab indexing
            medianv = (arraysorted[i//2] + arraysorted[i //2+1])/2                 
        else:
            medianv = arraysorted[i//2]
            
        # calculate skewness
        skewness = (meannew - medianv)/standarddev

        # threshold 6 sigma criterium
        threshold = (meannew + threshold_param * standarddev)
        if  threshold < np.abs(values[i]): 
             pot_index.append(i)
             pot_extreme.append(np.abs(values[i]))
                        
        # append results
        res_m[i] = meannew
        res_rms[i] = rmsnew
        res_std[i] = standarddev
        res_med[i] = medianv
        res_skew[i] = skewness
        res_thres[i] = threshold
        
    toc = timeit.default_timer()
    
    print('Elapsed time for get_pot_runtime function evaluation: ', toc-tic ,'s\n')

    return res_m, res_rms, res_std, res_med, res_skew, res_thres, pot_index, pot_extreme
    
#===============================================================================
## End functions
