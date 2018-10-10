#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS18-19 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        Statistical evaluation of a generated signal

Author: mate.pentek@tum.de, kodakkal.anoop@tum.de 
         
Description: This is a script for run-time evaluation and maxima extraction.
        Examples and methodology based upon various PhD projects and current research.

Proposed exercises: See Jupyter Notebook file

Note:   ...
        Module dependencies:
            python
            matplotlib
            numpy
            custom_utilities.py

Created on:  30.11.2015
Last update: 30.10.2017
'''
#===============================================================================
## Display which python version is used
print('Python version:') 
import sys
print(sys.version)

#===================== BLOCK 1 =======================
# import python modules
import numpy as np
import scipy
from matplotlib import pyplot as plt

# import own modules
import custom_utilities as c_utils

#===================== BLOCK 2 =======================
# start time
start_time = 0.0
# end time
end_time = 10.0
# steps 
n_steps = 10000
# time step
delta_time = end_time / (n_steps-1)
# time series
# generate grid size vector (array) 1D
time_series = np.arange(start_time, end_time + delta_time, delta_time)

#===================== BLOCK 3 =======================
# frequency of the cosine
cos_freq = 10
# amplitude of the cosine
cos_ampl = 1
# series of the cosine
cos_series = cos_ampl * np.cos( 2*np.pi* cos_freq * time_series)

#===================== BLOCK 4 =======================
plt.figure(num=1, figsize=(15, 4))
plt.plot(time_series, cos_series)
plt.ylabel('Amplitude')
plt.xlabel('Time [s]')
plt.title('1. Cosine signal')
plt.grid(True)

#===================== BLOCK 5 =======================
# amplitude of the constant
const_ampl = 10
# series of the constant
const_series = const_ampl * np.ones(len(time_series))

#===================== BLOCK 6 =======================
plt.figure(num=2, figsize=(15, 4))
plt.plot(time_series, const_series)
plt.ylabel('Amplitude')
plt.xlabel('Time [s]')
plt.title('2. Constant signal')
plt.grid(True)

#===================== BLOCK 7 =======================
# random signal 
# assuming nomarl distribution
# with given mean m = 0 and standard deviation std = 0.25
rand_m = 0.0
rand_std = 0.25
# series of the random
rand_series = np.random.normal(rand_m, rand_std, len(time_series))

#===================== BLOCK 8 =======================
plt.figure(num=3, figsize=(15, 4))
plt.plot(time_series, rand_series)
plt.ylabel('Amplitude')
plt.xlabel('Time [s]')
plt.title('3. Random signal')
plt.grid(True)

#===================== BLOCK 9 =======================
#rand_series = np.random.lognormal(0, 0.25, len(time_series))
#rand_series = np.random.beta(1, 0.25, len(time_series))
#rand_series = np.random.rand(len(time_series))
#rand_series = np.random.uniform(0,1,len(time_series))

#===================== BLOCK 10 ======================
# coefs -> weighting factors for the respective series of signals
coef_signal1 = 1
coef_signal2 = 0.25
coef_signal3 = 1
superposed_series = coef_signal1 * const_series + coef_signal2 * cos_series + coef_signal3 * rand_series

#===================== BLOCK 11 ======================
plt.figure(num=4, figsize=(15, 4))
plt.plot(time_series, superposed_series)
plt.ylabel('Amplitude')
plt.xlabel('Time [s]')
plt.title('4. Superposed signal')
plt.grid(True)

#===================== BLOCK 12 ======================
# computing statistical quantitites (scalar values) and "converting" to an array for later plotting
cos_series_m = np.mean(cos_series) * np.ones(len(time_series))
cos_series_std = np.std(cos_series) * np.ones(len(time_series))
cos_series_rms = np.sqrt(np.mean(np.square(cos_series)))  * np.ones(len(time_series))
# printing statistical quantitites (scalar values) to the console
print('Mean: ', np.mean(cos_series))
print('STD: ', np.std(cos_series))
print('RMS: ', np.sqrt(np.mean(np.square(cos_series))))
print('Median: ', np.median(cos_series))
print('Skewness: ',(np.mean(cos_series) - np.median(cos_series))/np.std(cos_series))

#===================== BLOCK 13 ======================
const_series_m = np.mean(const_series) * np.ones(len(time_series))
const_series_std = np.std(const_series) * np.ones(len(time_series))
const_series_rms = np.sqrt(np.mean(np.square(const_series))) * np.ones(len(time_series))
print('Mean: ', np.mean(const_series))
print('STD: ', np.std(const_series))
print('RMS: ', np.sqrt(np.mean(np.square(const_series))))
print('Median: ', np.median(const_series))
print('Skewness: ', (np.mean(const_series) - np.median(const_series))/np.std(const_series))

#===================== BLOCK 14 ======================
rand_series_m = np.mean(rand_series) * np.ones(len(time_series))
rand_series_std = np.std(rand_series) * np.ones(len(time_series))
rand_series_rms = np.sqrt(np.mean(np.square(rand_series))) * np.ones(len(time_series))
print('Mean: ', np.mean(rand_series))
print('STD: ', np.std(rand_series))
print('RMS: ', np.sqrt(np.mean(np.square(rand_series))))
print('Median: ', np.median(rand_series))
print('Skewness: ', (np.mean(rand_series) - np.median(rand_series))/np.std(rand_series))

#===================== BLOCK 15 ======================
superposed_series_m = np.mean(superposed_series) * np.ones(len(time_series))
superposed_series_std = np.std(superposed_series) * np.ones(len(time_series))
superposed_series_rms = np.sqrt(np.mean(np.square(superposed_series))) * np.ones(len(time_series))
print('Mean: ', np.mean(superposed_series))
print('STD: ', np.std(superposed_series))
print('RMS: ', np.sqrt(np.mean(np.square(superposed_series))))
print('Median: ', np.median(superposed_series))
print('Skewness: ', (np.mean(superposed_series) - np.median(superposed_series))/np.std(superposed_series))

#===================== BLOCK 16 ======================
plt.rcParams["figure.figsize"] = (15,15)
plt.figure(num=5)

# const
plt.subplot(4, 2, 1)
plt.plot(time_series, const_series,
         time_series, const_series_m,
         time_series, const_series_m - const_series_std,
         time_series, const_series_m + const_series_std,
         time_series, const_series_rms)
plt.ylabel('Amplitude')
plt.title('Time series')
plt.grid(True)

bins = 100
plt.subplot(4, 2, 2)
plt.hist(const_series, bins)
plt.title('Histogram of ' + str(n_steps) +' values')
plt.ylabel('Frequency of oc.')
plt.grid(True)

# cos
plt.subplot(4, 2, 3)
plt.plot(time_series, cos_series) 
plt.plot(time_series, cos_series_m, label = 'Mean') 
plt.plot(time_series, cos_series_m - cos_series_std, label = 'Mean-SD') 
plt.plot(time_series, cos_series_m + cos_series_std,label = 'Mean+SD') 
plt.plot(time_series, cos_series_rms, label = 'RMS') 
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)

plt.subplot(4, 2, 4)
plt.hist(cos_series, bins)
plt.ylabel('Frequency of oc.')
plt.grid(True)

# rand
plt.subplot(4, 2, 5)
plt.plot(time_series, rand_series,
         time_series, rand_series_m,
         time_series, rand_series_m - rand_series_std,
         time_series, rand_series_m + rand_series_std,
         time_series, rand_series_rms)
plt.ylabel('Amplitude')
plt.grid(True)

plt.subplot(4, 2, 6)
plt.hist(rand_series, bins)
plt.ylabel('Frequency of oc.')
plt.grid(True)

# superposed
plt.subplot(4, 2, 7)
plt.plot(time_series, superposed_series,
         time_series, superposed_series_m,
         time_series, superposed_series_m -  superposed_series_std,
         time_series, superposed_series_m +  superposed_series_std,
         time_series, superposed_series_rms)
plt.ylabel('Amplitude')
plt.xlabel('Time [s]')
plt.grid(True)

plt.subplot(4, 2, 8)
plt.hist(superposed_series, bins)
plt.ylabel('Frequency of oc.')
plt.xlabel('Amplitude')
plt.grid(True)

#===================== BLOCK 17 ======================
iqr = scipy.stats.iqr(superposed_series)
q75, q25 = np.percentile(superposed_series, [75 ,25])
print('Interquartile range = ',iqr, 'Interquantile range computed = ', q75-q25)

#===================== BLOCK 18 ======================
plt.rcParams["figure.figsize"] = (8,8)
plt.figure(num=6)
plt.boxplot(superposed_series)
plt.show()

#===================== BLOCK 19 ======================
# const
[const_pdf_x, const_pdf_y] = c_utils.get_pdf(const_series,'Normal')

# cos
[cos_pdf_x, cos_pdf_y] = c_utils.get_pdf(cos_series)

# rand
[rand_pdf_x, rand_pdf_y] = c_utils.get_pdf(rand_series)

# superposed
[superposed_pdf_x, superposed_pdf_y] = c_utils.get_pdf(superposed_series)

#===================== BLOCK 20 ======================
# sampling frequency the same in this case for all time series
sampling_freq = 1/delta_time

# const
[const_freq_half, const_series_fft] = c_utils.get_fft(const_series, sampling_freq)

# cos
[cos_freq_half, cos_series_fft] = c_utils.get_fft(cos_series, sampling_freq)

# rand
[rand_freq_half, rand_series_fft] = c_utils.get_fft(rand_series, sampling_freq)

# superposed
[superposed_freq_half, superposed_series_fft] = c_utils.get_fft(superposed_series, sampling_freq)

#===================== BLOCK 21 ======================
# pdf, cdf and frequency domain
plt.rcParams["figure.figsize"] = (15,4)

# const
plt.figure(num=7)

plt.subplot(1,3,1)
plt.plot(const_pdf_x, const_pdf_y)
plt.xlabel(' ')
plt.ylabel('PDF(Amplitude)')
plt.title('PDF')
plt.grid(True)

plt.subplot(1,3,2)
plt.plot(const_pdf_x, np.cumsum(const_pdf_y))
plt.ylabel('CDF(Amplitude)')
plt.title('Empirical CDF')
plt.grid(True)

plt.subplot(1,3,3)
plt.plot(const_freq_half, const_series_fft)
plt.xlim([1, 25])
plt.ylabel('|Amplitude|')
plt.title('Frequency domain using FFT')
plt.grid(True)
plt.show()

# cos
plt.figure(num=8)

plt.subplot(1,3,1)
plt.plot(cos_pdf_x, cos_pdf_y)
plt.xlabel(' ')
plt.ylabel('PDF(Amplitude)')
plt.grid(True)

cos_ecdf = c_utils.get_ecdf(cos_pdf_x, cos_pdf_y)

plt.subplot(1,3,2)
plt.plot(cos_pdf_x, cos_ecdf)
plt.ylabel('CDF(Amplitude)')
plt.grid(True)

plt.subplot(1,3,3)
plt.plot(cos_freq_half, cos_series_fft)
plt.xlim([1, 25]);
plt.ylabel('|Amplitude|')
plt.grid(True)
plt.show()

# rand
plt.figure(num=9)
plt.subplot(1,3,1)
plt.plot(rand_pdf_x, rand_pdf_y)
plt.xlabel(' ')
plt.ylabel('PDF(Amplitude)')
plt.grid(True)

rand_ecdf = c_utils.get_ecdf(rand_pdf_x, rand_pdf_y)

plt.subplot(1,3,2)
plt.plot(rand_pdf_x, rand_ecdf)
plt.ylabel('CDF(Amplitude)')
plt.grid(True)

plt.subplot(1,3,3)
plt.plot(rand_freq_half, rand_series_fft)
plt.xlim([1, 25]);
plt.ylabel('|Amplitude|')
plt.grid(True)

plt.show()

# superposed
plt.figure(num=10)
plt.subplot(1,3,1)
plt.plot(superposed_pdf_x, superposed_pdf_y)
plt.xlabel(' ')
plt.ylabel('PDF(Amplitude)')
plt.xlabel('Amplitude')
plt.grid(True)

superposed_ecdf = c_utils.get_ecdf(superposed_pdf_x, superposed_pdf_y)

plt.subplot(1,3,2)
plt.plot(superposed_pdf_x, superposed_ecdf)
plt.ylabel('CDF(Amplitude)')
plt.xlabel('Amplitude')
plt.grid(True)

plt.subplot(1,3,3)
plt.plot(superposed_freq_half, superposed_series_fft)
plt.ylim([0, 0.4])
plt.xlim([1, 25])
plt.xlabel('Frequency [Hz]')
plt.ylabel('|Amplitude|')
plt.grid(True)

plt.show()

#===================== BLOCK 22 ======================
# here give the value for given_series
# as you have 4 series at hand already generated, you could
# choose one of const_series, cos_series, random_series, superposed_series
given_series = rand_series

#===================== BLOCK 23 ======================
block_size = np.round(len(given_series)/20) #  /20 -> 5% parent size for around
# 0.2% from parent distribution to be in tails 
[bm_index, bm_extreme_values] = c_utils.get_bm(given_series, block_size)
[bm_pdf_x, bm_pdf_y] = c_utils.get_pdf(bm_extreme_values)

#===================== BLOCK 24 ======================
# plotting the initial time series and selected signal series - as a line plot
plt.figure(num=11)
plt.plot(time_series, given_series, label  ='signal')

# plotting the extracted bm - as a scatter plot with round red markers
plt.scatter(time_series[bm_index], given_series[bm_index], marker = 'o', color = 'r', label='BM')
plt.ylabel('Amplitude')
plt.title('Block Maxima')
plt.xlabel('Time [s]')

# add a verticle yellow dashed line to mark the separation between blocks used for extraction
for idx in np.arange(len(bm_index)):
    plt.axvline(x=time_series[np.int(block_size * idx)], color='y', linestyle='--')
plt.axvline(x=time_series[-1], color='y', linestyle='--')

plt.legend()
plt.grid(True)

#===================== BLOCK 25 ======================
# importing additional necessary modules
from scipy.stats import genextreme as gev

# getting the fitting parameters shape, location and scale for the bm_extreme_values based upon a certain GEV fitting
bm_shape, bm_loc, bm_scale = gev.fit(bm_extreme_values)

bm_pdf_x2 = np.linspace(np.min(bm_extreme_values), np.max(bm_extreme_values), 100)
bm_pdf_y2 = gev.pdf(bm_pdf_x2, bm_shape, bm_loc, bm_scale)

plt.figure(num=12)
plt.rcParams["figure.figsize"] = (10,6)

# PDF calculated using the get_pdf from custom_function_utilities
plt.plot(bm_pdf_x, bm_pdf_y, label = 'PDF of BM')
# PDF generated as a fitted curve using generalized extreme distribution
plt.plot(bm_pdf_x2, bm_pdf_y2, label = 'PDF from the fitted GEV')

plt.xlabel('BM values')
plt.ylabel('PDF(Amplitude)')
plt.legend()
plt.show()

#===================== BLOCK 26 ======================
series_m = np.mean(given_series)
series_std = np.std(given_series)

threshold_param = 2.5
threshold_value = series_m + threshold_param * series_std # for around 0.25% from parent 

# distribution to be in tails
# here end_time means values extracted after the whole given_series is available
[pot_endtime_index, pot_endtime_extreme_values] = c_utils.get_pot(given_series, threshold_value)
[pot_endtime_pdf_x, pot_endtime_pdf_y] = c_utils.get_pdf(pot_endtime_extreme_values)

print("POT: Threshold value: ", threshold_value)

#===================== BLOCK 27 ======================
plt.figure(num=13)
# plotting the initial time series and selected signal series - as a line plot
# for this case the whole series is available
# it represents a signal being made available at the end of a mearsuremen or simulation
plt.plot(time_series, given_series, label = 'signal')

# plotting the extracted pot - as a scatter plot with round red markers
plt.scatter(time_series[pot_endtime_index], given_series[pot_endtime_index], marker ='s', color = 'r', label = 'POT')
plt.ylabel('Amplitude')
plt.title('Peak Over Threshold')

# add a horizontal yellow dashed line to mark the the two trehsholds (upper and lower) used for extraction
plt.axhline(y=threshold_value, color='y', linestyle='--')
plt.axhline(y=-threshold_value, color='y', linestyle='--')

plt.xlabel('Time [s]')
plt.legend()
plt.grid(True)

#===================== BLOCK 28 ======================
[res_m, res_rms, res_std, res_med, res_skew,
 res_thres, pot_runtime_index, pot_runtime_extreme_values] = \
                                c_utils.get_pot_runtime(given_series, threshold_param)
[pot_runtime_pdf_x, pot_runtime_pdf_y] = c_utils.get_pdf(pot_runtime_extreme_values)

#===================== BLOCK 29 ======================
plt.figure(num=14)
plt.rcParams["figure.figsize"] = (10,6)
plt.plot(time_series, res_m,
         time_series, res_std,
         time_series, res_rms,
         time_series, res_med,
         time_series, res_thres)
plt.legend(['mean, et = '+ str(np.mean(given_series)),
            'std, et = ' + str(np.std(given_series)),
            'rms, et = ' + str(np.sqrt(np.mean(np.square(given_series)))),
            'median, et = '+ str(np.median(given_series)),
            'pot thres, et = '+ str(threshold_value)])
plt.ylabel('Amplitude')
plt.title('Run-time values')
plt.xlabel('Time [s]')
plt.grid(True)

#===================== BLOCK 30 ======================
plt.figure(num=15)
plt.rcParams["figure.figsize"] = (15,4)
# plotting the initial time series and selected signal series - as a line plot
# for this case the whole series is not available, but it is being made available one time step at a time
# it represents a signal being made available in run-time
plt.plot(time_series, given_series, label = 'signal')

# plotting the extracted pot - as a scatter plot with round red markers
plt.scatter(time_series[pot_runtime_index], given_series[pot_runtime_index], marker ='s', color = 'r', label = 'POT')
plt.ylabel('Amplitude')
plt.title('Peak Over Threshold run-time evaluation')

# add a horizontal yellow dashed line to mark the the two trehsholds (upper and lower) used for extraction
# note that these are not totally straight lines from the beginning until the end, but vary slightly in time
plt.plot(time_series, res_m + threshold_param * res_std, label = 'signal', color='y', linestyle='--')
plt.plot(time_series, -res_m - threshold_param * res_std, label = 'signal', color='y', linestyle='--')

plt.xlabel('Time [s]')
plt.legend()
plt.grid(True)

#===================== BLOCK 31 ======================
# importing additional necessary modules
from scipy.stats import genpareto as gp

# getting the fitting parameters shape, location and scale for the bm_extreme_values based upon a certain GEV fitting
pot_shape, pot_loc, pot_scale = gp.fit(pot_endtime_extreme_values, 0 , loc = threshold_value , scale = 1)

pot_endtime_pdf_x2 = np.linspace(0.99 * np.min(pot_endtime_extreme_values), np.max(pot_endtime_extreme_values), 100)
pot_endtime_pdf_y2 = gp.pdf(pot_endtime_pdf_x2, pot_shape, pot_loc, pot_scale)

plt.figure(num=16)
plt.rcParams["figure.figsize"] = (10,6)

# PDF calculated using the get_pdf from custom_function_utilities
# for endtime and runtime
plt.plot(pot_endtime_pdf_x, pot_endtime_pdf_y, label = 'PDF of endtime POT')
plt.plot(pot_runtime_pdf_x, pot_runtime_pdf_y, label = 'PDF of runtime POT')

# PDF generated as a fitted curve using generalized extreme distribution
plt.plot(pot_endtime_pdf_x2, pot_endtime_pdf_y2, label = 'PDF of endtime from the fitted GEV')

plt.xlabel('POT values')
plt.ylabel('PDF(Amplitude)')
plt.legend()
plt.show()

#===================== BLOCK 32 ======================
#file_name = 'given_data1.dat'  # has 5350 values for each column
#time_series = numpy.loadtxt(file_name, skiprows=0, usecols = (0,)) # in [s]
#bending_moment_series = numpy.loadtxt(file_name, skiprows=0, usecols = (1,)) # in [kNm]

#file_name = 'given_data2.dat' # has 53491 values for each column
#time_series = (numpy.loadtxt(file_name, skiprows=0, usecols = (0,))-1000)/100  # shift 1000 [cs] then divide by 100 to get [s]
#bending_moment_series = numpy.loadtxt(file_name, skiprows=0, usecols = (1,))/1000 # to get [kNm] from [Nm] divide by 1000
