#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS17-18 
        Chair of Structural Analysis @ TUM - A. Michalski, R. Wuchner, M. Pentek
        
        Function for the example of statistics, run-time evaluation and maxima 
        extraction

Author: mate.pentek@tum.de 
        Examples and methodology also based upon the ongoing PhD project and 
        M. Pentek: Wind-induced response of structures: contribution to statistical 
        quantification and passive mitigation techniques. Master Thesis, Chair 
        of Structural Analysis @TUM, 2014.  
         
Description: This is a script containing necessary function definitions for examples 
        and tasks.

Note:   It has been written and tested with Python 2.7.9. Tested and works also with Python 
        3.4.3,3.4.4,3.5.1
        Module dependencies:
            python
            matplotlib
            numpy
            OwnFunctionDef.py

Created on:  30.11.2015
Last update: 30.10.2017
'''
#===============================================================================
import matplotlib.mlab as mlab
import numpy #import max,min,std,mean,abs,floor,fft
from scipy.stats import gaussian_kde
import timeit

#===============================================================================
## Start functions
def CalcForPDF_KDE(signalGiven):
    '''
    The function CalcForPDF_KDE evaluates the probability distribution function (pdf)
    of the samples by using a non-parametric estimation technique called Kernal Desnity 
    Estimation (KDE). More details can be found at 
    https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.gaussian_kde.html.
    '''

    signalMax = numpy.max(signalGiven)
    signalMin = numpy.min(signalGiven)
    kde = gaussian_kde(signalGiven)
    pdfX = numpy.linspace(signalMin, signalMax, 1000)
    pdfY = kde(pdfX)
    return pdfX,pdfY

def CalcForPDF_Normal(signalGiven):
    '''
    The function CalcForPDF_Normal estimates the normal pdf of the signal from the mean 
    and standard deviation of the samples. Recall the fact that a Normal distribution 
    can be entirely defined by two parameters, namely the mean and standard deviation. 
    More details about the function mlab.normpdf can be found at 
    https://matplotlib.org/api/mlab_api.html. 
    '''    

    signalMax = numpy.max(signalGiven)
    signalSD= numpy.std(signalGiven)
    signalMean = numpy.mean(signalGiven)
    signalMin = numpy.min(signalGiven)
    signalStep = (signalMax - signalMin)/ 1000
    if signalSD ==0.0:
        pdfX = numpy.zeros(len(signalGiven))
        pdfY = numpy.zeros(len(signalGiven))
    else:
        signalPDF = mlab.normpdf(numpy.arange(signalMin, signalMax+signalStep, signalStep), signalMean, signalSD)
        pdfX = numpy.arange(signalMin, signalMax+signalStep, signalStep)
        pdfY = signalPDF
    return pdfX,pdfY

def CalcForPDF(signalGiven, case='KDE'):
    if case == 'KDE':
        return CalcForPDF_KDE(signalGiven)
    elif case == 'Normal':
        return CalcForPDF_Normal(signalGiven)
    else:
        raise Error("PDF type not implemented, choose either KDE or Normal")
    
def CalcForFFT(signalGiven,samplingFreq):

    signalLength=len(signalGiven)

    freqHalf =  numpy.arange(0, samplingFreq/2-samplingFreq/signalLength + samplingFreq/signalLength, samplingFreq/signalLength)
    #freqHalf = freqHalf'

    # single sided fourier
    signalFFT = numpy.fft.fft(signalGiven)
    signalFFT = numpy.abs(signalFFT[0:int(numpy.floor(signalLength/2))])/numpy.floor(signalLength/2)  
    
    maxLength = len(freqHalf)
    if maxLength < len(signalFFT):
        maxLength = len(signalFFT)
    
    freqHalf = freqHalf[:maxLength-1]
    signalFFT = signalFFT[:maxLength-1]
    
    return freqHalf,signalFFT
    
def CalcForECDF(signalPDFXGiven,signalPDF):
    
    # set up data
    dx = signalPDFXGiven[1] - signalPDFXGiven[0]    
    Y = signalPDF
    
    #normalize data
    Y /= (dx*signalPDF).sum()
    
    # compute ecdf
    CY = numpy.cumsum(Y*dx)
    
    return CY
    
def ExtractPOT(signalGiven,thresholdValue):
    
    indexPOT = []
    extremePOT = []

    for i in range(len(signalGiven)): 
         if  thresholdValue < numpy.abs(signalGiven[i]): 
             indexPOT.append(i)
             extremePOT.append(numpy.abs(signalGiven[i]))
   
    return indexPOT,extremePOT
    
def ExtractBM(signalGiven,blockSize):

    blockMax = numpy.abs(signalGiven[0])
    
    indexBM = []
    extremeBM = [] 

    for i in range(1,len(signalGiven)):
         if  blockMax < numpy.abs(signalGiven[i]):
             blockMax = numpy.abs(signalGiven[i])
             blockMaxIdx = i

         if  (i+1) % blockSize == 0: 
             indexBM.append(blockMaxIdx)
             extremeBM.append(blockMax)
            
             blockMax = 0
             
    return indexBM,extremeBM
    
def RunTimeEvalOwn(signalGiven,thresholdParam):
    
    tic=timeit.default_timer()    
    
    # renaming
    values = signalGiven
    allValues = len(signalGiven)
    
    #setting up necessary variables
    meannew = values[0]
    meanold = meannew
    rmsnew = numpy.abs(values[0])
    rmsold = rmsnew
    part1 = values[0]*values[0]
    part2 = values[0]

    arraysorted = numpy.zeros(allValues) #preallocate with zeros
    arraysorted[0] = values[0] # create sorted array with new values
    
    indexPOT = []
    extremePOT = []
    
    resMean = numpy.zeros(allValues)
    resMean[0] = meannew
    resRMS = numpy.zeros(allValues)
    resRMS[0] = rmsnew

    resSTD = numpy.zeros(allValues)
    
    resMedian = numpy.zeros(allValues)
    resMedian[0] = values[0]    
    
    resSkew = numpy.zeros(allValues)

    resThres = numpy.zeros(allValues)
    resThres[0] = meannew + thresholdParam * 0
    
    
    # starting calculation loop
    for i in range(1,allValues): 
        # calculate mean
        meannew = (meanold * (i-1) + values[i])/i
        meanold = meannew
        
        # calculate rms
        rmsnew = numpy.sqrt(((i-1)*rmsold*rmsold + values[i]*values[i])/i)
        rmsold = rmsnew
        
        # calculate standard deviation, initally summing up parts
        part1 = part1 + values[i] * values[i]
        part2 = part2 + values[i]
        standarddev = numpy.sqrt((part1 - 2* meannew * part2 + meannew*meannew* i)/(i-1))
        
        
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
        threshold = (meannew + thresholdParam * standarddev)
        if  threshold < numpy.abs(values[i]): 
             indexPOT.append(i)
             extremePOT.append(numpy.abs(values[i]))
             
                     
        # append results
        resMean[i] = meannew
        resRMS[i] = rmsnew
        resSTD[i] = standarddev
        resMedian[i] = medianv
        resSkew[i] = skewness
        resThres[i] = threshold
        
    toc=timeit.default_timer()
    
    print('Elapsed time for RunTimeEvalOwn function evaluation: ', toc-tic ,'s')
    print(' ')

    return resMean,resRMS,resSTD,resMedian,resSkew,resThres,indexPOT,extremePOT
    
#===============================================================================
## End functions
