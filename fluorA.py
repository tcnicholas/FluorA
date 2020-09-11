""" 11.09.20 - thomas.c.nicholas

e-mail: thomas.nicholas@chem.ox.ac.uk
twitter: @thomascnicholas

Short Python module for reading in fluoresence data (.TXT format) within a specified
directory, applying basic data correction (subtracting buffer intensity from each
frame), and plotting.  Option to extract intensities at a particular wavelength also
supported.

"""

import glob, os
import warnings
import re
import csv

import numpy as np
import matplotlib.pyplot as plt

####################################################################################
####################################################################################

def readNstrip(filePath):
    """ Reads files and returns wavelengths and intensity counts.
    Strips out any 'Not available' intensities. """
    
    # Read in buffer.txt file as numpy array.
    lambdas, I = np.genfromtxt(filePath, delimiter=';', usecols=(0,1), unpack=True, missing_values=('Not Available'), skip_header=8)

    # Strip out 'Not Available' intensities.
    nanIndex = np.argwhere(np.isnan(I))
    lambdas = np.delete(lambdas, nanIndex)
    I = np.delete(I, nanIndex)
    
    return lambdas, I


def get_timeArray(timeInterval,Is):
    """ Returns time array in sete time intervals.
    
    Parameters:
    
    timeInterval = seconds between experiments.
    Is = frames to plot. """
    
    timeInterval /= 60

    intervals = len(Is)
    timeArray = np.arange(start=0,stop=intervals*timeInterval,step=timeInterval)
    return timeArray


def iterable(obj):
    """ Returns True if object is iterable. """
    
    try:
        iter(obj)
    except Exception:
        return False
    else:
        return True

####################################################################################
####################################################################################

class fluorAnalysis:
    """ Normalises each frame in specefied directory using 'buffer.txt', retrieving the
    wavelengths and intensity counts. """
    
    def __init__(self, directoryPath):
        self.name = directoryPath.split('/')[-1]
        self.directoryPath = directoryPath
        
        self.frameNames = None
        self.lambdas = None
        self.bufferI = None
        self.frameIs = None
        
        self.extractWavelength = None
        self.extractIs = None
        
    
    def get_buffer(self):
        """ Retrieves buffer.txt (wavelengths, intensities) as tuple.
        Updates the object buffer attributes. """
        
        bufferLambda, bufferI = readNstrip(self.directoryPath+'/buffer.txt')
        self.lambdas = bufferLambda
        self.bufferI = bufferI
        return list(zip(bufferLambda, bufferI))
    
    
    def get_frames(self,extractWavelength=None):
        """ Returns frameIs (a list of lists of frame Is), and optional
        extracted wavelength for each frame.
        
        Parameters:
        
        extractWavelength (optional): real or int. """
        
        if not self.frameIs:
            self.__normaliseFrames()
            
        if extractWavelength:
            self.extractWavelength = extractWavelength
            findWavelength = np.argwhere(self.lambdas == extractWavelength)[0][0]
            extractedIs=[]
            for frame in self.frameIs:
                extractedIs.append(frame[findWavelength])
            self.extractIs = extractedIs
            return self.frameIs, extractedIs
        return self.frameIs
    
    
    def export_lambdaIs(self,timeInterval):
        """ Writes to file, frame--wavelength--intensity(normalised)--time. """
        
        if not self.frameIs:
            self.__normaliseFrames()
        
        numFrames=len(self.frameNames)
        numWavelengths=len(self.lambdas)
        
        wavelengths = [self.lambdas] * numFrames
        wavelengths = np.concatenate(wavelengths).ravel().tolist()
        intensities = np.concatenate(self.frameIs).ravel().tolist()
        
        times = get_timeArray(timeInterval,self.frameIs)
        
        frameNames=[]
        timesRepeat=[]
        for i,name in enumerate(self.frameNames):
            frameNames+=[name]*numWavelengths
            timesRepeat+=[times[i]]*numWavelengths

        data = zip(frameNames,wavelengths,intensities,timesRepeat)
        with open(self.name+'_lambdaIs.csv', 'w') as export:
            writer = csv.writer(export)
            for row in data:
                writer.writerow(row)
    
    
    def export_iTime(self,extractWavelength,timeInterval):
        """Writes to file, frame--time--intensity.
        
        Parameters:
        
        extractWavelength: real/int, wavelength value to extract for each frame.
        timeIntervals: time intervals between each frame.
        
        """
        
        if self.extractWavelength:
            if not self.extractWavelength == extractWavelength:
                self.get_frames(extractWavelength=extractWavelength)
        else:
            self.get_frames(extractWavelength=extractWavelength)
            
        times = get_timeArray(timeInterval,self.frameIs)
        data = zip(self.frameNames,times,self.extractIs)
        with open(self.name+f'_iTime_{extractWavelength}nm.csv', 'w') as export:
            writer = csv.writer(export)
            for row in data:
                writer.writerow(row)
            
    
    def __normaliseFrames(self):
        """ Normalises all frames in directoryPath.
        Updates self.framesI with list of lists of frame intensity counts. """
        
        if not self.lambdas:
            self.get_buffer()
        
        frameIs=[]; frameNames=[]; IDretired=0
        for infile in sorted(glob.glob(self.directoryPath+'/*.TXT')):
            if infile != 'log_29May20_114926/buffer.TXT':
                
                # Check frames have been read in correct order.
                frameNumber = re.search(r'_\d{4}_', infile)
                ID = int((frameNumber.group(0) if frameNumber else '') .strip('_'))
                ID_name = str(ID).zfill(4)
                frameNames.append(ID_name)
                if ID <= IDretired:
                    warnings.warn(f'Frame order comprimised at {ID}!')
                IDretired=ID
                
                frameLambda, frameI = readNstrip(infile)
                if frameLambda.all() != self.lambdas.all():
                    warnings.warn(f'Frame {ID} wavelengths do not match buffer file!')
                
                frameI -= self.bufferI # Normalise by subtracting buffer.
                frameIs.append(frameI)
        
        self.frameNames = frameNames
        self.frameIs = frameIs
        
    ####################################################################################
        
    def plot_Ilambda(self,plotColouring=None,figsize=(18.5, 10.5),grid=True,
                     xlim=None,ylim=None,title=None,
                     xlabel='Wavelength \u03BB (nm)',
                     ylabel='Intensity (arb. units)',
                     saveFig=False):
        
        """ Intensity versus Wavelength plot.
        
        Parameters:
        
        plotColouring: (a) default None; (b) 'rainbow'; (c) (frameNumTarget,frameNumClamp,frameNumFuel).
        figsize: tuple, default (18.5, 10.5).
        grid: bool, default True.
        xlim: tuple, default auto. (lowerbound,upperbound).
        ylim: tuple, default auto. (lowerbound,upperbound).
        title: string, default None.
        xlabel: string, default 'Wavelength \u03BB (nm)'.
        ylabel: string, default 'Intensity (arb. units)'.
        saveFig: string, default empty.
        
        """
        
        if not self.frameIs:
            self.__normaliseFrames()
            
        if plotColouring == 'rainbow':
            evenSpace = np.linspace(0, 1, len(self.frameIs))
            rainbow = [plt.cm.rainbow_r(x) for x in evenSpace]
            
        plt.figure(1,figsize=figsize)
        for i,Is in enumerate(self.frameIs):
            if plotColouring == 'rainbow':
                plt.plot(self.lambdas,Is,color=rainbow[i])
            elif isinstance(plotColouring, tuple) or isinstance(plotColouring, list):
                if len(plotColouring) != 3:
                    warnings.warn('plotColouring tuple incorrect length.  Expected 3 values (target,clamp,fuel)!')
                    plt.plot(self.lambdas,Is)
                if not bool(np.product(plotColouring)) or (np.product(plotColouring) < 0):
                    warnings.warn('plotColouring tuple cannot start on frame 0 or negatives!')
                    plt.plot(self.lambdas,Is)
                else:
                    colours = (['k']*(plotColouring[0]-1)) + (['r']*(plotColouring[1]-plotColouring[0])) + (['b']*(plotColouring[2]-plotColouring[1])) + (['g']*(1+len(Is)-plotColouring[2]))
                    plt.plot(self.lambdas,Is,color=colours[i])
            else:
                plt.plot(data.lambdas,Is)
        
        
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        
        if title:
            plt.title(title)
        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        if grid:
            plt.grid('on')
        plt.show(block=True)
            
        if saveFig:
            plt.savefig(saveFig,dpi=300,transparent=True)
    
    
    def plot_iTime(self,extractWavelength,timeIntervals,
                        plotColouring=None,grid=True,
                        figsize=(18.5, 10.5),marker='x',
                        markersize=50,
                        xlim=None,ylim=None,
                        xlabel='Time (min)',
                        ylabel='Intensity (arb. units)',
                        saveFig=False):
        
        """ Extracted wavelength intensity against time (scatter).
        
        Parameters:
        
        extractWavelength: real/int, wavelength value to extract for each frame.
        timeIntervals: time intervals between each frame.
        plotColouring: (a) default None; (b) 'rainbow'; (c) (frameNumTarget,frameNumClamp,frameNumFuel).
        figsize: tuple, default (18.5, 10.5).
        grid: bool, default True.
        marker: string, default 'x'.
        markersize: real/int, default 50.
        xlim: tuple, default auto. (lowerbound,upperbound)
        ylim: tuple, default auto. (lowerbound,upperbound)
        title: string, default 'Intensity at {extractWavelength} nm'.
        xlabel: string, default 'Time (min)'.
        ylabel: string, default 'Intensity (arb. units)'
        saveFig: string, default empty.
        
        """
        
        if self.extractWavelength:
            if not self.extractWavelength == extractWavelength:
                self.get_frames(extractWavelength=extractWavelength)
        else:
            self.get_frames(extractWavelength=extractWavelength)
        
        plt.figure(2,figsize=(18.5, 10.5))
        times = get_timeArray(timeIntervals,self.extractIs)
        
        if plotColouring == 'rainbow':
            plt.scatter(times,self.extractIs,marker=marker,s=markersize,color=rainbow)
        elif isinstance(plotColouring, tuple) or isinstance(plotColouring, list):
            if len(plotColouring) != 3:
                warnings.warn('plotColouring tuple incorrect length.  Expected 3 values (target,clamp,fuel)!')
                plt.scatter(times,self.extractIs,marker=marker,s=markersize,color='k')
            if not bool(np.product(plotColouring)) or (np.product(plotColouring) < 0):
                warnings.warn('plotColouring tuple cannot start on frame 0 or negatives!')
                plt.scatter(times,self.extractIs,marker=marker,s=markersize,color='k')
            else:
                colours = (['k']*(plotColouring[0]-1)) + (['r']*(plotColouring[1]-plotColouring[0])) + (['b']*(plotColouring[2]-plotColouring[1])) + (['g']*(1+len(self.frameIs)-plotColouring[2]))
                for i in range(len(times)):
                    plt.scatter(times[i],self.extractIs[i],marker=marker,s=markersize,color=colours[i])
        else:
            plt.scatter(times,self.extractIs,marker=marker,s=markersize,color='k')
        
        plt.title(f'Intensity at {extractWavelength} nm')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        
        plt.grid('on')
        plt.show(block=True)
        
        if saveFig:
            plt.savefig(saveFig,dpi=300,transparent=True)
            
####################################################################################
####################################################################################
