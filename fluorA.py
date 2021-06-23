#!/usr/bin/env python3
"""
21.02.21
@tcnicholas
fluorA : fluoresence experiment analysis.
"""

import re
import os
import csv
import copy
import pathlib
import warnings

import numpy as np
from collections import OrderedDict

import matplotlib.pyplot as plt
import matplotlib.animation as animation


class processData:
    def __init__(self,directoryPath,pathToBuffer="AutoSaveBuffer"):
        self.directoryPath = directoryPath
        self.buffer = None
        self._buffer(pathToBuffer,auto=True) # Initial search for buffer file.
        self.scanNumbers = set()
        self.scans = OrderedDict()
        self.getScans() # Run initial search for scans.

        # Figure objects for watching experiment live:
        self.fig = None
        self.ax = None


    def addBuffer(self,pathToBuffer="AutoSaveBuffer"):
        """ Gets buffer file. """
        self.buffer = self._buffer(pathToBuffer)


    def watchExpt(self, directoryPath=None, updateTime=10000, wavelength=520,
                    timeInterval=30, subtractDark=False, subtractBuffer=False):
        """
        Creates matplotlib graph that automatically updates until window
        is closed.
        """
        directoryPath = self.directoryPath if directoryPath is None else directoryPath
        # Create global figure attributes.
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(1, 1, 1)
        animate = animation.FuncAnimation(self.fig, self._watchExpt,
                        interval=updateTime,
                        fargs=(directoryPath, float(wavelength),
                        float(timeInterval), subtractDark, subtractBuffer))
        plt.show()


    def getScans(self,directoryPath=None):
        """ Get new scan from directory. """
        # Unless otherwise states, will search same directory for files.
        directoryPath = self.directoryPath if directoryPath is None else directoryPath
        files, _ = getFiles(directoryPath)

        # For any new files in the directory, parase and add to self.scans dictionary.
        newFiles = set(files.keys()) - self.scanNumbers
        if newFiles:
            self.scans.update({f:parseFile(directoryPath+"/"+files[f])
                                                            for f in newFiles})
            self.scanNumbers |= newFiles


    def countsTime(self, plot=False, write=False, get=False, wavelength=520,
                    timeInterval=30, fileName=""):
        """ Gets counts vs time for a given wavelength. """
        # Check at least one action is specified.
        if not (plot or write or get):
            raise ValueError("Specify either 'plot', 'write', or 'get'!")

        # Extract w:wavelength, counts, counts-dark, counts-buffer.
        # Nb: if no buffer file, countsSubBuffer will be empty dict.
        w, counts, countsSubDark, countsSubBuffer = self._getWavelengthCounts(
                                                wavelength=float(wavelength))
        times = timeArray(float(timeInterval),len(counts))

        # Check if filename given.
        fileName = fileName if fileName else f"countsTime{w}nm.csv"

        # Incase filename does not have (the correct) file format specified.
        sFileName = fileName.split(".")
        if sFileName == 1 or sFileName[-1] != ".csv":
            fileName = sFileName[0] + ".csv"

        # Define headers for csv file columns.
        header = ["Time (mins)","Intensity - raw (counts)",
                    "Intensity - subDark (counts)"]
        data = zip(times, counts, countsSubDark)

        # Add buffer subtraction value if available.
        if self.buffer:
            header.append("Intensity - subBuffer (counts)")
            data = zip(times, counts, countsSubDark, countsSubBuffer)

        if write:
            # Export to file using csv writer.
            with open(pathlib.Path(fileName),"w+") as file:
                writer = csv.writer(file,delimiter=",")
                writer.writerow(h for h in header)
                for row in data:
                    writer.writerow(row)

            # Confirm file written.
            print(f">>> Data written to file: '{pathlib.Path(fileName)}'")

        elif plot:
            pass

        print(f"Headers:\n{'; '.join(header)}")
        return tuple(zip(*data))


    def getCountsTime(self,wavelength=520,timeInterval=30,subtractBuffer=False,
                        subtractDark=False,scanSegments=[],timeSegments=[],
                        removeScans=[],removeTimes=[],asDict=False):
        """
        Gets OrderedDict of {scans:{"counts" : counts}}
        """
        # Check sorted and convert to numpy array.
        # Will be used to add a colour tag for each scan for plotting.
        if scanSegments:
            segments = np.array(sorted(scanSegments))
        elif timeSegments:
            segments = np.array(sorted(timeSegments))

        w, counts, countsSubDark, countsSubBuffer = self._getWavelengthCounts(
                                                wavelength=float(wavelength))
        times = timeArray(float(timeInterval),len(counts))
        w = re.sub(r"\(",r".",w)
        w = re.sub(r"\)",r"",w)

        if subtractBuffer:
            useCount = copy.copy(countsSubBuffer)
        elif subtractDark:
            useCount = copy.copy(countsSubDark)
        else:
            useCount = copy.copy(counts)

        # Loop through all counts and create plottable lists.
        scanNumbers = []; finalTimes = []; finalCounts = []; colourTags = []
        for i,x in enumerate(useCount):
            if (i+1) not in removeScans and times[i] not in removeTimes:
                finalCounts.append(x)
                finalTimes.append(times[i])
                scanNumbers.append(i+1)
                if scanSegments:
                    colourTags.append(segments.searchsorted(i+1,'right') - 1)
                elif timeSegments:
                    colourTags.append(segments.searchsorted(times[i],
                                                                'right') - 1)
                else:
                    colourTags.append(0)

        return w, scanNumbers, finalTimes, finalCounts, colourTags


    def getIntensityWavelength(self,subtractBuffer=False,subtractDark=False,
                                    scanSegments=[],timeSegments=[],
                                    timeInterval=30):
        """
        Gets OrderedDict of {scans:{"wavelengths" : wavelengths,
                                        "counts" : counts}}
        """
        # Check sorted and convert to numpy array.
        # Will be used to add a colour tag for each scan for plotting.
        if scanSegments:
            segments = np.array(sorted(scanSegments))
        elif timeSegments:
            times = timeArray(float(timeInterval),len(self.scans))
            segments = np.array(sorted(timeSegments))

        if subtractBuffer:
            # Get buffer first if needed, rather than repeating in loop.
            # Store as numpy array so it can be easily subtracted.
            assert self.buffer is not None, \
            f"No buffer file found. Add using 'addBuffer()'."
            bufferCounts = np.array(self.buffer["Sample"]["values"])

        # For each scan, get all wavelengths and counts, subtract either dark or
        # buffer as specified, and save to OrderedDict.
        # Note: wavelengths could be stored once (given they should all be the
        # same for all scans, but will make it easier to plot directly later).
        intensityWavelength = OrderedDict()
        for s in self.scans:
            wavelengths = self.scans[s]["Wave"]["values"]
            counts = np.array(self.scans[s]["Sample"]["values"])
            if subtractDark:
                dark = np.array(self.scans[s]["dark"]["values"])
                counts -= dark
            elif subtractBuffer:
                counts -= bufferCounts

            # Create dictionary entry.
            intensityWavelength[s] = {"wavelengths" : wavelengths,
                                        "counts" : counts}

            # If splitting up segments by scan numbers, also add colour tag.
            if scanSegments:
                seg = segments.searchsorted(s,'right') - 1
                intensityWavelength[s].update({"segment":seg})
            elif timeSegments:
                seg = segments.searchsorted(times[s-1],'right') - 1
                intensityWavelength[s].update({"segment":seg})
            else:
                intensityWavelength[s].update({"segment":0})

        return intensityWavelength


    def _buffer(self,pathToBuffer,auto=False):
        """ Gets buffer file. """
        try:
            # Search buffer directory for .txt files.
            files, otherFiles = getFiles(pathToBuffer)

            # Combine all files if present.
            if files and otherFiles:
                files = OrderedDict(list(files.items()) + list(otherFiles.items()))
            elif otherFiles:
                files = otherFiles

            if len(files) > 1:
                # If multiple .txt files found in buffer folder, asks which file to use.
                # Nb: this is case sensitive!
                keys = [str(k) for k in files]
                print(f"Found mulitiple potential buffer files:\n{', '.join(keys)}")
                bufferFile = ""
                keys += ["none"]
                while bufferFile not in keys:
                    bufferFile = input(f"Which one? (or 'none')").strip()
                if bufferFile != "none":
                    i = [k for k in files][keys.index(bufferFile)]
                    return parseFile(pathToBuffer+"/"+files[i])
            else:
                return parseFile(pathToBuffer+"/"+list(files.values())[0])
        except:
            if not auto:
                warnings.warn(f"No buffer file found in '{pathToBuffer}/' location.")


    def _getWavelengthCounts(self,wavelength):
        """ Extract counts at a given wavelength. """
        counts = []; countsSubDark = []; countsSubBuffer = []
        for s in self.scans:
            wavelengths = self.scans[s]["Wave"]["values"]
            try:
                # Try getting counts at specified wavelength.
                w = wavelengths.index(wavelength)
                usedW = str(wavelength)
            except:
                # If that wavelength is not in the list, finds closest
                # wavelength to that specified.
                close = nearest(wavelengths,wavelength)
                warnings.warn(f"{wavelength} nm not in scans. Using {close} nm instead.")

                # Store new wavelength for file output name.
                usedW = str(close).replace(".","(") + ")"
                w = wavelengths.index(close)

            # Then get all intensity count (with subtractionns).
            count = self.scans[s]["Sample"]["values"][w]
            counts.append(count)

            # Get counts - dark.
            dark = self.scans[s]["Dark"]["values"][w]
            countsSubDark.append(count-dark)

            # If buffer file available, also subtract buffer.
            if self.buffer:
                buffer = self.buffer["Sample"]["values"][w]
                countsSubBuffer.append(count-buffer)

        return usedW, counts, countsSubDark, countsSubBuffer


    def _watchExpt(self, i, directoryPath, wavelength, timeInterval,
                    subtractDark, subtractBuffer):
        """ Automatically updating figure for watching data live. """

        # First get any new files.
        self.getScans(directoryPath)

        # Must have at least one scan.
        try:
            assert len(self.scans) >= 1
        except:
            raise AssertionError("Need at least one scan before plotting.")

        # Extract w:wavelength, counts, counts-dark, counts-buffer.
        # Nb: if no buffer file, countsSubBuffer will be empty dict.
        w, counts, countsSubDark, countsSubBuffer = self._getWavelengthCounts(wavelength=wavelength)
        times = timeArray(timeInterval,len(counts))

        # Choose correct counts as specified above.
        # Nb. If both subtractDark and subtractBuffer are chosen, only
        # dark will be subtracted.
        if subtractDark:
            c = countsSubDark
        elif subtractBuffer:
            if not self.buffer:
                raise ValueError("No buffer file found. Add one first, or don't ask for it!")
            c = countsSubBuffer
        else:
            c = counts

        # Update axis.
        self.ax.clear()
        self.ax.scatter(times,c,marker="x",color="k")

        # Format y-label.
        if subtractDark:
            y = " - dark"
        elif subtractBuffer:
            y = " - buffer"
        else:
            y = ""

        # Format so decimal figures appear after a dot,
        # rather than in brackets.
        w = w.replace("(",".")
        w = w.replace(")","")

        # Set plot and axis titles.
        self.ax.set_title(f"Intensity of {w} nm wavelength vs time.")
        self.ax.set_xlabel("Time (mins)")
        self.ax.set_ylabel(f"Intensity{y} (counts)")


def sqBrackets(s):
    """ Checks if string contains [word in square brackets]. """
    return bool(re.search(r'\[\w+\]', s))


def wordInSqBrackets(s):
    """ Extracts word from square brackets. """
    result = re.search(r"\[([A-Za-z0-9_]+)\]", s)
    return result.group(1)


def parseFile(filePath):
    """ Parses single file. Returns headers, and data columns in Python dictionary. """
    filePath = pathlib.Path(filePath) # Create Path object to avoid mac/windows issues.
    fileDict = {}
    vals = []
    with open(filePath, "r") as f:
        lines = f.readlines()
        # Go through line by line.
        while lines:
            l = lines.pop(0).strip()
            if l:
                header = l.split(":") # Headers have a key:value format.
                if len(header) == 2:
                    h = header[0].strip() # Key.
                    v = header[1].strip() # Value.
                    fileDict[h] = v # Store as dictionary.
                else:
                    del header
                    data = l.split(";") # Data separated by semi-colons.
                    try:
                        # See if numbers present in line.
                        vals.append(list(np.array(data).astype(np.float)))
                    except:
                        # If not numbers, assume column headings (or units).
                        # This block assumes that the column headings are given
                        # before the column units (calls the dataCols list).
                        isUnits = all([sqBrackets(v) for v in data])
                        if isUnits:
                            # For each unit, create dictionary entry with column header.
                            for i,unit in enumerate(data):
                                u = wordInSqBrackets(unit.strip())
                                fileDict[dataCols[i]].update({"units":u})
                        else:
                            # Create empty dictionary entry for each column header.
                            dataCols = [x.strip() for x in data]
                            for c in dataCols:
                                fileDict[c] = {}

    # Then append the column values to the appropriate column header in fileDict.
    for i,v in enumerate(zip(*vals)):
        fileDict[dataCols[i]].update({"values" : v})

    return fileDict


def getFiles(directoryPath,fileFormat="txt"):
    """ Gets list of files in directory. """
    # Create Path object to avoid mac/windows issues.
    search = pathlib.Path(directoryPath)

    # Get all files in the directory with the correct fileFormat.
    files = sorted([f for f in os.listdir(search)
                    if os.path.isfile(os.path.join(search, f))
                    and f.split('.')[-1].lower()==fileFormat])

    # Then generate filesPathDict (to help order the files correctly).
    # This stores {fileNum:pathToFile} as {key:value} pairs.
    filesPathDict = {}
    otherFiles = {}
    for file in files:
        f = file.split(".")[0]
        try:
            num = int(f.split("_")[-1])
            filesPathDict[num] = file
        except:
            # If file name does not have "_\d+".
            otherFiles[f] = file

    # Return as ordered dictionary.
    return OrderedDict(filesPathDict), OrderedDict(otherFiles)


def timeArray(timeInterval,numScans):
    """ Returns time array in set time intervals.
    Args:
        timeInterval (float): seconds between experiments.
        numScans (int): number of scans.
    """
    timeInterval /= 60
    return np.arange(start=0,stop=numScans*timeInterval,step=timeInterval)


def nearest(array, value):
    """ Find nearest element value in array to specified value. """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
