#!/usr/bin/env python3
"""
09.03.21
@tcnicholas
fluorA : fluoresence experiment analysis.

Command line prompts for running analysis during an experiment.
"""

from fluorA import *

import re
import sys


def main():
    script = sys.argv[0]
    action = sys.argv[1]
    kwargs = sys.argv[2:]
    assert action in ["--watchExpt","--countsTime"], \
        f"Action is not supported: '{action}'"
    execute(action,kwargs)


def parseKWARG(kwargs,prop):
    """ Identifies kwargs from system argvs. """
    for x in kwargs:
        if re.match(fr"{prop}_",x):
            return re.sub(fr"{prop}_","",x)


def kwargDict(kwargs,pkwargs):
    """ Generates kwargs dictionary from potential list. """
    return {pkwarg:parseKWARG(kwargs,pkwarg) for pkwarg in pkwargs
        if parseKWARG(kwargs,pkwarg)}


def watch(kwargDict):
    """ Setup the routine for watching an experiment occur live. """
    # The main fluorA data processing class.
    expt = processData(kwargDict["directoryPath"])
    del kwargDict["directoryPath"]
    if "subtractBuffer" in kwargDict:
        expt.addBuffer(kwargDict["bufferPath"])
        del kwargDict["bufferPath"]
    # Run the experiment.
    expt.watchExpt(**kwargDict)


def cTime(kwargDict):
    """ Write counts-time csv file. """
    # The main fluorA data processing class.
    expt = processData(kwargDict["directoryPath"])
    del kwargDict["directoryPath"]
    if "bufferPath" in kwargDict:
        expt.addBuffer(kwargDict["bufferPath"])
        del kwargDict["bufferPath"]
    # Write the file.
    expt.countsTime(**kwargDict)


def execute(action,kwargs):
    """ Functions that can be executed from the command line. """
    if action == "--watchExpt":
        pkwargs = ["directoryPath", "updateTime", "wavelength", "timeInterval",
            "subtractDark", "subtractBuffer","bufferPath"]
        keys = kwargDict(kwargs,pkwargs)
        assert "directoryPath" in keys, \
        "Directory path must be given as argument! (directoryPath_<myDirectoryPath>)"
        if "subtractBuffer" in keys:
            assert "bufferPath" in keys, \
            "Cannot subtract buffer if buffer path not given! (bufferPath_<myBufferPath>)"
        watch(keys)
    elif action == "--countsTime":
        pkwargs = ["directoryPath","plot","write","wavelength",
            "timeInterval","fileName","bufferPath"]
        keys = kwargDict(kwargs,pkwargs)
        assert "directoryPath" in keys, \
        "Directory path must be given as argument! (directoryPath_<myDirectoryPath>)"
        cTime(keys)

if __name__ == '__main__':
   main()
