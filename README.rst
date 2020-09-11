
.. image:: FluorAnalysis.jpg
           :alt: fluorA logo

FluorA
=============================

Fluorescence Analysis Toolkit.

A basic Python package for reading, applying simple corrections, and plotting fluorescence data.

Requirements
------------

* Python_ 3.6 or later
* NumPy_ (vector and matrix algebra)
* MatPlotLib_ (for plotting)

Example:
------------

Reading in datafile and plotting:

* Intensity versus Wavelength.
* Intensity of a given wavelength versus time.

>>> from fluorA import *
>>> timeIntervals = 45
>>> extractWavelength = 520

>>> expt = fluorAnalysis('directoryName') # read data file.
>>> expt.plot_Ilambda(plotColouring='rainbow', grid=False, saveFig='image.pdf') # plots I vs lambda. exports fig image.
>>> expt.plot_iTime(extractWavelength, timeIntervals) # plots I vs time (for 520 nm).

>>> expt.export_lambdaIs(timeIntervals) # exports I vs lamda data.
>>> expt.export_iTime(extractWavelength, timeIntervals) # exports I vs time (for 520 nm) data.


.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _MatPlotLib: https://matplotlib.org

