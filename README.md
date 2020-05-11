# Mars-Crater-Intensity-Analysis

README

cratercurves.py is a file that will accept a clip of a crater taken from the Mars CTX Mosaic (link: http://murray-lab.caltech.edu/CTX/) in .tif format and will plot the frequency of pixel intensity (aka. darkness) as a curve. This curve will then be broken down into multiple Gaussian curves using non linear least squares curve fitting (particularly scipy.optimize.least_squares). Then, each Gaussian curve is integrated to find the total area for different intensity ranges, essentially describing the surface area of sediment that has been actively moving/recently deposited into the crater.     

CONTACT

If you have any problems, questions, or suggestions, please email saking3@uw.edu

WEBSITE

Please visit the linked github repository for updates and downloads: https://github.com/saking3/Mars-Crater-Intensity-Analysis
