import numpy
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.pex.exceptions as pexExcept

#######
# Helper functions used by all ip_diffim packages
#######

def findIqr(myArray):
    # Naive way to compute the interquartile range.  This should be
    # implemented thoroughly in afwMath
    
    aSort = numpy.sort(myArray)

    idx25 = int(0.25 * len(aSort))
    idx50 = int(0.50 * len(aSort))
    idx75 = int(0.75 * len(aSort))
    d25   = aSort[idx25]
    d50   = aSort[idx50]
    d75   = aSort[idx75]

    median         = d50
    Iqr            = (d75 - d25)
    sigma_mean     = 0.741 * Iqr
    sigma_median   = numpy.sqrt(numpy.pi / 2.) * sigma_mean
    
    return median, Iqr, sigma_mean, sigma_median


def fitFunction(function, values, errors, cols, rows, policy):
    # Fit a function to data
    nSigmaSq = policy.get('nSigmaSq')
    stepsize = policy.get('stepsize')

    # initialize fit parameters
    nPar = function.getNParameters()
    pars = numpy.zeros(nPar)       # start with no spatial variation
    pars[0]  = numpy.mean(values)  # except for the constant term
    stepsize = stepsize * numpy.ones(nPar)
    fit      = afwMath.minimize(function,
                                pars,
                                stepsize,
                                values,
                                errors,
                                cols,
                                rows,
                                nSigmaSq)
    if not fit.isValid:
        raise RuntimeError('Spatial fit fails')

    return fit

#######
# Expansions of functionality found in lsst.afw.image.testUtils
#######

def vectorFromImage(im, dtype=float):
    vec = numpy.zeros(im.getWidth()*im.getHeight(), dtype=dtype)
    idx = 0
    for row in range(im.getHeight()):
        for col in range(im.getWidth()):
            vec[idx] = im.get(col, row)
            idx     += 1
    return vec

def imageFromVector(vec, width, height, retType=afwImage.ImageF):
    im  = retType(width, height)
    idx = 0
    for row in range(height):
        for col in range(width):
            # need to cast numpy.float64 as float
            im.set(col, row, float(vec[idx]))
            idx     += 1
    return im

