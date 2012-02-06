# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import diffimLib 
import lsst.pex.config as pexConfig
import numpy as num
sigma2fwhm = 2. * num.sqrt(2. * num.log(2.))

def makeKernelBasisList(config, templateFwhmPix = None, scienceFwhmPix = None):
    if config.kernelBasisSet == "alard-lupton":
        return generateAlardLuptonBasisList(config, templateFwhmPix, scienceFwhmPix)
    elif config.kernelBasisSet == "delta-function":
        kernelSize = config.kernelSize
        return diffimLib.makeDeltaFunctionBasisList(kernelSize, kernelSize)
    else:
        raise ValueError("Cannot generate %s basis set" % (config.kernelBasisSet))

def generateAlardLuptonBasisList(config, templateFwhmPix = None, scienceFwhmPix = None, minSigma = 0.4, minRatio = 1.25):
    if config.kernelBasisSet != "alard-lupton":
        raise RuntimeError("Cannot generate %s basis within generateAlardLuptonBasisList" % (config.kernelBasisSet))

    kernelSize    = config.kernelSize
    alardNGauss   = config.alardNGauss
    alardDegGauss = config.alardDegGauss
    alardSigGauss = config.alardSigGauss

    if len(alardDegGauss) != alardNGauss:
        raise ValueError("len(alardDegGauss) != alardNGauss : %d vs %d" % (len(alardDegGauss), alardNGauss))
    if len(alardSigGauss) != alardNGauss:
        raise ValueError("len(alardSigGauss) != alardNGauss : %d vs %d" % (len(alardSigGauss), alardNGauss))
    if (kernelSize % 2) != 1:
        raise ValueError("Only odd-sized Alard-Lupton bases allowed")
        
    if not config.scaleByFwhm:
        return diffimLib.makeAlardLuptonBasisList(kernelSize//2, alardNGauss, alardSigGauss, alardDegGauss)
    
    if (templateFwhmPix == None) or (scienceFwhmPix == None):
        return diffimLib.makeAlardLuptonBasisList(kernelSize//2, alardNGauss, alardSigGauss, alardDegGauss)
        

    # Modify the size of Alard Lupton kernels based upon the images FWHM
    #
    # Note the operation is : template.x.kernel = science
    #
    # Assuming the template and science image Psfs are Gaussians with
    # the Fwhm above, Fwhm_T **2 + Fwhm_K **2 = Fwhm_S **2
    #
    if templateFwhmPix == scienceFwhmPix:
        # Leave defaults as-is
        pass
    elif float(scienceFwhmPix) / float(templateFwhmPix) > 2.0:
        # Extreme convolution; central Gaussian is at the template
        # scale, outer Gaussian is at the scale to match the two
        # cores, central Gaussian is geometric mean.
        kernelCoreSigma   = templateFwhmPix / sigma2fwhm
        kernelOuterSigma  = num.sqrt(scienceFwhmPix**2 - templateFwhmPix**2) / sigma2fwhm

        # Minimum ratio
        ratio             = num.sqrt(kernelOuterSigma / kernelCoreSigma)
        if ratio < minRatio:
            kernelOuterSigma = minRatio**2 * kernelCoreSigma
        kernelMiddleSigma = num.sqrt(kernelCoreSigma * kernelOuterSigma)

        alardNGauss   = 3
        alardSigGauss = (kernelCoreSigma, kernelMiddleSigma, kernelOuterSigma)
        alardDegGauss = alardDegGauss[:3]

    elif scienceFwhmPix > templateFwhmPix:
        # Normal convolution; put the bulk of the power in the Gaussian
        # that matches the core Fwhms.  Outer gaussian corresponds to
        # the science image's Fwhm.  Middle is geometric mean to create
        # geometric progression of Gaussian sizes
        kernelCoreFwhm    = num.sqrt(scienceFwhmPix**2 - templateFwhmPix**2)
        kernelCoreSigma   = max(minSigma, kernelCoreFwhm / sigma2fwhm)
        kernelOuterSigma  = scienceFwhmPix / sigma2fwhm

        # Minimum ratio
        ratio             = num.sqrt(kernelOuterSigma / kernelCoreSigma)
        if ratio < minRatio:
            kernelOuterSigma = minRatio**2 * kernelCoreSigma
        kernelMiddleSigma = num.sqrt(kernelCoreSigma * kernelOuterSigma)

        alardNGauss   = 3
        alardSigGauss = (kernelCoreSigma, kernelMiddleSigma, kernelOuterSigma)
        alardDegGauss = alardDegGauss[:3]

    else:
        # Deconvolution; put the smallest Gaussian at the smallest
        # allowed scale, and define the progression of Gaussians using
        # a method used to derive a deconvolution sum-of-Gaussians
        # from its convolution counterpart.
        #
        # http://iopscience.iop.org/0266-5611/26/8/085002  Equation 40
        kernelDeconvSigma = minSigma
        kernelCoreSigma   = minSigma
        kernelOuterSigma  = templateFwhmPix / sigma2fwhm

        # Minimum ratio
        ratio             = num.sqrt(kernelOuterSigma / kernelCoreSigma)
        if ratio < minRatio:
            kernelOuterSigma = minRatio**2 * kernelCoreSigma
        kernelMiddleSigma = num.sqrt(kernelCoreSigma * kernelOuterSigma)

        alardNGauss   = config.alardNGaussDeconv
        alardDegGauss = config.alardDegGaussDeconv
        alardSigGauss = []
        alardSigGauss.append(kernelDeconvSigma)
        for n in range(3):
            for j in range(n):
                sigma2jn  = (n - j) * kernelMiddleSigma**2
                sigma2jn += j * kernelOuterSigma**2
                sigma2jn -= (n + 1) * kernelCoreSigma**2
                sigmajn   = num.sqrt(sigma2jn)
                alardSigGauss.append(sigmajn)

    return diffimLib.makeAlardLuptonBasisList(kernelSize//2, alardNGauss, alardSigGauss, alardDegGauss)

