#
# LSST Data Management System
# Copyright 2008-2016 LSST Corporation.
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
from __future__ import absolute_import, division, print_function

from builtins import range

from . import diffimLib
from lsst.log import Log
import numpy as np
sigma2fwhm = 2. * np.sqrt(2. * np.log(2.))


def makeKernelBasisList(config, targetFwhmPix=None, referenceFwhmPix=None,
                        basisDegGauss=None, metadata=None):
    """Generate the appropriate Kernel basis based on the Config"""
    if config.kernelBasisSet == "alard-lupton":
        return generateAlardLuptonBasisList(config, targetFwhmPix=targetFwhmPix,
                                            referenceFwhmPix=referenceFwhmPix,
                                            basisDegGauss=basisDegGauss,
                                            metadata=metadata)
    elif config.kernelBasisSet == "delta-function":
        kernelSize = config.kernelSize
        return diffimLib.makeDeltaFunctionBasisList(kernelSize, kernelSize)
    else:
        raise ValueError("Cannot generate %s basis set" % (config.kernelBasisSet))


def generateAlardLuptonBasisList(config, targetFwhmPix=None, referenceFwhmPix=None,
                                 basisDegGauss=None, metadata=None):
    """Generate an Alard-Lupton kernel basis based upon the Config and
    the input FWHM of the science and template images"""

    if config.kernelBasisSet != "alard-lupton":
        raise RuntimeError("Cannot generate %s basis within generateAlardLuptonBasisList" %
                           config.kernelBasisSet)

    kernelSize = config.kernelSize
    fwhmScaling = config.kernelSizeFwhmScaling
    basisNGauss = config.alardNGauss
    basisSigmaGauss = config.alardSigGauss
    basisGaussBeta = config.alardGaussBeta
    basisMinSigma = config.alardMinSig
    if basisDegGauss is None:
        basisDegGauss = config.alardDegGauss

    if len(basisDegGauss) != basisNGauss:
        raise ValueError("len(basisDegGauss) != basisNGauss : %d vs %d" % (len(basisDegGauss), basisNGauss))
    if len(basisSigmaGauss) != basisNGauss:
        raise ValueError("len(basisSigmaGauss) != basisNGauss : %d vs %d" %
                         (len(basisSigmaGauss), basisNGauss))
    if (kernelSize % 2) != 1:
        raise ValueError("Only odd-sized Alard-Lupton bases allowed")

    if (targetFwhmPix is None) or (referenceFwhmPix is None) or (not config.scaleByFwhm):
        if metadata is not None:
            metadata.add("ALBasisNGauss", basisNGauss)
            metadata.add("ALBasisDegGauss", basisDegGauss)
            metadata.add("ALBasisSigGauss", basisSigmaGauss)
            metadata.add("ALKernelSize", kernelSize)

        return diffimLib.makeAlardLuptonBasisList(kernelSize//2, basisNGauss, basisSigmaGauss, basisDegGauss)

    targetSigma = targetFwhmPix / sigma2fwhm
    referenceSigma = referenceFwhmPix / sigma2fwhm
    logger = Log.getLogger("lsst.ip.diffim.generateAlardLuptonBasisList")
    logger.debug("Generating matching bases for sigma %.2f pix -> %.2f pix", targetSigma, referenceSigma)

    # Modify the size of Alard Lupton kernels based upon the images FWHM
    #
    # Note the operation is : template.x.kernel = science
    #
    # Assuming the template and science image Psfs are Gaussians with
    # the Fwhm above, Fwhm_T **2 + Fwhm_K **2 = Fwhm_S **2
    #
    if targetSigma == referenceSigma:
        # Leave defaults as-is
        pass
    elif referenceSigma > targetSigma:
        # Normal convolution

        # First Gaussian has the sigma that comes from the convolution
        # of two Gaussians : Sig_S**2 = Sig_T**2 + Sig_K**2
        #
        # If it's larger than basisMinSigma * basisGaussBeta, make it the
        # second kernel.  Else make it the smallest kernel.  Unless
        # only 1 kernel is asked for.
        kernelSigma = np.sqrt(referenceSigma**2 - targetSigma**2)
        if kernelSigma < basisMinSigma:
            kernelSigma = basisMinSigma

        basisSigmaGauss = []
        if basisNGauss == 1:
            basisSigmaGauss.append(kernelSigma)
            nAppended = 1
        else:
            if (kernelSigma/basisGaussBeta) > basisMinSigma:
                basisSigmaGauss.append(kernelSigma/basisGaussBeta)
                basisSigmaGauss.append(kernelSigma)
                nAppended = 2
            else:
                basisSigmaGauss.append(kernelSigma)
                nAppended = 1

        # Any other Gaussians above basisNGauss=1 come from a scaling
        # relationship: Sig_i+1 / Sig_i = basisGaussBeta
        for i in range(nAppended, basisNGauss):
            basisSigmaGauss.append(basisSigmaGauss[-1]*basisGaussBeta)

        kernelSize = int(fwhmScaling * basisSigmaGauss[-1])
        kernelSize += 0 if kernelSize%2 else 1 # Make sure it's odd
        kernelSize = min(config.kernelSizeMax, max(kernelSize, config.kernelSizeMin))

    else:
        # Deconvolution; Define the progression of Gaussians using a
        # method to derive a deconvolution sum-of-Gaussians from it's
        # convolution counterpart.  Only use 3 since the algorithm
        # assumes 3 components.
        #
        # http://iopscience.iop.org/0266-5611/26/8/085002  Equation 40

        # Use specializations for deconvolution
        basisNGauss = config.alardNGaussDeconv
        basisMinSigma = config.alardMinSigDeconv

        kernelSigma = np.sqrt(targetSigma**2 - referenceSigma**2)
        if kernelSigma < basisMinSigma:
            kernelSigma = basisMinSigma

        basisSigmaGauss = []
        if (kernelSigma/basisGaussBeta) > basisMinSigma:
            basisSigmaGauss.append(kernelSigma/basisGaussBeta)
            basisSigmaGauss.append(kernelSigma)
            nAppended = 2
        else:
            basisSigmaGauss.append(kernelSigma)
            nAppended = 1

        for i in range(nAppended, basisNGauss):
            basisSigmaGauss.append(basisSigmaGauss[-1]*basisGaussBeta)

        kernelSize = int(fwhmScaling * basisSigmaGauss[-1])
        kernelSize += 0 if kernelSize%2 else 1 # Make sure it's odd
        kernelSize = min(config.kernelSizeMax, max(kernelSize, config.kernelSizeMin))

        # Now build a deconvolution set from these sigmas
        sig0 = basisSigmaGauss[0]
        sig1 = basisSigmaGauss[1]
        sig2 = basisSigmaGauss[2]
        basisSigmaGauss = []
        for n in range(1, 3):
            for j in range(n):
                sigma2jn = (n - j)*sig1**2
                sigma2jn += j * sig2**2
                sigma2jn -= (n + 1)*sig0**2
                sigmajn = np.sqrt(sigma2jn)
                basisSigmaGauss.append(sigmajn)

        basisSigmaGauss.sort()
        basisNGauss = len(basisSigmaGauss)
        basisDegGauss = [config.alardDegGaussDeconv for x in basisSigmaGauss]

    if metadata is not None:
        metadata.add("ALBasisNGauss", basisNGauss)
        metadata.add("ALBasisDegGauss", basisDegGauss)
        metadata.add("ALBasisSigGauss", basisSigmaGauss)
        metadata.add("ALKernelSize", kernelSize)

    return diffimLib.makeAlardLuptonBasisList(kernelSize//2, basisNGauss, basisSigmaGauss, basisDegGauss)
