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
import lsst.pex.logging as pexLog
import numpy as np
sigma2fwhm = 2. * np.sqrt(2. * np.log(2.))

def makeKernelBasisList(config, targetFwhmPix = None, referenceFwhmPix = None, 
                        alardDegGauss = None, metadata = None):
    if config.kernelBasisSet == "alard-lupton":
        return generateAlardLuptonBasisList(config, targetFwhmPix=targetFwhmPix, 
                                            referenceFwhmPix=referenceFwhmPix, 
                                            alardDegGauss=alardDegGauss,
                                            metadata=metadata)
    elif config.kernelBasisSet == "delta-function":
        kernelSize = config.kernelSize
        return diffimLib.makeDeltaFunctionBasisList(kernelSize, kernelSize)
    else:
        raise ValueError("Cannot generate %s basis set" % (config.kernelBasisSet))

def generateAlardLuptonBasisList(config, targetFwhmPix = None, referenceFwhmPix = None, 
                                 alardDegGauss = None, metadata = None):
    if config.kernelBasisSet != "alard-lupton":
        raise RuntimeError("Cannot generate %s basis within generateAlardLuptonBasisList" % (
                config.kernelBasisSet))

    kernelSize     = config.kernelSize
    fwhmScaling    = config.kernelSizeFwhmScaling
    alardNGauss    = config.alardNGauss
    alardSigGauss  = config.alardSigGauss
    alardGaussBeta = config.alardGaussBeta
    alardMinSig    = config.alardMinSig
    if alardDegGauss == None:
        alardDegGauss = config.alardDegGauss

    if len(alardDegGauss) != alardNGauss:
        raise ValueError("len(alardDegGauss) != alardNGauss : %d vs %d" % (len(alardDegGauss), alardNGauss))
    if len(alardSigGauss) != alardNGauss:
        raise ValueError("len(alardSigGauss) != alardNGauss : %d vs %d" % (len(alardSigGauss), alardNGauss))
    if (kernelSize % 2) != 1:
        raise ValueError("Only odd-sized Alard-Lupton bases allowed")
        
    if (targetFwhmPix == None) or (referenceFwhmPix == None) or (not config.scaleByFwhm):
        if metadata is not None:
            metadata.add("ALBasisNGauss", alardNGauss)
            metadata.add("ALBasisDegGauss", alardDegGauss)
            metadata.add("ALBasisSigGauss", alardSigGauss)
            metadata.add("ALKernelSize", kernelSize)

        return diffimLib.makeAlardLuptonBasisList(kernelSize//2, alardNGauss, alardSigGauss, alardDegGauss)

    targetSigma    = targetFwhmPix / sigma2fwhm
    referenceSigma = referenceFwhmPix / sigma2fwhm
    pexLog.Trace("lsst.ip.diffim.generateAlardLuptonBasisList", 2,
                 "Generating matching bases for sigma %.2f pix -> %.2f pix" % (targetSigma, referenceSigma))

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
        # If its larger than alardMinSig * alardGaussBeta, make it the
        # second kernel.  Else make it the smallest kernel.  Unless
        # only 1 kernel is asked for.
        kernelSigma = np.sqrt(referenceSigma**2 - targetSigma**2)
        if kernelSigma < alardMinSig:
            kernelSigma = alardMinSig

        alardSigGauss = []
        if alardNGauss == 1:
            alardSigGauss.append(kernelSigma)
            loc = 1
        else:
            if (kernelSigma/alardGaussBeta) > alardMinSig:
                alardSigGauss.append(kernelSigma/alardGaussBeta)
                alardSigGauss.append(kernelSigma)
                loc = 2
            else:
                alardSigGauss.append(kernelSigma)
                loc = 1

        # Any other Gaussians above alardNGauss=1 come from a scaling
        # relationship: Sig_i+1 / Sig_i = alardGaussBeta
        for i in range(loc,alardNGauss):
            alardSigGauss.append(alardSigGauss[-1]*alardGaussBeta)

        kernelSize  = int(fwhmScaling * alardSigGauss[-1])
        kernelSize += 0 if kernelSize%2 else 1 # Make sure its odd
        kernelSize  = min(config.kernelSizeMax, max(kernelSize, config.kernelSizeMin))

    else:
        # Deconvolution; Define the progression of Gaussians using a
        # method to derive a deconvolution sum-of-Gaussians from its
        # convolution counterpart.  Only use 3 since the algorithm
        # assumes 3 components.
        #
        # http://iopscience.iop.org/0266-5611/26/8/085002  Equation 40
        alardNGauss = 3
        kernelSigma = np.sqrt(targetSigma**2 - referenceSigma**2)
        if kernelSigma < alardMinSig:
            kernelSigma = alardMinSig

        alardSigGauss = []
        if (kernelSigma/alardGaussBeta) > alardMinSig:
            alardSigGauss.append(kernelSigma/alardGaussBeta)
            alardSigGauss.append(kernelSigma)
            loc = 2
        else:
            alardSigGauss.append(kernelSigma)
            loc = 1

        for i in range(loc,alardNGauss):
            alardSigGauss.append(alardSigGauss[-1]*alardGaussBeta)

        kernelSize  = int(fwhmScaling * alardSigGauss[-1])
        kernelSize += 0 if kernelSize%2 else 1 # Make sure its odd
        kernelSize  = min(config.kernelSizeMax, max(kernelSize, config.kernelSizeMin))

        # Now build a deconvolution set from these sigmas
        sig0 = alardSigGauss[0]
        sig1 = alardSigGauss[1]
        sig2 = alardSigGauss[2]
        alardSigGauss = []
        for n in range(3):
            for j in range(n):
                sigma2jn  = (n - j) * sig1**2
                sigma2jn += j * sig2**2
                sigma2jn -= (n + 1) * sig0**2
                sigmajn   = np.sqrt(sigma2jn)
                alardSigGauss.append(sigmajn)
        alardSigGauss.sort()
        alardNGauss = len(alardSigGauss)
        alardDegGauss = [config.alardDegGaussDeconv for x in alardSigGauss]

    if metadata is not None:
        metadata.add("ALBasisNGauss", alardNGauss)
        metadata.add("ALBasisDegGauss", alardDegGauss)
        metadata.add("ALBasisSigGauss", alardSigGauss)
        metadata.add("ALKernelSize", kernelSize)

    return diffimLib.makeAlardLuptonBasisList(kernelSize//2, alardNGauss, alardSigGauss, alardDegGauss)

