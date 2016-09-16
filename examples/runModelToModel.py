from __future__ import print_function
import lsst.ip.diffim as ipDiffim

import sys
import re
import os
import numpy as num

import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersist
import lsst.afw.detection as afwDet
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.pex.policy as pexPolicy
import lsst.pex.logging as pexLog
import lsst.meas.algorithms as measAlg
import lsst.afw.display.ds9 as ds9


def pcapsf_read_boost(fn):
    import lsst.meas.algorithms as measAlgorithms # needed to register pcaPsf formatter

    print('# Reading', fn, end=' ')
    loc = dafPersist.LogicalLocation(fn)
    storageList = dafPersist.StorageList()
    additionalData = dafBase.PropertySet()
    persistence = dafPersist.Persistence.getPersistence(pexPolicy.Policy())
    storageList.append(persistence.getRetrieveStorage("BoostStorage", loc))
    psfptr = persistence.unsafeRetrieve("Psf", storageList, additionalData)
    psf = afwDet.Psf.swigConvert(psfptr)
    return psf

if __name__ == '__main__':
    import lsst.ip.diffim.diffimTools as diffimTools

    pexLog.Trace_setVerbosity("lsst.ip.diffim", 5)

    calexpPath = sys.argv[1]
    boostPath = sys.argv[2]
    sigGauss = float(sys.argv[3])

    calexp = afwImage.ExposureF(calexpPath)
    psf = pcapsf_read_boost(boostPath)
    if not calexp.hasPsf():
        calexp.setPsf(psf)

    # match to this
    gaussPsf = measAlg.DoubleGaussianPsf(psf.getKernel().getWidth(),
                                         psf.getKernel().getHeight(), sigGauss)

    config = ipDiffim.ModelPsfMatchTask.ConfigClass()
    subconfig = config.kernel
    psfMatch = ipDiffim.ModelPsfMatchTask(subconfig)
    results = psfMatch.run(calexp, gaussPsf)
    cim = results.psfMatchedExposure
    sk = results.psfMatchingKernel
    kcs = results.kernelCellSet

    if 0:
        diffimTools.displaySpatialKernelQuality(kcs, sk, sb, frame=1)
    else:
        ds9.setMaskPlaneVisibility("DETECTED", False)
        ds9.mtv(calexp, frame=1)
        diffimTools.displaySpatialKernelMosaic(psf.getKernel(), calexp.getWidth(),
                                               calexp.getHeight(),
                                               frame=2)
        diffimTools.displaySpatialKernelMosaic(gaussPsf.getKernel(), calexp.getWidth(),
                                               calexp.getHeight(),
                                               frame=3)
        diffimTools.displayKernelMosaic(kcs, frame=4)
        diffimTools.displayCandidateMosaic(kcs, frame=5)
        diffimTools.displaySpatialKernelMosaic(sk, calexp.getWidth(), calexp.getHeight(), frame=6)
        ds9.mtv(cim, frame=7)

# python examples/runModelToModel.py
# ~/LSST/PT1/psfMatch/wp_trunk_2011_0420_195756/update/calexp/v856880811-fg/R22/S20.fits
# ~/LSST/PT1/psfMatch/wp_trunk_2011_0420_195756/update/psf/v856880811-fg/R22/S20.boost
# 4.0
