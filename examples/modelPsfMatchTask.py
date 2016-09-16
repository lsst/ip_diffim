#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2014 LSST Corporation.
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
from __future__ import print_function
import os
import sys

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.daf.base as dafBase
import lsst.meas.algorithms as measAlg
from lsst.ip.diffim import ModelPsfMatchTask


class MyModelPsfMatchTask(ModelPsfMatchTask):
    """An override for ModelPsfMatchTask"""

    def __init__(self, *args, **kwargs):
        ModelPsfMatchTask.__init__(self, *args, **kwargs)

    def run(self, templateExp, scienceExp):
        return ModelPsfMatchTask.run(self, scienceExp, templateExp.getPsf())


def createImageAndKernel(sigma, psfSize, image):
    function = afwMath.GaussianFunction2D(sigma, sigma)
    kernel = afwMath.AnalyticKernel(psfSize, psfSize, function)
    psf = measAlg.KernelPsf(kernel)

    cim = afwImage.ImageF(image.getDimensions())
    afwMath.convolve(cim, image, kernel, True)

    # Trim off the border pixels
    bbox = kernel.shrinkBBox(cim.getBBox(afwImage.LOCAL))
    cim = afwImage.ImageF(cim, bbox, afwImage.LOCAL)
    cim.setXY0(0, 0)

    return cim, psf


def generateFakeData():
    cellSize = 128
    nCell = 3
    psfSize = 21
    counts = 1e4
    border = psfSize // 2
    totalSize = nCell * cellSize + 2 * border
    templateIm = afwImage.ImageF(afwGeom.Extent2I(totalSize, totalSize))
    scienceIm = afwImage.ImageF(afwGeom.Extent2I(totalSize, totalSize))
    for x in range(nCell):
        for y in range(nCell):
            templateIm.set(x * cellSize + cellSize // 2 + border - 1,
                           y * cellSize + cellSize // 2 + border - 1,
                           counts)
            scienceIm.set(x * cellSize + cellSize // 2 + border - 1,
                          y * cellSize + cellSize // 2 + border - 1,
                          counts)

    templateImage, templatePsf = createImageAndKernel(3.0, psfSize, templateIm)
    scienceImage, sciencePsf = createImageAndKernel(2.0, psfSize, scienceIm)

    templateExp = afwImage.ExposureF(afwImage.MaskedImageF(templateImage))
    scienceExp = afwImage.ExposureF(afwImage.MaskedImageF(scienceImage))

    templateExp.setPsf(templatePsf)
    scienceExp.setPsf(sciencePsf)

    # Note here that the template image contains a reference Psf, that the science image gets matched to.
    return templateExp, scienceExp


def run(args):
    #
    # Create the Config and use sum of gaussian basis
    #
    config = ModelPsfMatchTask.ConfigClass()
    config.kernel.active.scaleByFwhm = False

    # Run the requested method of the Task
    if args.template is not None and args.science is not None:
        if not os.path.isfile(args.template):
            raise Exception("Template image %s does not exist" % (args.template))
        if not os.path.isfile(args.science):
            raise Exception("Science image %s does not exist" % (args.science))

        try:
            templateExp = afwImage.ExposureF(args.template)
        except pexExcept.LsstCppException as e:
            raise Exception("Cannot read template image %s" % (args.template))
        try:
            scienceExp = afwImage.ExposureF(args.science)
        except pexExcept.LsstCppException as e:
            raise Exception("Cannot read science image %s" % (args.science))
    else:
        templateExp, scienceExp = generateFakeData()
        config.kernel.active.sizeCellX = 128
        config.kernel.active.sizeCellY = 128

    if args.debug:
        ds9.mtv(templateExp, frame=1, title="Example script: Input Template")
        ds9.mtv(scienceExp, frame=2, title="Example script: Input Science Image")

    # Create the Task
    psfMatchTask = MyModelPsfMatchTask(config=config)

    # Run the Task
    result = psfMatchTask.run(templateExp, scienceExp)

    if args.debug:
        # See if the LSST debug has incremented the frame number; if not start with frame 3
        try:
            frame = debug.lsstDebug.frame+1
        except Exception:
            frame = 3
        ds9.mtv(result.psfMatchedExposure, frame=frame, title="Example script: Matched Science Image")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Demonstrate the use of ModelPsfMatchTask")

    parser.add_argument("--debug", "-d", action="store_true", help="Load debug.py?", default=False)
    parser.add_argument("--template", "-t", help="Template Exposure to use", default=None)
    parser.add_argument("--science", "-s", help="Science Exposure to use", default=None)

    args = parser.parse_args()

    if args.debug:
        try:
            import debug
            # Since I am displaying 2 images here, set the starting frame number for the LSST debug LSST
            debug.lsstDebug.frame = 3
        except ImportError as e:
            print(e, file=sys.stderr)

    run(args)
