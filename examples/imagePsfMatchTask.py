#!/usr/bin/env python

# This file is part of ip_diffim.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import sys

import lsst.afw.display as afwDisplay
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.meas.algorithms as measAlg
from lsst.afw.geom import makeSkyWcs
from lsst.ip.diffim import ImagePsfMatchTask, diffimTools


class MyImagePsfMatchTask(ImagePsfMatchTask):
    """An override for ImagePsfMatchTask"""

    def __init__(self, *args, **kwargs):
        ImagePsfMatchTask.__init__(self, *args, **kwargs)

    def run(self, templateExp, scienceExp, mode):
        if mode == "matchExposures":
            return self.matchExposures(templateExp, scienceExp)
        elif mode == "subtractExposures":
            return self.subtractExposures(templateExp, scienceExp)


def generateFakeWcs(offset=0):
    metadata = dafBase.PropertySet()
    metadata.set("SIMPLE", "T")
    metadata.set("BITPIX", -32)
    metadata.set("NAXIS", 2)
    metadata.set("NAXIS1", 425)
    metadata.set("NAXIS2", 425)
    metadata.set("RADESYS", 'FK5')
    metadata.set("EQUINOX", 2000.)
    metadata.setDouble("CRVAL1", 215.604025685476)
    metadata.setDouble("CRVAL2", 53.1595451514076)
    metadata.setDouble("CRPIX1", 1109.99981456774 + offset)
    metadata.setDouble("CRPIX2", 560.018167811613 + offset)
    metadata.set("CTYPE1", 'RA---SIN')
    metadata.set("CTYPE2", 'DEC--SIN')
    metadata.setDouble("CD1_1", 5.10808596133527E-05)
    metadata.setDouble("CD1_2", 1.85579539217196E-07)
    metadata.setDouble("CD2_2", -5.10281493481982E-05)
    metadata.setDouble("CD2_1", -8.27440751733828E-07)
    return makeSkyWcs(metadata)


def generateFakeImages():
    tSigma = 1.5
    tMi, sMi, sK, kcs, confake = diffimTools.makeFakeKernelSet(tGaussianWidth=tSigma, bgValue=200)
    sSigma = 2.5
    tWcs = generateFakeWcs()
    sWcs = generateFakeWcs()
    tExp = afwImage.ExposureF(tMi, tWcs)
    sExp = afwImage.ExposureF(sMi, sWcs)
    tPsf = measAlg.DoubleGaussianPsf(21, 21, tSigma)
    sPsf = measAlg.DoubleGaussianPsf(21, 21, sSigma)
    tExp.setPsf(tPsf)
    sExp.setPsf(sPsf)
    return tExp, sExp


def run(args):
    #
    # Create the Config and use sum of gaussian basis
    #
    config = ImagePsfMatchTask.ConfigClass()
    config.kernel.name = "AL"
    config.kernel.active.fitForBackground = True
    config.kernel.active.spatialKernelOrder = 1
    config.kernel.active.spatialBgOrder = 0

    # Run the requested method of the Task
    if args.template is not None and args.science is not None:
        if not os.path.isfile(args.template):
            raise Exception("Template image %s does not exist" % (args.template))
        if not os.path.isfile(args.science):
            raise Exception("Science image %s does not exist" % (args.science))

        try:
            templateExp = afwImage.ExposureF(args.template)
        except Exception:
            raise Exception("Cannot read template image %s" % (args.template))
        try:
            scienceExp = afwImage.ExposureF(args.science)
        except Exception:
            raise Exception("Cannot read science image %s" % (args.science))
    else:
        templateExp, scienceExp = generateFakeImages()
        config.kernel.active.sizeCellX = 128
        config.kernel.active.sizeCellY = 128

    if args.debug:
        afwDisplay.Display(frame=1).mtv(templateExp, title="Example script: Input Template")
        afwDisplay.Display(frame=2).mtv(scienceExp, title="Example script: Input Science Image")

    # Create the Task
    psfMatchTask = MyImagePsfMatchTask(config=config)

    # Run the Task
    result = psfMatchTask.run(templateExp, scienceExp, args.mode)

    if args.debug:
        # See if the LSST debug has incremented the frame number; if not start with frame 3
        try:
            frame = debug.lsstDebug.frame + 1
        except Exception:
            frame = 3
        afwDisplay.Display(frame=frame).mtv(result.matchedExposure,
                                            title="Example script: Matched Template Image")
        if "subtractedExposure" in result.getDict():
            afwDisplay.Display(frame=frame + 1).mtv(result.subtractedExposure,
                                                    title="Example script: Subtracted Image")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Demonstrate the use of ImagePsfMatchTask")

    parser.add_argument("--debug", "-d", action="store_true", help="Load debug.py?", default=False)
    parser.add_argument("--template", "-t", help="Template Exposure to use", default=None)
    parser.add_argument("--science", "-s", help="Science Exposure to use", default=None)
    parser.add_argument("--mode", choices=["matchExposures", "subtractExposures"],
                        default="subtractExposures", help="Which method of ImagePsfMatchTask to invoke")

    args = parser.parse_args()

    if args.debug:
        try:
            import debug
            # Since I am displaying 2 images here, set the starting frame number for the LSST debug LSST
            debug.lsstDebug.frame = 3
        except ImportError as e:
            print(e, file=sys.stderr)

    run(args)
