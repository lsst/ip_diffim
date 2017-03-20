from __future__ import absolute_import, division, print_function
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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
import sys
import optparse

import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim
import lsst.log.utils as logUtils


def subtractSnaps(snap1, snap2, subconfig, doWarping=False):
    psfmatch = ipDiffim.SnapPsfMatchTask(subconfig)
    results = psfmatch.run(snap1, snap2, "subtractExposures", doWarping=doWarping)
    return results.subtractedImage


def main():
    defVerbosity = 5

    usage = """usage: %%prog [options] snap1 snap2 snapdiff
    
    Notes:
    - image arguments are paths to Expsosure (calexp) fits files
    - snap1 is convolved
    """
    parser = optparse.OptionParser(usage)
    parser.add_option('--s1', help='snap1')
    parser.add_option('--s2', help='snap2')
    parser.add_option('--sdiff', help='snap2 - snap1.x.kernel')
    parser.add_option('--warp', action='store_true', default=False, help='astrometrically warp snap1')
    parser.add_option('-v', '--verbosity', type=int, default=defVerbosity,
                      help='verbosity of Trace messages')

    (options, args) = parser.parse_args()
    if options.s1 == None or options.s2 == None or options.sdiff == None:
        parser.print_help()
        sys.exit(1)

    print('Verbosity =', options.verbosity)
    logUtils.traceSetAt("ip.diffim", options.verbosity)

    snap1Exp = afwImage.ExposureF(options.s1)
    snap2Exp = afwImage.ExposureF(options.s2)

    config = ipDiffim.SnapPsfMatchTask.ConfigClass()
    config.kernel.name = "AL"
    subconfig = config.kernel.active

    snapDiff = subtractSnaps(snap1Exp, snap2Exp, subconfig, doWarping=options.warp)
    snapDiff.writeFits(options.sdiff)


def run():
    main()

if __name__ == '__main__':
    run()

# For debugging script:
#
# python examples/runSubtractSnaps.py --s1
# ~/LSST/becker_2012_0209_181253/update/calexp/v886894611-fr/R22/S11.fits
# --s2
# ~/LSST/becker_2012_0209_181253/update/calexp/v886264371-fr/R22/S11.fits
# --sdiff sdiff.fits --warp
