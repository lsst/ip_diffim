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

import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim
from lsst.pex.logging import Log, Trace


def subtractSnaps(snap1, snap2, policy, doWarping = False):
    psfmatch = ipDiffim.ImagePsfMatch(policy)
    results  = psfmatch.subtractExposures(snap1, snap2, doWarping = doWarping)
    snapDiff, kernelModel, bgModel, kernelCellSet = results
    return snapDiff
    
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
    parser.add_option('--policy', help='user override policy file')
    parser.add_option('--warp', action='store_true', default=False, help='astrometrically warp snap1')
    parser.add_option('-v', '--verbosity', type=int, default=defVerbosity,
                      help='verbosity of Trace messages')

    (options, args) = parser.parse_args()
    if options.s1 == None or options.s2 == None or options.sdiff == None:
        parser.print_help()
        sys.exit(1)

    print 'Verbosity =', options.verbosity
    Trace.setVerbosity('lsst.ip.diffim', options.verbosity)
         
    snap1Exp   = afwImage.ExposureF(options.s1)
    snap2Exp   = afwImage.ExposureF(options.s2)
    policy     = ipDiffim.makeDefaultPolicy(mergePolicy = options.policy)
    snapPolicy = ipDiffim.modifyForSnapSubtraction(policy)

    snapDiff   = subtractSnaps(snap1Exp, snap2Exp, snapPolicy, doWarping = options.warp)
    snapDiff.writeFits(options.sdiff)
    
def run():
    Log.getDefaultLog()
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), 'Objects leaked:'
        print dafBase.Citizen_census(dafBase.cout, memId0)

if __name__ == '__main__':
    run()

# For debugging script:
# 
# python examples/subtractSnaps.py --s1    /home/becker/LSST/PT1/psfMatch/wp_trunk_2011_0601_171743/update/calexp/v885335911-fr/R22/S11.fits --s2  /home/becker/LSST/PT1/psfMatch/wp_trunk_2011_0601_171743/update/calexp/v886257211-fr/R22/S11.fits --sdiff sdiff.fits --warp
