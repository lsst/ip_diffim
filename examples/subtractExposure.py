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

import eups
import sys
import os
import optparse

import lsst.daf.base as dafBase
import lsst.afw.image as afwImage

from lsst.pex.logging import Trace
from lsst.pex.logging import Log

from lsst.ip.diffim import subtractExposure, generateDefaultPolicy
import lsst.ip.diffim.diffimTools as diffimTools

def main():
    defDataDir   = eups.productDir('afwdata')
    if defDataDir == None:
        print 'Error: afwdata not set up'
        sys.exit(1)
    imageProcDir = eups.productDir('ip_diffim')
    if imageProcDir == None:
        print 'Error: could not set up ip_diffim'
        sys.exit(1)

    defSciencePath  = os.path.join(defDataDir, 'CFHT', 'D4', 'cal-53535-i-797722_1')
    defTemplatePath = os.path.join(defDataDir, 'CFHT', 'D4', 'cal-53535-i-797722_1_tmpl')
    defPolicyPath   = os.path.join(imageProcDir, 'policy', 'ImageSubtractStageDictionary.paf')
    defOutputPath   = 'diffExposure'
    defVerbosity    = 0
    
    usage = """usage: %%prog [options] [scienceExposure [templateExposure [outputExposure]]]]

Notes:
- image arguments are paths to Expoure fits files
- image arguments must NOT include the final _img.fits
- the result is science image - template image
- the template exposure is convolved, the science exposure is not
- default scienceExposure=%s
- default templateExposure=%s
- default outputExposure=%s 
- default --policy=%s
""" % (defSciencePath, defTemplatePath, defOutputPath, defPolicyPath)
    
    parser = optparse.OptionParser(usage)
    parser.add_option('-p', '--policy', default=defPolicyPath, help='policy file')
    parser.add_option('-v', '--verbosity', type=int, default=defVerbosity,
                      help='verbosity of Trace messages')
    parser.add_option('-d', '--display', action='store_true', default=False,
                      help='display the images')
    parser.add_option('-b', '--bg', action='store_true', default=False,
                      help='subtract backgrounds')

    (options, args) = parser.parse_args()
    
    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue
    
    sciencePath  = getArg(0, defSciencePath)
    templatePath = getArg(1, defTemplatePath)
    outputPath   = getArg(2, defOutputPath)
    policyPath   = options.policy
    
    print 'Science exposure: ', sciencePath
    print 'Template exposure:', templatePath
    print 'Output exposure:  ', outputPath
    print 'Policy file:      ', policyPath
    
    templateExposure = afwImage.ExposureF(templatePath)
    scienceExposure  = afwImage.ExposureF(sciencePath)
    policy           = generateDefaultPolicy(policyPath)

    display = False
    if options.display:
        print 'Display =', options.display
        display = True

    bgSub = False
    if options.bg:
        print 'Background subtract =', options.bg
        bgSub = True

    if options.verbosity > 0:
        print 'Verbosity =', options.verbosity
        Trace.setVerbosity('lsst.ip.diffim', options.verbosity)

    if bgSub:
        diffimTools.backgroundSubtract(policy, [templateExposure.getMaskedImage(),
                                                scienceExposure.getMaskedImage()])

    results = subtractExposure(templateExposure,
                               scienceExposure,
                               policy,
                               display)
    differenceExposure = results[0]
    differenceExposure.writeFits(outputPath)

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
