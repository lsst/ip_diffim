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

import time
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.logging as pexLog

def warpTemplateExposure(templateExposure, scienceExposure, policy):
    t0 = time.time()

    # The destination Wcs is in scienceExposure
    # Create the warping Kernel according to policy
    warpingKernelName = policy.getString("warpingKernelName")
    warpingKernel     = afwMath.makeWarpingKernel(warpingKernelName)

    # create a blank exposure to hold the remapped template exposure
    remappedTemplateExposure = templateExposure.Factory(
        scienceExposure.getBBox(afwImage.PARENT), 
        scienceExposure.getWcs())

    # warp the template exposure
    afwMath.warpExposure(remappedTemplateExposure, 
                         templateExposure, 
                         warpingKernel)
        
    t1 = time.time()
    pexLog.Trace("lsst.ip.diffim.warpTemplateExposure", 1,
                 "Total time for warping : %.2f s" % (t1-t0))

    return remappedTemplateExposure
    
