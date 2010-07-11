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

import lsst.afw.image as afwImage
import lsst.afw.math  as afwMath
import sys

tRoot  = sys.argv[1]
iRoot  = sys.argv[2]
kernel = sys.argv[3]

tMi   = afwImage.MaskedImageF(tRoot)
iMi   = afwImage.MaskedImageF(iRoot)
kImg  = afwImage.ImageD(kernel)
k     = afwMath.FixedKernel(kImg)

cMi   = afwImage.MaskedImageF(tMi.getDimensions())
afwMath.convolve(cMi, tMi, k, False)

bbox      = afwImage.BBox(afwImage.PointI(k.getCtrX(),
                                          k.getCtrY()) ,
                          afwImage.PointI(cMi.getWidth() - (k.getWidth() - k.getCtrX()),
                                          cMi.getHeight() - (k.getHeight() - k.getCtrY())))
cMi2 = afwImage.MaskedImageF(cMi, bbox)
cMi2.writeFits('conv')
