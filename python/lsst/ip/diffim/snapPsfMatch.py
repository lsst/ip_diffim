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
import lsst.pex.config as pexConfig
from psfMatch import PsfMatchConfigDF, PsfMatchConfigAL
from imagePsfMatch import ImagePsfMatch

class SnapPsfMatchConfigDF(PsfMatchConfigDF):
    """Version of Psf Matching optimized for snap subtraction"""
    def __init__(self):
        PsfMatchConfigDF.__init__(self)
        
        # No spatial variation in model
        self.spatialKernelOrder = 0
        
        # Don't fit for differential background
        self.fitForBackground = False

        # Small kernel size
        self.kernelSize = 7

        # With zero spatial order don't worry about spatial clipping
        self.spatialKernelClipping = False

        # No regularization
        self.useRegularization = False
            
        
class SnapPsfMatchConfigAL(PsfMatchConfigAL):
    """Version of Psf Matching optimized for snap subtraction"""
    def __init__(self):
        PsfMatchConfigAL.__init__(self)
        
        # No spatial variation in model
        self.spatialKernelOrder = 0
        
        # Don't fit for differential background
        self.fitForBackground = False

        # Small kernel size
        self.kernelSize = 7

        # With zero spatial order don't worry about spatial clipping
        self.spatialKernelClipping = False

        # Simple basis set
        self.kernelBasisSet = "alard-lupton"
        self.alardNGauss = 2
        self.alardDegGauss = (4, 2)
        self.alardSigGauss = (1.0, 2.5)

class SnapPsfMatchConfig(pexConfig.Config):
    kernel = pexConfig.ConfigChoiceField(
        doc = "kernel type",
        typemap = dict(
            AL = SnapPsfMatchConfigAL,
            DF = SnapPsfMatchConfigDF
        ),
        default = "AL",
    )
            
class SnapPsfMatch(ImagePsfMatch):
    ConfigClass = SnapPsfMatchConfig
