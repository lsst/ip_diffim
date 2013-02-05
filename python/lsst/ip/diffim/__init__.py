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

from .version import *

# c wrapper
from .diffimLib import *

# python code
from psfMatch import *
from imagePsfMatch import *
from modelPsfMatch import *
from snapPsfMatch import *
from diaSourceAnalysis import *
from makeKernelBasisList import *
from diaCatalogSourceSelector import *
from diffimTools import *

# automatically register ip_diffim Algorithms
import lsst.meas.algorithms
lsst.meas.algorithms.AlgorithmRegistry.register("centroid.dipole.naive", NaiveDipoleCentroidControl)
lsst.meas.algorithms.AlgorithmRegistry.register("flux.dipole.naive", NaiveDipoleFluxControl)
lsst.meas.algorithms.AlgorithmRegistry.register("flux.dipole.psf", PsfDipoleFluxControl)
del lsst # cleanup namespace
