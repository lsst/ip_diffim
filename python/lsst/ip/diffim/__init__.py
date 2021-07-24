#
# LSST Data Management System
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
# See the COPYRIGHT file
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#

# C++ wrapper
# hoist symbols lsst.ip.diffim.detail up into lsst.ip.diffim
from .diffimLib import *

# Python code
from .dcrModel import *
from .psfMatch import *
from .imagePsfMatch import *
from .modelPsfMatch import *
from .snapPsfMatch import *
from .makeKernelBasisList import *
from .diaCatalogSourceSelector import *
from .dipoleMeasurement import *
from .diffimTools import *
from .kernelCandidateQa import *
from .getTemplate import *
from .diaCatalogSourceSelector import *
from lsst.meas.base import wrapSimpleAlgorithm
from .dipoleFitTask import *
from .imageDecorrelation import *
from .imageMapReduce import *
from .zogy import *
from .version import *

# automatically register ip_diffim Algorithms
wrapSimpleAlgorithm(NaiveDipoleCentroid, Control=DipoleCentroidControl, executionOrder=0.0)
wrapSimpleAlgorithm(NaiveDipoleFlux, Control=DipoleFluxControl, executionOrder=2.0)
wrapSimpleAlgorithm(PsfDipoleFlux, Control=PsfDipoleFluxControl, executionOrder=2.0)

