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
from .computeSpatiallySampledMetrics import *
from .dcrModel import *
from .psfMatch import *
from .modelPsfMatch import *
from .makeKernel import *
from .makeKernelBasisList import *
from .dipoleMeasurement import *
from .getTemplate import *
from .dipoleFitTask import *
from .detectAndMeasure import *
from .imageDecorrelation import *
from .imageMapReduce import *
from .scoreMeasurement import *
from .subtractImages import *
from .version import *

# automatically register ip_diffim Algorithms;
# CENTROID_ORDER=0.0, FLUX_ORDER==2.0
from lsst.meas.base import wrapSimpleAlgorithm
wrapSimpleAlgorithm(PsfDipoleFlux, Control=PsfDipoleFluxControl, executionOrder=2.0)
