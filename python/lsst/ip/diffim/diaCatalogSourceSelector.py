#
# LSST Data Management System
# Copyright 2008-2016 LSST Corporation.
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

import numpy as np

from lsst.afw.table import SourceCatalog
from lsst.pipe.base import Struct
import lsst.pex.config as pexConfig
import lsst.afw.display.ds9 as ds9
import lsst.meas.algorithms as measAlg

__all__ = ["DiaCatalogSourceSelectorConfig", "DiaCatalogSourceSelectorTask"]


class DiaCatalogSourceSelectorConfig(measAlg.BaseStarSelectorConfig):
    # Selection cuts on the input source catalog
    fluxLim = pexConfig.Field(
        doc="specify the minimum psfFlux for good Kernel Candidates",
        dtype=float,
        default=0.0,
        check=lambda x: x >= 0.0,
    )
    fluxMax = pexConfig.Field(
        doc="specify the maximum psfFlux for good Kernel Candidates (ignored if == 0)",
        dtype=float,
        default=0.0,
        check=lambda x: x >= 0.0,
    )
    # Selection cuts on the reference catalog
    selectStar = pexConfig.Field(
        doc="Select objects that are flagged as stars",
        dtype=bool,
        default=True
    )
    selectGalaxy = pexConfig.Field(
        doc="Select objects that are flagged as galaxies",
        dtype=bool,
        default=False
    )
    includeVariable = pexConfig.Field(
        doc="Include objects that are known to be variable",
        dtype=bool,
        default=False
    )
    grMin = pexConfig.Field(
        doc="Minimum g-r color for selection (inclusive)",
        dtype=float,
        default=0.0
    )
    grMax = pexConfig.Field(
        doc="Maximum g-r color for selection (inclusive)",
        dtype=float,
        default=3.0
    )

    def setDefaults(self):
        measAlg.BaseStarSelectorConfig.setDefaults(self)
        self.badFlags = [
            "base_PixelFlags_flag_edge",
            "base_PixelFlags_flag_interpolatedCenter",
            "base_PixelFlags_flag_saturatedCenter",
            "slot_Centroid_flag",
        ]


class CheckSource(object):
    """A functor to check whether a source has any flags set that should cause it to be labeled bad."""

    def __init__(self, table, fluxLim, fluxMax, badFlags):
        self.keys = [table.getSchema().find(name).key for name in badFlags]
        self.fluxLim = fluxLim
        self.fluxMax = fluxMax

    def __call__(self, source):
        for k in self.keys:
            if source.get(k):
                return False
        if self.fluxLim is not None and source.getPsfFlux() < self.fluxLim:  # ignore faint objects
            return False
        if self.fluxMax != 0.0 and source.getPsfFlux() > self.fluxMax:  # ignore bright objects
            return False
        return True

## @addtogroup LSST_task_documentation
## @{
## @page DiaCatalogSourceSelectorTask
## @ref DiaCatalogSourceSelectorTask_ "DiaCatalogSourceSelectorTask"
## @copybrief DiaCatalogSourceSelectorTask
## @}


class DiaCatalogSourceSelectorTask(measAlg.BaseStarSelectorTask):
    """!Select sources for Kernel candidates

    @anchor DiaCatalogSourceSelectorTask_

    @section ip_diffim_diaCatalogSourceSelector_Contents  Contents

     - @ref ip_diffim_diaCatalogSourceSelector_Purpose
     - @ref ip_diffim_diaCatalogSourceSelector_Initialize
     - @ref ip_diffim_diaCatalogSourceSelector_IO
     - @ref ip_diffim_diaCatalogSourceSelector_Config
     - @ref ip_diffim_diaCatalogSourceSelector_Debug

    @section ip_diffim_diaCatalogSourceSelector_Purpose  Description

    A naive star selector based on second moments. Use with caution.

    @section ip_diffim_diaCatalogSourceSelector_Initialize  Task initialisation

    @copydoc \_\_init\_\_

    @section ip_diffim_diaCatalogSourceSelector_IO  Invoking the Task

    Like all star selectors, the main method is `run`.

    @section ip_diffim_diaCatalogSourceSelector_Config  Configuration parameters

    See @ref DiaCatalogSourceSelectorConfig

    @section ip_diffim_diaCatalogSourceSelector_Debug  Debug variables

    DiaCatalogSourceSelectorTask has a debug dictionary with the following keys:
    <dl>
    <dt>display
    <dd>bool; if True display debug information
    <dt>displayExposure
    <dd>bool; if True display exposure
    <dt>pauseAtEnd
    <dd>bool; if True wait after displaying everything and wait for user input
    </dl>

    For example, put something like:
    @code{.py}
        import lsstDebug
        def DebugInfo(name):
            di = lsstDebug.getInfo(name)  # N.b. lsstDebug.Info(name) would call us recursively
            if name.endswith("catalogStarSelector"):
                di.display = True

            return di

        lsstDebug.Info = DebugInfo
    @endcode
    into your `debug.py` file and run your task with the `--debug` flag.
    """
    ConfigClass = DiaCatalogSourceSelectorConfig
    usesMatches = True  # selectStars uses (requires) its matches argument

    def selectStars(self, exposure, sourceCat, matches=None):
        """Select sources for Kernel candidates

        @param[in] exposure  the exposure containing the sources
        @param[in] sourceCat  catalog of sources that may be stars (an lsst.afw.table.SourceCatalog)
        @param[in] matches  a match vector as produced by meas_astrom; required
                            (defaults to None to match the StarSelector API and improve error handling)

        @return an lsst.pipe.base.Struct containing:
        - starCat  a list of sources to be used as kernel candidates
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayExposure = lsstDebug.Info(__name__).displayExposure
        pauseAtEnd = lsstDebug.Info(__name__).pauseAtEnd

        if matches is None:
            raise RuntimeError("DiaCatalogSourceSelector requires matches")

        mi = exposure.getMaskedImage()

        if display:
            if displayExposure:
                ds9.mtv(mi, title="Kernel candidates", frame=lsstDebug.frame)
        #
        # Look for flags in each Source
        #
        isGoodSource = CheckSource(sourceCat, self.config.fluxLim, self.config.fluxMax, self.config.badFlags)

        #
        # Go through and find all the acceptable candidates in the catalogue
        #
        starCat = SourceCatalog(sourceCat.schema)

        if display and displayExposure:
            symbs = []
            ctypes = []

        doColorCut = True

        refSchema = matches[0][0].schema
        rRefFluxField = measAlg.getRefFluxField(refSchema, "r")
        gRefFluxField = measAlg.getRefFluxField(refSchema, "g")
        for ref, source, d in matches:
            if not isGoodSource(source):
                if display and displayExposure:
                    symbs.append("+")
                    ctypes.append(ds9.RED)
            else:
                isStar = not ref.get("resolved")
                isVar = not ref.get("photometric")
                gMag = None
                rMag = None
                if doColorCut:
                    try:
                        gMag = -2.5 * np.log10(ref.get(gRefFluxField))
                        rMag = -2.5 * np.log10(ref.get(rRefFluxField))
                    except KeyError:
                        self.log.warn("Cannot cut on color info; fields 'g' and 'r' do not exist")
                        doColorCut = False
                        isRightColor = True
                    else:
                        isRightColor = (gMag-rMag) >= self.config.grMin and (gMag-rMag) <= self.config.grMax

                isRightType = (self.config.selectStar and isStar) or (self.config.selectGalaxy and not isStar)
                isRightVar = (self.config.includeVariable) or (self.config.includeVariable is isVar)
                if isRightType and isRightVar and isRightColor:
                    starCat.append(source)
                    if display and displayExposure:
                        symbs.append("+")
                        ctypes.append(ds9.GREEN)
                elif display and displayExposure:
                    symbs.append("o")
                    ctypes.append(ds9.BLUE)

        if display and displayExposure:
            with ds9.Buffering():
                for (ref, source, d), symb, ctype in zip(matches, symbs, ctypes):
                    if display and displayExposure:
                        ds9.dot(symb, source.getX() - mi.getX0(), source.getY() - mi.getY0(),
                                size=4, ctype=ctype, frame=lsstDebug.frame)

        if display:
            lsstDebug.frame += 1
            if pauseAtEnd:
                input("Continue? y[es] p[db] ")

        return Struct(
            starCat=starCat,
        )


measAlg.starSelectorRegistry.register("diacatalog", DiaCatalogSourceSelectorTask)
