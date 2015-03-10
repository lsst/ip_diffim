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
import numpy as np
import lsst.pex.config as pexConfig
import lsst.afw.display.ds9 as ds9
import lsst.meas.algorithms as measAlg
import lsst.pex.logging as pexLog

class DiaCatalogSourceSelectorConfig(pexConfig.Config):
    # Selection cuts on the input source catalog
    fluxLim = pexConfig.Field(
        doc = "specify the minimum psfFlux for good Kernel Candidates",
        dtype = float,
        default = 0.0,
        check = lambda x: x >= 0.0,
    )
    fluxMax = pexConfig.Field(
        doc = "specify the maximum psfFlux for good Kernel Candidates (ignored if == 0)",
        dtype = float,
        default = 0.0,
        check = lambda x: x >= 0.0,
    )
    badPixelFlags = pexConfig.ListField(
        doc = "Kernel candidate objects may not have any of these bits set",
        dtype = str,
        default = ["flags.pixel.edge", "flags.pixel.interpolated.center", "flags.pixel.saturated.center", "flags.badcentroid"],
        )
    # Selection cuts on the reference catalog
    selectStar = pexConfig.Field(
        doc = "Select objects that are flagged as stars",
        dtype = bool,
        default = True
    )
    selectGalaxy = pexConfig.Field(
        doc = "Select objects that are flagged as galaxies",
        dtype = bool,
        default = False
    )
    includeVariable = pexConfig.Field(
        doc = "Include objects that are known to be variable",
        dtype = bool,
        default = False
    )
    grMin = pexConfig.Field(
        doc = "Minimum g-r color for selection (inclusive)",
        dtype = float,
        default = 0.0
    )
    grMax = pexConfig.Field(
        doc = "Maximum g-r color for selection (inclusive)",
        dtype = float,
        default = 3.0
    )

class CheckSource(object):
    """A functor to check whether a source has any flags set that should cause it to be labeled bad."""

    def __init__(self, table, fluxLim, fluxMax, badPixelFlags):
        self.keys = [table.getSchema().find(name).key for name in badPixelFlags]
        self.fluxLim = fluxLim
        self.fluxMax = fluxMax

    def __call__(self, source):
        for k in self.keys:
            if source.get(k):
                return False
        if self.fluxLim != None and source.getPsfFlux() < self.fluxLim: # ignore faint objects
            return False
        if self.fluxMax != 0.0 and source.getPsfFlux() > self.fluxMax: # ignore bright objects
            return False
        return True

class DiaCatalogSourceSelector(object):
    ConfigClass = DiaCatalogSourceSelectorConfig

    def __init__(self, config=None):
        """Construct a source selector that uses a reference catalog
        
        @param[in] config: An instance of ConfigClass
        """
        if not config:
            config = DiaCatalogSourceSelector.ConfigClass()
        self.config = config
        self.log = pexLog.Log(pexLog.Log.getDefaultLog(),
                              'lsst.ip.diffim.DiaCatalogSourceSelector', pexLog.Log.INFO)

    def selectSources(self, exposure, sources, matches=None):
        """Return a list of Sources for Kernel candidates 
        
        @param[in] exposure: the exposure containing the sources
        @param[in] sources: a source list containing sources that may be candidates
        @param[in] matches: a match vector as produced by meas_astrom; not optional
                            (passing None just allows us to handle the exception better here
                            than in calling code)
        
        @return kernelCandidateSourceList: a list of sources to be used as kernel candidates
 
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayExposure = lsstDebug.Info(__name__).displayExposure
        pauseAtEnd = lsstDebug.Info(__name__).pauseAtEnd

        if matches is None:
            raise RuntimeError(
                "Cannot use catalog source selector without running astrometry."
                )

        mi = exposure.getMaskedImage()
        
        if display:
            if displayExposure:
                ds9.mtv(mi, title="Kernel candidates", frame=lsstDebug.frame)
        #
        # Look for flags in each Source
        #
        isGoodSource = CheckSource(sources, self.config.fluxLim, self.config.fluxMax, self.config.badPixelFlags)

        #
        # Go through and find all the acceptable candidates in the catalogue
        #
        kernelCandidateSourceList = []

        doColorCut = True
        with ds9.Buffering():
            for ref, source, d in matches:
                if not isGoodSource(source):
                    symb, ctype = "+", ds9.RED
                else:
                    isStar = ref.get("stargal")
                    isVar = not ref.get("photometric")
                    gMag = None
                    rMag = None
                    if doColorCut:
                        try:
                            gMag = -2.5 * np.log10(ref.get("g"))
                            rMag = -2.5 * np.log10(ref.get("r"))
                        except KeyError:
                            self.log.warn("Cannot cut on color info; fields 'g' and 'r' do not exist")
                            doColorCut = False
                            isRightColor = True
                        else:
                            isRightColor = (gMag-rMag) >= self.config.grMin and (gMag-rMag) <= self.config.grMax
                        
                    isRightType  = (self.config.selectStar and isStar) or (self.config.selectGalaxy and not isStar)
                    isRightVar   = (self.config.includeVariable) or (self.config.includeVariable is isVar)
                    if isRightType and isRightVar and isRightColor:
                        kernelCandidateSourceList.append(source)
                        symb, ctype = "+", ds9.GREEN
                    else:
                        symb, ctype = "o", ds9.BLUE

                if display and displayExposure:
                    ds9.dot(symb, source.getX() - mi.getX0(), source.getY() - mi.getY0(),
                            size=4, ctype=ctype, frame=lsstDebug.frame)

        if display:
            lsstDebug.frame += 1
            if pauseAtEnd:
                raw_input("Continue? y[es] p[db] ")

        return kernelCandidateSourceList

measAlg.starSelectorRegistry.register("diacatalog", DiaCatalogSourceSelector)

