# 
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
import numpy as np
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetect
import lsst.pipe.base as pipeBase
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConfig
import lsst.meas.algorithms as measAlg
import lsst.meas.deblender.baseline as deblendBaseline
from lsst.meas.algorithms import SourceMeasurementTask, SourceMeasurementConfig

class DipoleClassificationConfig(pexConfig.Config):
    minSn = pexConfig.Field(
        doc="Minimum quadrature sum of positive+negative lobe S/N to be considered a dipole",
        dtype=float, default=np.sqrt(2) * 5.0, 
        )
    maxFluxRatio = pexConfig.Field(
        doc = "Maximum flux ratio in either lobe to be considered a dipole",
        dtype = float, default = 0.65
        )

class DipoleMeasurementConfig(SourceMeasurementConfig):
    classification = pexConfig.ConfigField(
        dtype=DipoleClassificationConfig,
        doc="Dipole classification config"
        )

    def setDefaults(self):
        self.algorithms.names.add("centroid.dipole.naive")
        self.algorithms.names.add("flux.dipole.naive")
        self.algorithms.names.add("flux.dipole.psf")
        self.doReplaceWithNoise = False

class DipoleMeasurementTask(SourceMeasurementTask):
    ConfigClass = DipoleMeasurementConfig
    _DefaultName = "dipoleMeasurement"
    _ClassificationFlag = "classification.dipole"

    def __init__(self, schema, algMetadata=None, **kwds):
        """Create the task, adding necessary fields to the given schema.

        @param[in,out] schema        Schema object for measurement fields; modified in-place.
        @param[in,out] algMetadata   Passed to MeasureSources object to be filled with 
                                     metadata by algorithms (e.g. radii for aperture photometry).
        @param         **kwds        Passed to Task.__init__.
        """
        SourceMeasurementTask.__init__(self, schema, algMetadata, **kwds)
        self.dipoleAnalysis = DipoleAnalysis()

    @pipeBase.timeMethod
    def classify(self, sources):
        self.log.log(self.log.INFO, "Classifying %d sources" % len(sources))
        if not sources:
            return

        ctrl = self.config.classification
        try:
            key = sources[0].getSchema().find(self._ClassificationFlag).getKey()
        except:
            self.log.warn("Key %s not found in table" % (self._ClassificationFlag))            
            return

        for source in sources:
            passesSn = self.dipoleAnalysis.getSn(source) > ctrl.minSn

            negFlux   = np.abs(source.get("flux.dipole.psf.neg"))
            posFlux   = np.abs(source.get("flux.dipole.psf.pos"))
            totalFlux = negFlux + posFlux
            passesFluxNeg = (negFlux / (negFlux + posFlux)) < ctrl.maxFluxRatio
            passesFluxPos = (posFlux / (negFlux + posFlux)) < ctrl.maxFluxRatio

            if (passesSn and passesFluxPos and passesFluxNeg):
                val = 1.0
            else:
                val = 0.0

            source.set(key, val)

    def run(self, exposure, sources, **kwds):
        SourceMeasurementTask.run(self, exposure, sources, **kwds)
        self.classify(sources)

#########
# Other Support classs
#########

class SourceFlagChecker(object):
    """A functor to check whether a difference image source has any
    flags set that should cause it to be labeled bad."""
    def __init__(self, sources, badFlags=['flags.pixel.edge', 
                                          'flags.pixel.interpolated.center', 
                                          'flags.pixel.saturated.center']):
        self.keys = [sources.getSchema().find(name).key for name in badFlags]
        self.keys.append(sources.table.getCentroidFlagKey())

    def __call__(self, source):
        for k in self.keys:
            if source.get(k):
                return False
        return True

class DipoleAnalysis(object):
    """A functor to check for dipoles in difference image source tables."""
    def __init__(self):
        pass

    def __call__(self, source):
        return self.getSn(source), self.getCentroid(source), self.getOrientation(source)

    def getSn(self, source):
        posflux = source.get("flux.dipole.psf.pos")
        posfluxErr = source.get("flux.dipole.psf.pos.err")
        negflux = source.get("flux.dipole.psf.neg")
        negfluxErr = source.get("flux.dipole.psf.neg.err")

        # Not a dipole!
        if (posflux < 0) is (negflux < 0):
            return 0

        return np.sqrt((posflux/posfluxErr)**2 + (negflux/negfluxErr)**2)

    def getCentroid(self, source):
        negCen = source.get("flux.dipole.psf.neg.centroid")
        posCen = source.get("flux.dipole.psf.pos.centroid")
        if (False in np.isfinite(negCen)) or (False in np.isfinite(posCen)):
            return None
        
        center = afwGeom.Point2D(0.5*(negCen[0]+posCen[0]),
                                 0.5*(negCen[1]+posCen[1]))
        return center

    def getOrientation(self, source):
        negCen = source.get("flux.dipole.psf.neg.centroid")
        posCen = source.get("flux.dipole.psf.pos.centroid")
        if (False in np.isfinite(negCen)) or (False in np.isfinite(posCen)):
            return None

        dx, dy = posCen[0]-negCen[0], posCen[1]-negCen[1]
        angle  = afwGeom.Angle(np.arctan2(dx, dy), afwGeom.radians)
        return angle

    def displayDipoles(self, exposure, sources, frame=1):
        import lsst.afw.display.ds9 as ds9                
        ds9.mtv(exposure, frame=frame)
        with ds9.Buffering():
            for source in sources:
                cen = source.get("flux.dipole.psf.centroid")
                if (False in np.isfinite(cen)):
                    cen = source.getCentroid()

                isdipole = source.get("classification.dipole")
                if isdipole:
                    ctype= "green"
                    print "DIPOLE: %.1f,%.1f %.1f,%.1f" % (
                        cen.getX(), cen.getY(), 
                        source.get("flux.dipole.psf.neg"), 
                        source.get("flux.dipole.psf.pos"))
                else:
                    ctype = "red"
                    print "NOT DIPOLE: %.1f,%.1f" % (cen.getX(), cen.getY())

                ds9.dot("o", cen.getX(), cen.getY(), size=2, ctype=ctype, frame=frame)

                negCen = source.get("flux.dipole.psf.neg.centroid")
                posCen = source.get("flux.dipole.psf.pos.centroid")
                if (False in np.isfinite(negCen)) or (False in np.isfinite(posCen)):
                    continue

                ds9.line([(negCen.getX(), negCen.getY()),(posCen.getX(), posCen.getY())], 
                         ctype="yellow", frame=frame)
            
        


class DipoleDeblender(object):
    """A functor to deblend a source as a dipole, and return a new
       source with deblended footprints.  This necessarily overrides
       some of the functionality from
       meas_algorithms/python/lsst/meas/algorithms/deblend.py since we
       need a single source that contains the blended peaks, not
       multiple children sources.  This directly calls the core
       deblending code deblendBaseline.deblend (optionally _fitPsf for
       debugging).

       Not actively being used, but there is a unit test for it in 
       dipoleAlgorithm.py.
    """
    def __init__(self):
        # Set up defaults to send to deblender

        # Always deblend as Psf
        self.psfChisqCut1 = self.psfChisqCut2 = self.psfChisqCut2b = np.inf 
        self.log = pexLog.Log(pexLog.Log.getDefaultLog(), 
                              'lsst.ip.diffim.DipoleDeblender', pexLog.Log.INFO)
        self.sigma2fwhm = 2. * np.sqrt(2. * np.log(2.))

    def __call__(self, source, exposure):
        fp     = source.getFootprint()
        peaks  = fp.getPeaks()
        peaksF = [pk.getF() for pk in peaks]
        fbb    = fp.getBBox()
        fmask  = afwImage.MaskU(fbb)
        fmask.setXY0(fbb.getMinX(), fbb.getMinY())
        afwDetect.setMaskFromFootprint(fmask, fp, 1)

        psf        = exposure.getPsf()
        psfSigPix  = psf.computeShape().getDeterminantRadius()
        psfFwhmPix = psfSigPix * self.sigma2fwhm 
        subimage   = afwImage.ExposureF(exposure, fbb, True)
        cpsf       = deblendBaseline.CachingPsf(psf)

        # if fewer than 2 peaks, just return a copy of the source
        if len(peaks) < 2:
            return source.getTable().copyRecord(source)

        # make sure you only deblend 2 peaks; take the brighest and faintest
        speaks = [(p.getPeakValue(), p) for p in peaks]
        speaks.sort() 
        dpeaks = [speaks[0][1], speaks[-1][1]]

        # and only set these peaks in the footprint (peaks is mutable)
        peaks.clear()
        for peak in dpeaks:
            peaks.push_back(peak)

        if True:
            # Call top-level deblend task
            fpres = deblendBaseline.deblend(fp, exposure.getMaskedImage(), psf, psfFwhmPix,
                                            log = self.log,
                                            psfChisqCut1 = self.psfChisqCut1,
                                            psfChisqCut2 = self.psfChisqCut2,
                                            psfChisqCut2b = self.psfChisqCut2b)
        else:
            # Call lower-level _fit_psf task

            # Prepare results structure
            fpres = deblendBaseline.PerFootprint()
            fpres.peaks = []
            for pki,pk in enumerate(dpeaks):
                pkres = deblendBaseline.PerPeak()
                pkres.peak = pk
                pkres.pki = pki
                fpres.peaks.append(pkres)
            
            for pki,(pk,pkres,pkF) in enumerate(zip(dpeaks, fpres.peaks, peaksF)):
                self.log.logdebug('Peak %i' % pki)
                deblendBaseline._fitPsf(fp, fmask, pk, pkF, pkres, fbb, dpeaks, peaksF, self.log, 
                                         cpsf, psfFwhmPix, 
                                         subimage.getMaskedImage().getImage(), 
                                         subimage.getMaskedImage().getVariance(), 
                                         self.psfChisqCut1, self.psfChisqCut2, self.psfChisqCut2b)


        deblendedSource = source.getTable().copyRecord(source)
        deblendedSource.setParent(source.getId())
        peakList        = deblendedSource.getFootprint().getPeaks()
        peakList.clear()

        for i, peak in enumerate(fpres.peaks):
            if peak.psfFitFlux > 0:
                suffix = "pos"
            else:
                suffix = "neg"
            c = peak.psfFitCenter
            self.log.info("deblended.centroid.dipole.psf.%s %f %f" % (
                suffix, c[0], c[1]))
            self.log.info("deblended.chi2dof.dipole.%s %f" % (
                suffix, peak.psfFitChisq / peak.psfFitDof))
            self.log.info("deblended.flux.dipole.psf.%s %f" % (
                suffix, peak.psfFitFlux * np.sum(peak.templateMaskedImage.getImage().getArray())))
            peakList.push_back(peak.peak)
        return deblendedSource
        
