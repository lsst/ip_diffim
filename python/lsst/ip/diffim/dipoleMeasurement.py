# 
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#

"""Dipole measurement plugins and analysis tools.

Configs for measurement on image differences are controlled by
imageDifferenceTask (processImageDifferenceTask)

Note: If we need Dipole Measurement to do things that
cannot be expressed as plugins, a lightweight wrapper class
DipoleMeasurementTask(lsst.pipe.base.Task) which runs
SingleFrameMeasurementTask may be added to this module.
"""

import numpy as np
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetect
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConfig
import lsst.meas.deblender.baseline as deblendBaseline
from lsst.meas.base.pluginRegistry import register
from lsst.meas.base import SingleFramePluginConfig, SingleFramePlugin, SingleFrameMeasurementConfig
import lsst.afw.display.ds9 as ds9

__all__ = ("DipoleAnalysis", "DipoleDeblender", "makeDipoleMeasurementConfig",
           "SourceFlagChecker", "ClassificationDipoleConfig", "ClassificationDipolePlugin")


def makeDipoleMeasurementConfig(config=None):
    """Return default config for dipole measurement
    """
    if config is None:
        config = SingleFrameMeasurementConfig()

    config.plugins.names = ["base_PsfFlux",
                            "ip_diffim_PsfDipoleFlux",
                            "ip_diffim_NaiveDipoleFlux",
                            "ip_diffim_NaiveDipoleCentroid",
                            "ip_diffim_ClassificationDipole",
                            "base_CircularApertureFlux",
                            "base_SkyCoord"]

    config.slots.calibFlux = None
    config.slots.modelFlux = None
    config.slots.instFlux = None
    config.slots.shape = None
    config.slots.centroid = "ip_diffim_NaiveDipoleCentroid"
    config.doReplaceWithNoise = False
    return config


class ClassificationDipoleConfig(SingleFramePluginConfig):
    """Configuration for classification of detected diaSources as dipole or not"""
    minSn = pexConfig.Field(
        doc="Minimum quadrature sum of positive+negative lobe S/N to be considered a dipole",
        dtype=float, default=np.sqrt(2) * 5.0,
        )
    maxFluxRatio = pexConfig.Field(
        doc="Maximum flux ratio in either lobe to be considered a dipole",
        dtype = float, default = 0.65
        )


@register("ip_diffim_ClassificationDipole")
class ClassificationDipolePlugin(SingleFramePlugin):
    """
    A binary measure of the whether a diasource is a dipole.
    """

    ConfigClass = ClassificationDipoleConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.CLASSIFY_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.dipoleAnalysis = DipoleAnalysis()
        self.keyProbability = schema.addField(name + "_value", type="D",
                                              doc="Set to 1 for dipoles, else 0.")
        self.keyFlag = schema.addField(name + "_flag", type="Flag", doc="Set to 1 for any fatal failure.")

    def measure(self, measRecord, exposure):
            passesSn = self.dipoleAnalysis.getSn(measRecord) > self.config.minSn
            negFlux = np.abs(measRecord.get("ip_diffim_PsfDipoleFlux_neg_flux"))
            negFluxFlag = measRecord.get("ip_diffim_PsfDipoleFlux_neg_flag")
            posFlux = np.abs(measRecord.get("ip_diffim_PsfDipoleFlux_pos_flux"))
            posFluxFlag = measRecord.get("ip_diffim_PsfDipoleFlux_pos_flag")

            if negFluxFlag or posFluxFlag:
                self.fail(measRecord)
                # continue on to classify

            totalFlux = negFlux + posFlux

            # If negFlux or posFlux are NaN, these evaluate to False
            passesFluxNeg = (negFlux / totalFlux) < self.config.maxFluxRatio
            passesFluxPos = (posFlux / totalFlux) < self.config.maxFluxRatio
            if (passesSn and passesFluxPos and passesFluxNeg):
                val = 1.0
            else:
                val = 0.0

            measRecord.set(self.keyProbability, val)

    def fail(self, measRecord, error=None):
        measRecord.set(self.keyFlag, True)


class SourceFlagChecker(object):
    """!Functor class to check whether a diaSource has flags set that should cause it to be labeled bad."""
    def __init__(self, sources, badFlags=None):
        """!Constructor

        @param sources     Sources that will be measured
        @param badFlags    A list of flags that will be used to determine if there was a measurement problem

        The list of badFlags will be used to make a list of keys to check for measurement flags on.  By
        default the centroid keys are added to this list"""

        self.badFlags = ['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter']
        if badFlags is not None:
            for flag in badFlags:
                self.badFlags.append(flag)
        self.keys = [sources.getSchema().find(name).key for name in self.badFlags]
        self.keys.append(sources.table.getCentroidFlagKey())

    def __call__(self, source):
        """!Call the source flag checker on a single Source

        @param source      Source that will be examined"""
        for k in self.keys:
            if source.get(k):
                return False
        return True

class DipoleAnalysis(object):
    """!Functor class that provides (S/N, position, orientation) of measured dipoles"""
    def __init__(self):
        """!Constructor"""
        pass

    def __call__(self, source):
        """!Parse information returned from dipole measurement

        @param source  The source that will be examined"""
        return self.getSn(source), self.getCentroid(source), self.getOrientation(source)

    def getSn(self, source):
        """!Get the total signal-to-noise of the dipole; total S/N is from positive and negative lobe

        @param source  The source that will be examined"""

        posflux = source.get("ip_diffim_PsfDipoleFlux_pos_flux")
        posfluxErr = source.get("ip_diffim_PsfDipoleFlux_pos_fluxSigma")
        negflux = source.get("ip_diffim_PsfDipoleFlux_neg_flux")
        negfluxErr = source.get("ip_diffim_PsfDipoleFlux_neg_fluxSigma")

        # Not a dipole!
        if (posflux < 0) is (negflux < 0):
            return 0

        return np.sqrt((posflux/posfluxErr)**2 + (negflux/negfluxErr)**2)

    def getCentroid(self, source):
        """!Get the centroid of the dipole; average of positive and negative lobe

        @param source  The source that will be examined"""

        negCenX = source.get("ip_diffim_PsfDipoleFlux_neg_centroid_x")
        negCenY = source.get("ip_diffim_PsfDipoleFlux_neg_centroid_y")
        posCenX = source.get("ip_diffim_PsfDipoleFlux_pos_centroid_x")
        posCenY = source.get("ip_diffim_PsfDipoleFlux_pos_centroid_y")
        if (np.isinf(negCenX) or np.isinf(negCenY) or np.isinf(posCenX) or np.isinf(posCenY)):
            return None

        center = afwGeom.Point2D(0.5*(negCenX+posCenX),
                                 0.5*(negCenY+posCenY))
        return center

    def getOrientation(self, source):
        """!Calculate the orientation of dipole; vector from negative to positive lobe

        @param source  The source that will be examined"""

        negCenX = source.get("ip_diffim_PsfDipoleFlux_neg_centroid_x")
        negCenY = source.get("ip_diffim_PsfDipoleFlux_neg_centroid_y")
        posCenX = source.get("ip_diffim_PsfDipoleFlux_pos_centroid_x")
        posCenY = source.get("ip_diffim_PsfDipoleFlux_pos_centroid_y")
        if (np.isinf(negCenX) or np.isinf(negCenY) or np.isinf(posCenX) or np.isinf(posCenY)):
            return None

        dx, dy = posCenX-negCenX, posCenY-negCenY
        angle  = afwGeom.Angle(np.arctan2(dx, dy), afwGeom.radians)
        return angle

    def displayDipoles(self, exposure, sources):
        """!Display debugging information on the detected dipoles

        @param exposure  Image the dipoles were measured on
        @param sources   The set of diaSources that were measured"""

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayDiaSources = lsstDebug.Info(__name__).displayDiaSources
        maskTransparency = lsstDebug.Info(__name__).maskTransparency
        if not maskTransparency:
            maskTransparency = 90
        ds9.setMaskTransparency(maskTransparency)
        ds9.mtv(exposure, frame=lsstDebug.frame)

        if display and displayDiaSources:
            with ds9.Buffering():
                for source in sources:
                    cenX, cenY = source.get("ipdiffim_DipolePsfFlux_centroid")
                    if np.isinf(cenX) or np.isinf(cenY):
                        cenX, cenY = source.getCentroid()

                    isdipole = source.get("classification.dipole")
                    if isdipole and np.isfinite(isdipole):
                        # Dipole
                        ctype= "green"
                    else:
                        # Not dipole
                        ctype = "red"

                    ds9.dot("o", cenX, cenY, size=2, ctype=ctype, frame=lsstDebug.frame)

                    negCenX = source.get("ip_diffim_PsfDipoleFlux_neg_centroid_x")
                    negCenY = source.get("ip_diffim_PsfDipoleFlux_neg_centroid_y")
                    posCenX = source.get("ip_diffim_PsfDipoleFlux_pos_centroid_x")
                    posCenY = source.get("ip_diffim_PsfDipoleFlux_pos_centroid_y")
                    if (np.isinf(negCenX) or np.isinf(negCenY) or np.isinf(posCenX) or np.isinf(posCenY)):
                        continue

                    ds9.line([(negCenX, negCenY), (posCenX, posCenY)], ctype="yellow", frame=lsstDebug.frame)

            lsstDebug.frame += 1



class DipoleDeblender(object):
    """!Functor to deblend a source as a dipole, and return a new source with deblended footprints.

       This necessarily overrides some of the functionality from
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
            peaks.append(peak)

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
                suffix, peak.psfFitFlux * np.sum(peak.templateImage.getArray())))
            peakList.append(peak.peak)
        return deblendedSource
