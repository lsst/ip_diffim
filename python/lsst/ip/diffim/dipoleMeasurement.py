# This file is part of ip_diffim.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np

import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.meas.deblender.baseline as deblendBaseline
from lsst.meas.base.pluginRegistry import register
from lsst.meas.base import SingleFrameMeasurementTask, SingleFrameMeasurementConfig, \
    SingleFramePluginConfig, SingleFramePlugin
import lsst.afw.display as afwDisplay
from lsst.utils.logging import getLogger

__all__ = ("DipoleMeasurementConfig", "DipoleMeasurementTask", "DipoleAnalysis", "DipoleDeblender",
           "SourceFlagChecker", "ClassificationDipoleConfig", "ClassificationDipolePlugin")


class ClassificationDipoleConfig(SingleFramePluginConfig):
    """Configuration for classification of detected diaSources as dipole or not"""
    minSn = pexConfig.Field(
        doc="Minimum quadrature sum of positive+negative lobe S/N to be considered a dipole",
        dtype=float, default=np.sqrt(2) * 5.0,
    )
    maxFluxRatio = pexConfig.Field(
        doc="Maximum flux ratio in either lobe to be considered a dipole",
        dtype=float, default=0.65
    )


@register("ip_diffim_ClassificationDipole")
class ClassificationDipolePlugin(SingleFramePlugin):
    """A plugin to classify whether a diaSource is a dipole.
    """

    ConfigClass = ClassificationDipoleConfig

    @classmethod
    def getExecutionOrder(cls):
        """
        Returns
        -------
        result : `callable`
        """
        return cls.APCORR_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.dipoleAnalysis = DipoleAnalysis()
        self.keyProbability = schema.addField(name + "_value", type="D",
                                              doc="Set to 1 for dipoles, else 0.")
        self.keyFlag = schema.addField(name + "_flag", type="Flag", doc="Set to 1 for any fatal failure.")

    def measure(self, measRecord, exposure):
        passesSn = self.dipoleAnalysis.getSn(measRecord) > self.config.minSn
        negFlux = np.abs(measRecord.get("ip_diffim_PsfDipoleFlux_neg_instFlux"))
        negFluxFlag = measRecord.get("ip_diffim_PsfDipoleFlux_neg_flag")
        posFlux = np.abs(measRecord.get("ip_diffim_PsfDipoleFlux_pos_instFlux"))
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


class DipoleMeasurementConfig(SingleFrameMeasurementConfig):
    """Measurement of detected diaSources as dipoles"""

    def setDefaults(self):
        SingleFrameMeasurementConfig.setDefaults(self)
        self.plugins = ["base_CircularApertureFlux",
                        "base_PixelFlags",
                        "base_SkyCoord",
                        "base_PsfFlux",
                        "ip_diffim_PsfDipoleFlux",
                        "ip_diffim_ClassificationDipole",
                        ]

        self.slots.calibFlux = None
        self.slots.modelFlux = None
        self.slots.gaussianFlux = None
        self.slots.shape = None
        self.slots.centroid = "ip_diffim_PsfDipoleFlux"
        self.doReplaceWithNoise = False


class DipoleMeasurementTask(SingleFrameMeasurementTask):
    """Measurement of Sources, specifically ones from difference images, for characterization as dipoles

    Parameters
    ----------
    sources : 'lsst.afw.table.SourceCatalog'
        Sources that will be measured
    badFlags : `list` of `dict`
        A list of flags that will be used to determine if there was a measurement problem

    """
    ConfigClass = DipoleMeasurementConfig
    _DefaultName = "dipoleMeasurement"


#########
# Other Support classs
#########

class SourceFlagChecker(object):
    """Functor class to check whether a diaSource has flags set that should cause it to be labeled bad."""

    def __init__(self, sources, badFlags=None):
        self.badFlags = ['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter',
                         'base_PixelFlags_flag_saturatedCenter']
        if badFlags is not None:
            for flag in badFlags:
                self.badFlags.append(flag)
        self.keys = [sources.getSchema().find(name).key for name in self.badFlags]
        self.keys.append(sources.table.getCentroidFlagKey())

    def __call__(self, source):
        """Call the source flag checker on a single Source

        Parameters
        ----------
        source :
            Source that will be examined
        """
        for k in self.keys:
            if source.get(k):
                return False
        return True


class DipoleAnalysis(object):
    """Functor class that provides (S/N, position, orientation) of measured dipoles"""

    def __init__(self):
        pass

    def __call__(self, source):
        """Parse information returned from dipole measurement

        Parameters
        ----------
        source : `lsst.afw.table.SourceRecord`
            The source that will be examined"""
        return self.getSn(source), self.getCentroid(source), self.getOrientation(source)

    def getSn(self, source):
        """Get the total signal-to-noise of the dipole; total S/N is from positive and negative lobe

        Parameters
        ----------
        source : `lsst.afw.table.SourceRecord`
            The source that will be examined"""

        posflux = source.get("ip_diffim_PsfDipoleFlux_pos_instFlux")
        posfluxErr = source.get("ip_diffim_PsfDipoleFlux_pos_instFluxErr")
        negflux = source.get("ip_diffim_PsfDipoleFlux_neg_instFlux")
        negfluxErr = source.get("ip_diffim_PsfDipoleFlux_neg_instFluxErr")

        # Not a dipole!
        if (posflux < 0) is (negflux < 0):
            return 0

        return np.sqrt((posflux/posfluxErr)**2 + (negflux/negfluxErr)**2)

    def getCentroid(self, source):
        """Get the centroid of the dipole; average of positive and negative lobe

        Parameters
        ----------
        source : `lsst.afw.table.SourceRecord`
            The source that will be examined"""

        negCenX = source.get("ip_diffim_PsfDipoleFlux_neg_centroid_x")
        negCenY = source.get("ip_diffim_PsfDipoleFlux_neg_centroid_y")
        posCenX = source.get("ip_diffim_PsfDipoleFlux_pos_centroid_x")
        posCenY = source.get("ip_diffim_PsfDipoleFlux_pos_centroid_y")
        if (np.isinf(negCenX) or np.isinf(negCenY) or np.isinf(posCenX) or np.isinf(posCenY)):
            return None

        center = geom.Point2D(0.5*(negCenX+posCenX),
                              0.5*(negCenY+posCenY))
        return center

    def getOrientation(self, source):
        """Calculate the orientation of dipole; vector from negative to positive lobe

        Parameters
        ----------
        source : `lsst.afw.table.SourceRecord`
            The source that will be examined"""

        negCenX = source.get("ip_diffim_PsfDipoleFlux_neg_centroid_x")
        negCenY = source.get("ip_diffim_PsfDipoleFlux_neg_centroid_y")
        posCenX = source.get("ip_diffim_PsfDipoleFlux_pos_centroid_x")
        posCenY = source.get("ip_diffim_PsfDipoleFlux_pos_centroid_y")
        if (np.isinf(negCenX) or np.isinf(negCenY) or np.isinf(posCenX) or np.isinf(posCenY)):
            return None

        dx, dy = posCenX-negCenX, posCenY-negCenY
        angle = geom.Angle(np.arctan2(dx, dy), geom.radians)
        return angle

    def displayDipoles(self, exposure, sources):
        """Display debugging information on the detected dipoles

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Image the dipoles were measured on
        sources : `lsst.afw.table.SourceCatalog`
            The set of diaSources that were measured"""

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayDiaSources = lsstDebug.Info(__name__).displayDiaSources
        maskTransparency = lsstDebug.Info(__name__).maskTransparency
        if not maskTransparency:
            maskTransparency = 90
        disp = afwDisplay.Display(frame=lsstDebug.frame)
        disp.setMaskTransparency(maskTransparency)
        disp.mtv(exposure)

        if display and displayDiaSources:
            with disp.Buffering():
                for source in sources:
                    cenX = source.get("ipdiffim_DipolePsfFlux_x")
                    cenY = source.get("ipdiffim_DipolePsfFlux_y")
                    if np.isinf(cenX) or np.isinf(cenY):
                        cenX, cenY = source.getCentroid()

                    isdipole = source.get("ip_diffim_ClassificationDipole_value")
                    if isdipole and np.isfinite(isdipole):
                        # Dipole
                        ctype = afwDisplay.GREEN
                    else:
                        # Not dipole
                        ctype = afwDisplay.RED

                    disp.dot("o", cenX, cenY, size=2, ctype=ctype)

                    negCenX = source.get("ip_diffim_PsfDipoleFlux_neg_centroid_x")
                    negCenY = source.get("ip_diffim_PsfDipoleFlux_neg_centroid_y")
                    posCenX = source.get("ip_diffim_PsfDipoleFlux_pos_centroid_x")
                    posCenY = source.get("ip_diffim_PsfDipoleFlux_pos_centroid_y")
                    if (np.isinf(negCenX) or np.isinf(negCenY) or np.isinf(posCenX) or np.isinf(posCenY)):
                        continue

                    disp.line([(negCenX, negCenY), (posCenX, posCenY)], ctype=afwDisplay.YELLOW)

            lsstDebug.frame += 1


class DipoleDeblender(object):
    """Functor to deblend a source as a dipole, and return a new source with deblended footprints.

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
        self.log = getLogger('lsst.ip.diffim.DipoleDeblender')
        self.sigma2fwhm = 2. * np.sqrt(2. * np.log(2.))

    def __call__(self, source, exposure):
        fp = source.getFootprint()
        peaks = fp.getPeaks()
        peaksF = [pk.getF() for pk in peaks]
        fbb = fp.getBBox()
        fmask = afwImage.Mask(fbb)
        fmask.setXY0(fbb.getMinX(), fbb.getMinY())
        fp.spans.setMask(fmask, 1)

        psf = exposure.getPsf()
        psfSigPix = psf.computeShape(psf.getAveragePosition()).getDeterminantRadius()
        psfFwhmPix = psfSigPix * self.sigma2fwhm
        subimage = afwImage.ExposureF(exposure, bbox=fbb, deep=True)
        cpsf = deblendBaseline.CachingPsf(psf)

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
                                            log=self.log,
                                            psfChisqCut1=self.psfChisqCut1,
                                            psfChisqCut2=self.psfChisqCut2,
                                            psfChisqCut2b=self.psfChisqCut2b)
        else:
            # Call lower-level _fit_psf task

            # Prepare results structure
            fpres = deblendBaseline.DeblenderResult(fp, exposure.getMaskedImage(), psf, psfFwhmPix, self.log)

            for pki, (pk, pkres, pkF) in enumerate(zip(dpeaks, fpres.deblendedParents[0].peaks, peaksF)):
                self.log.debug('Peak %i', pki)
                deblendBaseline._fitPsf(fp, fmask, pk, pkF, pkres, fbb, dpeaks, peaksF, self.log,
                                        cpsf, psfFwhmPix,
                                        subimage.image,
                                        subimage.variance,
                                        self.psfChisqCut1, self.psfChisqCut2, self.psfChisqCut2b)

        deblendedSource = source.getTable().copyRecord(source)
        deblendedSource.setParent(source.getId())
        peakList = deblendedSource.getFootprint().getPeaks()
        peakList.clear()

        for i, peak in enumerate(fpres.deblendedParents[0].peaks):
            if peak.psfFitFlux > 0:
                suffix = "pos"
            else:
                suffix = "neg"
            c = peak.psfFitCenter
            self.log.info("deblended.centroid.dipole.psf.%s %f %f",
                          suffix, c[0], c[1])
            self.log.info("deblended.chi2dof.dipole.%s %f",
                          suffix, peak.psfFitChisq / peak.psfFitDof)
            self.log.info("deblended.flux.dipole.psf.%s %f",
                          suffix, peak.psfFitFlux * np.sum(peak.templateImage.array))
            peakList.append(peak.peak)
        return deblendedSource
