#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
from __future__ import absolute_import, division, print_function

from builtins import zip
from builtins import object
import numpy as np

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetect
import lsst.pex.config as pexConfig
from lsst.log import Log
import lsst.meas.deblender.baseline as deblendBaseline
from lsst.meas.base.pluginRegistry import register
from lsst.meas.base import SingleFrameMeasurementTask, SingleFrameMeasurementConfig, \
    SingleFramePluginConfig, SingleFramePlugin
import lsst.afw.display.ds9 as ds9

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
        return cls.APCORR_ORDER

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


class DipoleMeasurementConfig(SingleFrameMeasurementConfig):
    """!Measurement of detected diaSources as dipoles"""

    def setDefaults(self):
        SingleFrameMeasurementConfig.setDefaults(self)
        self.plugins = ["base_CircularApertureFlux",
                        "base_PixelFlags",
                        "base_SkyCoord",
                        "base_PsfFlux",
                        "ip_diffim_NaiveDipoleCentroid",
                        "ip_diffim_NaiveDipoleFlux",
                        "ip_diffim_PsfDipoleFlux",
                        "ip_diffim_ClassificationDipole",
                        ]

        self.slots.calibFlux = None
        self.slots.modelFlux = None
        self.slots.instFlux = None
        self.slots.shape = None
        self.slots.centroid = "ip_diffim_NaiveDipoleCentroid"
        self.doReplaceWithNoise = False

## \addtogroup LSST_task_documentation
## \{
## \page DipoleMeasurementTask
## \ref DipoleMeasurementTask_ "DipoleMeasurementTask"
## \copybrief DipoleMeasurementTask
## \}


class DipoleMeasurementTask(SingleFrameMeasurementTask):
    """!
\anchor DipoleMeasurementTask_

\brief Measurement of Sources, specifically ones from difference images, for characterization as dipoles

\section ip_diffim_dipolemeas_Contents Contents

 - \ref ip_diffim_dipolemeas_Purpose
 - \ref ip_diffim_dipolemeas_Initialize
 - \ref ip_diffim_dipolemeas_IO
 - \ref ip_diffim_dipolemeas_Config
 - \ref ip_diffim_dipolemeas_Metadata
 - \ref ip_diffim_dipolemeas_Debug
 - \ref ip_diffim_dipolemeas_Example

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_dipolemeas_Purpose   Description

This class provides a default configuration for running Source measurement on image differences.

These default plugins include:
\dontinclude dipoleMeasurement.py
\skip class DipoleMeasurementConfig
@until self.doReplaceWithNoise

These plugins enabled by default allow the user to test the hypothesis that the Source is a dipole.
This includes a set of measurements derived from intermediate base classes
DipoleCentroidAlgorithm and DipoleFluxAlgorithm.  Their respective algorithm control classes are defined in
DipoleCentroidControl and DipoleFluxControl.  Each centroid and flux measurement will have _neg (negative)
and _pos (positive lobe) fields.

The first set of measurements uses a "naive" alrogithm for centroid and flux measurements, implemented in
NaiveDipoleCentroidControl and NaiveDipoleFluxControl.  The algorithm uses a naive 3x3 weighted moment around
the nominal centroids of each peak in the Source Footprint.  These algorithms fill the table fields
ip_diffim_NaiveDipoleCentroid* and ip_diffim_NaiveDipoleFlux*

The second set of measurements undertakes a joint-Psf model on the negative and positive lobe simultaneously.
This fit simultaneously solves for the negative and positive lobe centroids and fluxes using non-linear
least squares minimization.  The fields are stored in table elements ip_diffim_PsfDipoleFlux*.

Because this Task is just a config for SourceMeasurementTask, the same result may be acheived by manually
editing the config and running SourceMeasurementTask. For example:

\code
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

schema = afwTable.SourceTable.makeMinimalSchema()
task = SingleFrameMeasurementTask(schema, config=config)

task.run(sources, exposure)
\endcode


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_dipolemeas_Initialize    Task initialization

\copydoc \_\_init\_\_

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_dipolemeas_IO        Invoking the Task

\copydoc run

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_dipolemeas_Config       Configuration parameters

See \ref DipoleMeasurementConfig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_dipolemeas_Metadata   Quantities set in Metadata

No specific values are set in the Task metadata.  However, the Source schema are modified to store the
results of the dipole-specific measurements.


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_dipolemeas_Debug     Debug variables

The \link lsst.pipe.base.cmdLineTask.CmdLineTask command line task\endlink interface supports a
flag \c -d/--debug to import \b debug.py from your \c PYTHONPATH.  The relevant contents of debug.py
for this Task include:

\code{.py}
    import sys
    import lsstDebug
    def DebugInfo(name):
        di = lsstDebug.getInfo(name)
        if name == "lsst.ip.diffim.dipoleMeasurement":
            di.display = True                 # enable debug output
            di.maskTransparency = 90          # ds9 mask transparency
            di.displayDiaSources = True       # show exposure with dipole results
        return di
    lsstDebug.Info = DebugInfo
    lsstDebug.frame = 1
\endcode

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_dipolemeas_Example   A complete example of using DipoleMeasurementTask

This code is dipoleMeasTask.py in the examples directory, and can be run as \em e.g.
\code
examples/dipoleMeasTask.py
examples/dipoleMeasTask.py --debug
examples/dipoleMeasTask.py --debug --image /path/to/image.fits
\endcode

\dontinclude dipoleMeasTask.py
Start the processing by parsing the command line, where the user has the option of enabling debugging output
and/or sending their own image for demonstration (in case they have not downloaded the afwdata package).
\skip main
@until run

\dontinclude dipoleMeasTask.py
The processing occurs in the run function.  We first extract an exposure from disk or afwdata, displaying
it if requested:
\skip args
@until mtv

Create a default source schema that we will append fields to as we add more algorithms:
\skip makeMinimalSchema
@until makeMinimalSchema

Create the detection and measurement Tasks, with some minor tweaking of their configs:
\skip Create
@until measureTask

Having fully initialied the schema, we create a Source table from it:
\skip output
@until SourceTable

Run detection:
\skip Process
@until detectionTask

Because we are looking for dipoles, we need to merge the positive and negative detections:
\skip Merge
@until numNeg

Finally, perform measurement (both standard and dipole-specialized) on the merged sources:
\skip measureTask
@until measureTask

Optionally display debugging information:
\skip Display
@until displayDipoles
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    """
    ConfigClass = DipoleMeasurementConfig
    _DefaultName = "dipoleMeasurement"


#########
# Other Support classs
#########

class SourceFlagChecker(object):
    """!Functor class to check whether a diaSource has flags set that should cause it to be labeled bad."""

    def __init__(self, sources, badFlags=None):
        """!Constructor

        @param sources     Sources that will be measured
        @param badFlags    A list of flags that will be used to determine if there was a measurement problem

        The list of badFlags will be used to make a list of keys to check for measurement flags on.  By
        default the centroid keys are added to this list"""

        self.badFlags = ['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter',
                         'base_PixelFlags_flag_saturatedCenter']
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
        angle = afwGeom.Angle(np.arctan2(dx, dy), afwGeom.radians)
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
                        ctype = "green"
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
        self.log = Log.getLogger('ip.diffim.DipoleDeblender')
        self.sigma2fwhm = 2. * np.sqrt(2. * np.log(2.))

    def __call__(self, source, exposure):
        fp = source.getFootprint()
        peaks = fp.getPeaks()
        peaksF = [pk.getF() for pk in peaks]
        fbb = fp.getBBox()
        fmask = afwImage.MaskU(fbb)
        fmask.setXY0(fbb.getMinX(), fbb.getMinY())
        fp.spans.setMask(fmask, 1)

        psf = exposure.getPsf()
        psfSigPix = psf.computeShape().getDeterminantRadius()
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
                                        subimage.getMaskedImage().getImage(),
                                        subimage.getMaskedImage().getVariance(),
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
                          suffix, peak.psfFitFlux * np.sum(peak.templateImage.getArray()))
            peakList.append(peak.peak)
        return deblendedSource
