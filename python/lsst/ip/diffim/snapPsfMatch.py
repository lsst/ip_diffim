from __future__ import absolute_import
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
import lsst.pex.config as pexConfig
from .psfMatch import PsfMatchConfigDF, PsfMatchConfigAL
from .imagePsfMatch import ImagePsfMatchTask, ImagePsfMatchConfig


class SnapPsfMatchConfigDF(PsfMatchConfigDF):
    """Delta-function Psf-matching config optimized for snap subtraction"""

    def setDefaults(self):
        PsfMatchConfigDF.setDefaults(self)

        # No regularization
        self.useRegularization = False

        # Pca
        self.usePcaForSpatialKernel = True
        self.subtractMeanForPca = True
        self.numPrincipalComponents = 5


class SnapPsfMatchConfigAL(PsfMatchConfigAL):
    """Sum-of-Gaussian (Alard-Lupton) Psf-matching config optimized for snap subtraction"""

    def setDefaults(self):
        PsfMatchConfigAL.setDefaults(self)

        # Simple basis set
        self.alardNGauss = 2
        self.alardDegGauss = (4, 2)
        self.alardSigGauss = (1.0, 2.5)


class SnapPsfMatchConfig(ImagePsfMatchConfig):
    kernel = pexConfig.ConfigChoiceField(
        doc="kernel type",
        typemap=dict(
            AL=SnapPsfMatchConfigAL,
            DF=SnapPsfMatchConfigDF
        ),
        default="AL",
    )

    doWarping = pexConfig.Field(
        dtype=bool,
        doc="Warp the snaps?",
        default=False
    )

    def setDefaults(self):
        ImagePsfMatchConfig.setDefaults(self)

        # No spatial variation in model
        self.kernel.active.spatialKernelOrder = 0

        # Don't fit for differential background
        self.kernel.active.fitForBackground = False

        # Small kernel size
        self.kernel.active.kernelSize = 7

        # With zero spatial order don't worry about spatial clipping
        self.kernel.active.spatialKernelClipping = False

## \addtogroup LSST_task_documentation
## \{
## \page SnapPsfMatchTask
## \ref SnapPsfMatchTask_ "SnapPsfMatchTask"
## \copybrief SnapPsfMatchTask
## \}


class SnapPsfMatchTask(ImagePsfMatchTask):
    """!
\anchor SnapPsfMatchTask_

\brief Image-based Psf-matching of two subsequent snaps from the same visit

\section ip_diffim_snappsfmatch_Contents Contents

 - \ref ip_diffim_snappsfmatch_Purpose
 - \ref ip_diffim_snappsfmatch_Initialize
 - \ref ip_diffim_snappsfmatch_IO
 - \ref ip_diffim_snappsfmatch_Config
 - \ref ip_diffim_snappsfmatch_Metadata
 - \ref ip_diffim_snappsfmatch_Debug
 - \ref ip_diffim_snappsfmatch_Example

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_snappsfmatch_Purpose   Description

\copybrief SnapPsfMatchTask

This Task differs from ImagePsfMatchTask in that it matches two Exposures assuming that the images have
been acquired very closely in time.  Under this assumption, the astrometric misalignments and/or
relative distortions should be within a pixel, and the Psf-shapes should be very similar.  As a
consequence, the default configurations for this class assume a very simple solution.

 . The spatial variation in the kernel (SnapPsfMatchConfig.spatialKernelOrder) is assumed to be zero

 . With no spatial variation, we turn of the spatial clipping loops (SnapPsfMatchConfig.spatialKernelClipping)

 . The differential background is _not_ fit for (SnapPsfMatchConfig.fitForBackground)

 . The kernel is expected to be appx. a delta function, and has a small size (SnapPsfMatchConfig.kernelSize)

The sub-configurations for the Alard-Lupton (SnapPsfMatchConfigAL) and delta-function (SnapPsfMatchConfigDF)
bases also are designed to generate a small, simple kernel.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_snappsfmatch_Initialize    Task initialization

Initialization is the same as base class ImagePsfMatch.__init__, with the difference being that the Task's
ConfigClass is SnapPsfMatchConfig.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_snappsfmatch_IO        Invoking the Task

The Task is only configured to have a subtractExposures method, which in turn calls
ImagePsfMatchTask.subtractExposures.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_snappsfmatch_Config       Configuration parameters

See \ref SnapPsfMatchConfig, which uses either \ref SnapPsfMatchConfigDF and \ref SnapPsfMatchConfigAL
as its active configuration.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_snappsfmatch_Metadata   Quantities set in Metadata

See \ref ip_diffim_psfmatch_Metadata "PsfMatchTask"

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_snappsfmatch_Debug     Debug variables

The \link lsst.pipe.base.cmdLineTask.CmdLineTask command line task\endlink interface supports a
flag \c -d/--debug to import \b debug.py from your \c PYTHONPATH.  The relevant contents of debug.py
for this Task include:

\code{.py}
    import sys
    import lsstDebug
    def DebugInfo(name):
        di = lsstDebug.getInfo(name)
        if name == "lsst.ip.diffim.psfMatch":
            di.display = True                 # enable debug output
            di.maskTransparency = 80          # ds9 mask transparency
            di.displayCandidates = True       # show all the candidates and residuals
            di.displayKernelBasis = False     # show kernel basis functions
            di.displayKernelMosaic = True     # show kernel realized across the image
            di.plotKernelSpatialModel = False # show coefficients of spatial model
            di.showBadCandidates = True       # show the bad candidates (red) along with good (green)
        elif name == "lsst.ip.diffim.imagePsfMatch":
            di.display = True                 # enable debug output
            di.maskTransparency = 30          # ds9 mask transparency
            di.displayTemplate = True         # show full (remapped) template
            di.displaySciIm = True            # show science image to match to
            di.displaySpatialCells = True     # show spatial cells
            di.displayDiffIm = True           # show difference image
            di.showBadCandidates = True       # show the bad candidates (red) along with good (green)
        elif name == "lsst.ip.diffim.diaCatalogSourceSelector":
            di.display = False                # enable debug output
            di.maskTransparency = 30          # ds9 mask transparency
            di.displayExposure = True         # show exposure with candidates indicated
            di.pauseAtEnd = False             # pause when done
        return di
    lsstDebug.Info = DebugInfo
    lsstDebug.frame = 1
\endcode

Note that if you want addional logging info, you may add to your scripts:
\code{.py}
import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)
\endcode

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_snappsfmatch_Example   A complete example of using SnapPsfMatchTask

This code is snapPsfMatchTask.py in the examples directory, and can be run as \em e.g.
\code
examples/snapPsfMatchTask.py
examples/snapPsfMatchTask.py --debug
examples/snapPsfMatchTask.py --debug --template /path/to/templateExp.fits --science /path/to/scienceExp.fits
\endcode

\dontinclude snapPsfMatchTask.py
First, create a subclass of SnapPsfMatchTask that accepts two exposures.  Ideally these exposures would have
been taken back-to-back, such that the pointing/background/Psf does not vary substantially between the two:
\skip MySnapPsfMatchTask
\until return

And allow the user the freedom to either run the script in default mode, or point to their own images on disk.
Note that these images must be readable as an lsst.afw.image.Exposure:
\skip main
\until parse_args

We have enabled some minor display debugging in this script via the --debug option.  However, if you
have an lsstDebug debug.py in your PYTHONPATH you will get additional debugging displays.  The following
block checks for this script:
\skip args.debug
\until sys.stderr


\dontinclude snapPsfMatchTask.py
Finally, we call a run method that we define below.  First set up a Config and choose the basis set to use:
\skip run(args)
\until AL

Make sure the images (if any) that were sent to the script exist on disk and are readable.  If no images
are sent, make some fake data up for the sake of this example script (have a look at the code if you want
more details on generateFakeImages; as a detail of how the fake images were made, you do have to fit for a
differential background):
\skip requested
\until sizeCellY

Display the two images if --debug:
\skip args.debug
\until Science

Create and run the Task:
\skip Create
\until result

And finally provide optional debugging display of the Psf-matched (via the Psf models) science image:
\skip args.debug
\until result.subtractedExposure

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    """
    ConfigClass = SnapPsfMatchConfig

    # Override ImagePsfMatchTask.subtractExposures to set doWarping on config.doWarping
    def subtractExposures(self, templateExposure, scienceExposure,
                          templateFwhmPix=None, scienceFwhmPix=None,
                          candidateList=None):
        return ImagePsfMatchTask.subtractExposures(self,
                                                   templateExposure=templateExposure,
                                                   scienceExposure=scienceExposure,
                                                   templateFwhmPix=templateFwhmPix,
                                                   scienceFwhmPix=scienceFwhmPix,
                                                   candidateList=candidateList,
                                                   doWarping=self.config.doWarping,
                                                   )
