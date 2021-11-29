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

__all__ = ["SnapPsfMatchConfigDF", "SnapPsfMatchConfigAL", "SnapPsfMatchConfig", "SnapPsfMatchTask"]

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


class SnapPsfMatchTask(ImagePsfMatchTask):
    """Image-based Psf-matching of two subsequent snaps from the same visit

    Notes
    -----
    This Task differs from ImagePsfMatchTask in that it matches two Exposures assuming that the images have
    been acquired very closely in time.  Under this assumption, the astrometric misalignments and/or
    relative distortions should be within a pixel, and the Psf-shapes should be very similar.  As a
    consequence, the default configurations for this class assume a very simple solution.

    - The spatial variation in the kernel (SnapPsfMatchConfig.spatialKernelOrder) is assumed to be zero

    - With no spatial variation, we turn of the spatial
        clipping loops (SnapPsfMatchConfig.spatialKernelClipping)

    - The differential background is not fit for (SnapPsfMatchConfig.fitForBackground)

    - The kernel is expected to be appx.
        a delta function, and has a small size (SnapPsfMatchConfig.kernelSize)

    The sub-configurations for the Alard-Lupton (SnapPsfMatchConfigAL)
    and delta-function (SnapPsfMatchConfigDF)
    bases also are designed to generate a small, simple kernel.

    Task initialization

    Initialization is the same as base class ImagePsfMatch.__init__,
    with the difference being that the Task's
    ConfigClass is SnapPsfMatchConfig.

    Invoking the Task

    The Task is only configured to have a subtractExposures method, which in turn calls
    ImagePsfMatchTask.subtractExposures.

    Configuration parameters

    See SnapPsfMatchConfig, which uses either SnapPsfMatchConfigDF and SnapPsfMatchConfigAL
    as its active configuration.

    Debug variables

    The lsst.pipe.base.cmdLineTask.CmdLineTask command line task interface supports a
    flag -d/--debug to importdebug.py from your PYTHONPATH.  The relevant contents of debug.py
    for this Task include:

    .. code-block:: py

        import sys
        import lsstDebug
        def DebugInfo(name):
            di = lsstDebug.getInfo(name)
            if name == "lsst.ip.diffim.psfMatch":
                di.display = True                 # enable debug output
                di.maskTransparency = 80          # display mask transparency
                di.displayCandidates = True       # show all the candidates and residuals
                di.displayKernelBasis = False     # show kernel basis functions
                di.displayKernelMosaic = True     # show kernel realized across the image
                di.plotKernelSpatialModel = False # show coefficients of spatial model
                di.showBadCandidates = True       # show the bad candidates (red) along with good (green)
            elif name == "lsst.ip.diffim.imagePsfMatch":
                di.display = True                 # enable debug output
                di.maskTransparency = 30          # display mask transparency
                di.displayTemplate = True         # show full (remapped) template
                di.displaySciIm = True            # show science image to match to
                di.displaySpatialCells = True     # show spatial cells
                di.displayDiffIm = True           # show difference image
                di.showBadCandidates = True       # show the bad candidates (red) along with good (green)
            elif name == "lsst.ip.diffim.diaCatalogSourceSelector":
                di.display = False                # enable debug output
                di.maskTransparency = 30          # display mask transparency
                di.displayExposure = True         # show exposure with candidates indicated
                di.pauseAtEnd = False             # pause when done
            return di
        lsstDebug.Info = DebugInfo
        lsstDebug.frame = 1

    Note that if you want addional logging info, you may add to your scripts:

    .. code-block:: py

        import lsst.utils.logging as logUtils
        logUtils.trace_set_at("lsst.ip.diffim", 4)

    Examples
    --------
    This code is snapPsfMatchTask.py in the examples directory, and can be run as e.g.

    .. code-block:: py

        examples/snapPsfMatchTask.py
        examples/snapPsfMatchTask.py --debug
        examples/snapPsfMatchTask.py --debug --template /path/to/templateExp.fits
        --science /path/to/scienceExp.fits

    First, create a subclass of SnapPsfMatchTask that accepts two exposures.
    Ideally these exposures would have been taken back-to-back,
    such that the pointing/background/Psf does not vary substantially between the two:

    .. code-block:: py

        class MySnapPsfMatchTask(SnapPsfMatchTask):
            def __init__(self, *args, **kwargs):
                SnapPsfMatchTask.__init__(self, *args, **kwargs)
            def run(self, templateExp, scienceExp):
                return self.subtractExposures(templateExp, scienceExp)

    And allow the user the freedom to either run the script in default mode,
    or point to their own images on disk. Note that these images must be
    readable as an lsst.afw.image.Exposure

    .. code-block:: py

        if __name__ == "__main__":
            import argparse
            parser = argparse.ArgumentParser(description="Demonstrate the use of ImagePsfMatchTask")
            parser.add_argument("--debug", "-d", action="store_true", help="Load debug.py?", default=False)
            parser.add_argument("--template", "-t", help="Template Exposure to use", default=None)
            parser.add_argument("--science", "-s", help="Science Exposure to use", default=None)
            args = parser.parse_args()

    We have enabled some minor display debugging in this script via the –debug option. However,
    if you have an lsstDebug debug.in your PYTHONPATH you will get additional debugging displays.
    The following block checks for this script

    .. code-block:: py

        if args.debug:
            try:
                import debug
                # Since I am displaying 2 images here, set the starting frame number for the LSST debug LSST
                debug.lsstDebug.frame = 3
            except ImportError as e:
                print(e, file=sys.stderr)

    Finally, we call a run method that we define below.
    First set up a Config and choose the basis set to use:

    .. code-block:: py

        def run(args):
            #
            # Create the Config and use sum of gaussian basis
            #
            config = SnapPsfMatchTask.ConfigClass()
            config.doWarping = True
            config.kernel.name = "AL"

    Make sure the images (if any) that were sent to the script exist on disk and are readable.
    If no images are sent, make some fake data up for the sake of this example script
    (have a look at the code if you want more details on generateFakeImages;
    as a detail of how the fake images were made, you do have to fit for a differential background):

    .. code-block:: py

        # Run the requested method of the Task
        if args.template is not None and args.science is not None:
            if not os.path.isfile(args.template):
                raise FileNotFoundError("Template image %s does not exist" % (args.template))
            if not os.path.isfile(args.science):
                raise FileNotFoundError("Science image %s does not exist" % (args.science))
            try:
                templateExp = afwImage.ExposureF(args.template)
            except Exception as e:
                raise RuntimeError("Cannot read template image %s" % (args.template))
            try:
                scienceExp = afwImage.ExposureF(args.science)
            except Exception as e:
                raise RuntimeError("Cannot read science image %s" % (args.science))
        else:
            templateExp, scienceExp = generateFakeImages()
            config.kernel.active.fitForBackground = True
            config.kernel.active.spatialBgOrder = 0
            config.kernel.active.sizeCellX = 128
            config.kernel.active.sizeCellY = 128

    Display the two images if -debug

    .. code-block:: py

        if args.debug:
            afwDisplay.Display(frame=1).mtv(templateExp, title="Example script: Input Template")
            afwDisplay.Display(frame=2).mtv(scienceExp, title="Example script: Input Science Image")

    Create and run the Task

    .. code-block:: py

        # Create the Task
        psfMatchTask = MySnapPsfMatchTask(config=config)
        # Run the Task
        result = psfMatchTask.run(templateExp, scienceExp)

    And finally provide optional debugging display of the Psf-matched (via the Psf models) science image:

    .. code-block:: py

        if args.debug:
            # See if the LSST debug has incremented the frame number; if not start with frame 3
            try:
                frame = debug.lsstDebug.frame + 1
            except Exception:
                frame = 3
            afwDisplay.Display(frame=frame).mtv(result.matchedExposure,
                                            title="Example script: Matched Template Image")
            if "subtractedExposure" in result.getDict():
                afwDisplay.Display(frame=frame + 1).mtv(result.subtractedExposure,
                                                        title="Example script: Subtracted Image")

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
