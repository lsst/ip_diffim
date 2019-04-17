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

__all__ = ["DetectionConfig", "PsfMatchConfig", "PsfMatchConfigAL", "PsfMatchConfigDF", "PsfMatchTask"]

import time

import numpy as np

import lsst.afw.image as afwImage
import lsst.pex.config as pexConfig
import lsst.afw.math as afwMath
import lsst.afw.display as afwDisplay
import lsst.log as log
import lsst.pipe.base as pipeBase
from lsst.meas.algorithms import SubtractBackgroundConfig
from . import utils as diutils
from . import diffimLib


class DetectionConfig(pexConfig.Config):
    """Configuration for detecting sources on images for building a
    PSF-matching kernel

    Configuration for turning detected lsst.afw.detection.FootPrints into an
    acceptable (unmasked, high signal-to-noise, not too large or not too small)
    list of `lsst.ip.diffim.KernelSources` that are used to build the
    Psf-matching kernel"""

    detThreshold = pexConfig.Field(
        dtype=float,
        doc="Value of footprint detection threshold",
        default=10.0,
        check=lambda x: x >= 3.0
    )
    detThresholdType = pexConfig.ChoiceField(
        dtype=str,
        doc="Type of detection threshold",
        default="pixel_stdev",
        allowed={
            "value": "Use counts as the detection threshold type",
            "stdev": "Use standard deviation of image plane",
            "variance": "Use variance of image plane",
            "pixel_stdev": "Use stdev derived from variance plane"
        }
    )
    detOnTemplate = pexConfig.Field(
        dtype=bool,
        doc="""If true run detection on the template (image to convolve);
                 if false run detection on the science image""",
        default=True
    )
    badMaskPlanes = pexConfig.ListField(
        dtype=str,
        doc="""Mask planes that lead to an invalid detection.
                 Options: NO_DATA EDGE SAT BAD CR INTRP""",
        default=("NO_DATA", "EDGE", "SAT")
    )
    fpNpixMin = pexConfig.Field(
        dtype=int,
        doc="Minimum number of pixels in an acceptable Footprint",
        default=5,
        check=lambda x: x >= 5
    )
    fpNpixMax = pexConfig.Field(
        dtype=int,
        doc="""Maximum number of pixels in an acceptable Footprint;
                 too big and the subsequent convolutions become unwieldy""",
        default=500,
        check=lambda x: x <= 500
    )
    fpGrowKernelScaling = pexConfig.Field(
        dtype=float,
        doc="""If config.scaleByFwhm, grow the footprint based on
                 the final kernelSize.  Each footprint will be
                 2*fpGrowKernelScaling*kernelSize x
                 2*fpGrowKernelScaling*kernelSize.  With the value
                 of 1.0, the remaining pixels in each KernelCandiate
                 after convolution by the basis functions will be
                 equal to the kernel size itself.""",
        default=1.0,
        check=lambda x: x >= 1.0
    )
    fpGrowPix = pexConfig.Field(
        dtype=int,
        doc="""Growing radius (in pixels) for each raw detection
                 footprint.  The smaller the faster; however the
                 kernel sum does not converge if the stamp is too
                 small; and the kernel is not constrained at all if
                 the stamp is the size of the kernel.  The grown stamp
                 is 2 * fpGrowPix pixels larger in each dimension.
                 This is overridden by fpGrowKernelScaling if scaleByFwhm""",
        default=30,
        check=lambda x: x >= 10
    )
    scaleByFwhm = pexConfig.Field(
        dtype=bool,
        doc="Scale fpGrowPix by input Fwhm?",
        default=True,
    )


class PsfMatchConfig(pexConfig.Config):
    """Base configuration for Psf-matching

    The base configuration of the Psf-matching kernel, and of the warping, detection,
    and background modeling subTasks."""

    warpingConfig = pexConfig.ConfigField("Config for warping exposures to a common alignment",
                                          afwMath.warper.WarperConfig)
    detectionConfig = pexConfig.ConfigField("Controlling the detection of sources for kernel building",
                                            DetectionConfig)
    afwBackgroundConfig = pexConfig.ConfigField("Controlling the Afw background fitting",
                                                SubtractBackgroundConfig)

    useAfwBackground = pexConfig.Field(
        dtype=bool,
        doc="Use afw background subtraction instead of ip_diffim",
        default=False,
    )
    fitForBackground = pexConfig.Field(
        dtype=bool,
        doc="Include terms (including kernel cross terms) for background in ip_diffim",
        default=False,
    )
    kernelBasisSet = pexConfig.ChoiceField(
        dtype=str,
        doc="Type of basis set for PSF matching kernel.",
        default="alard-lupton",
        allowed={
            "alard-lupton": """Alard-Lupton sum-of-gaussians basis set,
                           * The first term has no spatial variation
                           * The kernel sum is conserved
                           * You may want to turn off 'usePcaForSpatialKernel'""",
            "delta-function": """Delta-function kernel basis set,
                           * You may enable the option useRegularization
                           * You should seriously consider usePcaForSpatialKernel, which will also
                             enable kernel sum conservation for the delta function kernels"""
        }
    )
    kernelSize = pexConfig.Field(
        dtype=int,
        doc="""Number of rows/columns in the convolution kernel; should be odd-valued.
                 Modified by kernelSizeFwhmScaling if scaleByFwhm = true""",
        default=21,
    )
    scaleByFwhm = pexConfig.Field(
        dtype=bool,
        doc="Scale kernelSize, alardGaussians by input Fwhm",
        default=True,
    )
    kernelSizeFwhmScaling = pexConfig.Field(
        dtype=float,
        doc="Multiplier of the largest AL Gaussian basis sigma to get the kernel bbox (pixel) size.",
        default=6.0,
        check=lambda x: x >= 1.0
    )
    kernelSizeMin = pexConfig.Field(
        dtype=int,
        doc="Minimum kernel bbox (pixel) size.",
        default=21,
    )
    kernelSizeMax = pexConfig.Field(
        dtype=int,
        doc="Maximum kernel bbox (pixel) size.",
        default=35,
    )
    spatialModelType = pexConfig.ChoiceField(
        dtype=str,
        doc="Type of spatial functions for kernel and background",
        default="chebyshev1",
        allowed={
            "chebyshev1": "Chebyshev polynomial of the first kind",
            "polynomial": "Standard x,y polynomial",
        }
    )
    spatialKernelOrder = pexConfig.Field(
        dtype=int,
        doc="Spatial order of convolution kernel variation",
        default=2,
        check=lambda x: x >= 0
    )
    spatialBgOrder = pexConfig.Field(
        dtype=int,
        doc="Spatial order of differential background variation",
        default=1,
        check=lambda x: x >= 0
    )
    sizeCellX = pexConfig.Field(
        dtype=int,
        doc="Size (rows) in pixels of each SpatialCell for spatial modeling",
        default=128,
        check=lambda x: x >= 32
    )
    sizeCellY = pexConfig.Field(
        dtype=int,
        doc="Size (columns) in pixels of each SpatialCell for spatial modeling",
        default=128,
        check=lambda x: x >= 32
    )
    nStarPerCell = pexConfig.Field(
        dtype=int,
        doc="Number of KernelCandidates in each SpatialCell to use in the spatial fitting",
        default=3,
        check=lambda x: x >= 1
    )
    maxSpatialIterations = pexConfig.Field(
        dtype=int,
        doc="Maximum number of iterations for rejecting bad KernelCandidates in spatial fitting",
        default=3,
        check=lambda x: x >= 1 and x <= 5
    )
    usePcaForSpatialKernel = pexConfig.Field(
        dtype=bool,
        doc="""Use Pca to reduce the dimensionality of the kernel basis sets.
                 This is particularly useful for delta-function kernels.
                 Functionally, after all Cells have their raw kernels determined, we run
                 a Pca on these Kernels, re-fit the Cells using the eigenKernels and then
                 fit those for spatial variation using the same technique as for Alard-Lupton kernels.
                 If this option is used, the first term will have no spatial variation and the
                 kernel sum will be conserved.""",
        default=False,
    )
    subtractMeanForPca = pexConfig.Field(
        dtype=bool,
        doc="Subtract off the mean feature before doing the Pca",
        default=True,
    )
    numPrincipalComponents = pexConfig.Field(
        dtype=int,
        doc="""Number of principal components to use for Pca basis, including the
                 mean kernel if requested.""",
        default=5,
        check=lambda x: x >= 3
    )
    singleKernelClipping = pexConfig.Field(
        dtype=bool,
        doc="Do sigma clipping on each raw kernel candidate",
        default=True,
    )
    kernelSumClipping = pexConfig.Field(
        dtype=bool,
        doc="Do sigma clipping on the ensemble of kernel sums",
        default=True,
    )
    spatialKernelClipping = pexConfig.Field(
        dtype=bool,
        doc="Do sigma clipping after building the spatial model",
        default=True,
    )
    checkConditionNumber = pexConfig.Field(
        dtype=bool,
        doc="""Test for maximum condition number when inverting a kernel matrix.
                 Anything above maxConditionNumber is not used and the candidate is set as BAD.
                 Also used to truncate inverse matrix in estimateBiasedRisk.  However,
                 if you are doing any deconvolution you will want to turn this off, or use
                 a large maxConditionNumber""",
        default=False,
    )
    badMaskPlanes = pexConfig.ListField(
        dtype=str,
        doc="""Mask planes to ignore when calculating diffim statistics
                 Options: NO_DATA EDGE SAT BAD CR INTRP""",
        default=("NO_DATA", "EDGE", "SAT")
    )
    candidateResidualMeanMax = pexConfig.Field(
        dtype=float,
        doc="""Rejects KernelCandidates yielding bad difference image quality.
                 Used by BuildSingleKernelVisitor, AssessSpatialKernelVisitor.
                 Represents average over pixels of (image/sqrt(variance)).""",
        default=0.25,
        check=lambda x: x >= 0.0
    )
    candidateResidualStdMax = pexConfig.Field(
        dtype=float,
        doc="""Rejects KernelCandidates yielding bad difference image quality.
                 Used by BuildSingleKernelVisitor, AssessSpatialKernelVisitor.
                 Represents stddev over pixels of (image/sqrt(variance)).""",
        default=1.50,
        check=lambda x: x >= 0.0
    )
    useCoreStats = pexConfig.Field(
        dtype=bool,
        doc="""Use the core of the footprint for the quality statistics, instead of the entire footprint.
                 WARNING: if there is deconvolution we probably will need to turn this off""",
        default=False,
    )
    candidateCoreRadius = pexConfig.Field(
        dtype=int,
        doc="""Radius for calculation of stats in 'core' of KernelCandidate diffim.
                 Total number of pixels used will be (2*radius)**2.
                 This is used both for 'core' diffim quality as well as ranking of
                 KernelCandidates by their total flux in this core""",
        default=3,
        check=lambda x: x >= 1
    )
    maxKsumSigma = pexConfig.Field(
        dtype=float,
        doc="""Maximum allowed sigma for outliers from kernel sum distribution.
                 Used to reject variable objects from the kernel model""",
        default=3.0,
        check=lambda x: x >= 0.0
    )
    maxConditionNumber = pexConfig.Field(
        dtype=float,
        doc="Maximum condition number for a well conditioned matrix",
        default=5.0e7,
        check=lambda x: x >= 0.0
    )
    conditionNumberType = pexConfig.ChoiceField(
        dtype=str,
        doc="Use singular values (SVD) or eigen values (EIGENVALUE) to determine condition number",
        default="EIGENVALUE",
        allowed={
            "SVD": "Use singular values",
            "EIGENVALUE": "Use eigen values (faster)",
        }
    )
    maxSpatialConditionNumber = pexConfig.Field(
        dtype=float,
        doc="Maximum condition number for a well conditioned spatial matrix",
        default=1.0e10,
        check=lambda x: x >= 0.0
    )
    iterateSingleKernel = pexConfig.Field(
        dtype=bool,
        doc="""Remake KernelCandidate using better variance estimate after first pass?
                 Primarily useful when convolving a single-depth image, otherwise not necessary.""",
        default=False,
    )
    constantVarianceWeighting = pexConfig.Field(
        dtype=bool,
        doc="""Use constant variance weighting in single kernel fitting?
                 In some cases this is better for bright star residuals.""",
        default=True,
    )
    calculateKernelUncertainty = pexConfig.Field(
        dtype=bool,
        doc="""Calculate kernel and background uncertainties for each kernel candidate?
                 This comes from the inverse of the covariance matrix.
                 Warning: regularization can cause problems for this step.""",
        default=False,
    )
    useBicForKernelBasis = pexConfig.Field(
        dtype=bool,
        doc="""Use Bayesian Information Criterion to select the number of bases going into the kernel""",
        default=False,
    )


class PsfMatchConfigAL(PsfMatchConfig):
    """The parameters specific to the "Alard-Lupton" (sum-of-Gaussian) Psf-matching basis"""

    def setDefaults(self):
        PsfMatchConfig.setDefaults(self)
        self.kernelBasisSet = "alard-lupton"
        self.maxConditionNumber = 5.0e7

    alardNGauss = pexConfig.Field(
        dtype=int,
        doc="Number of Gaussians in alard-lupton basis",
        default=3,
        check=lambda x: x >= 1
    )
    alardDegGauss = pexConfig.ListField(
        dtype=int,
        doc="Polynomial order of spatial modification of Gaussian basis functions. "
            "List length must be `alardNGauss`.",
        default=(4, 2, 2),
    )
    alardSigGauss = pexConfig.ListField(
        dtype=float,
        doc="Default sigma values in pixels of Gaussian kernel base functions. "
            "List length must be `alardNGauss`.",
        default=(0.7, 1.5, 3.0),
    )
    alardGaussBeta = pexConfig.Field(
        dtype=float,
        doc="""Multiplier of Gaussian sigmas between adjacent elements of the kernel basis list.""",
        default=2.0,
        check=lambda x: x >= 0.0,
    )
    alardMinSig = pexConfig.Field(
        dtype=float,
        doc="""Minimum Sigma (pixels) for Gaussians""",
        default=0.7,
        check=lambda x: x >= 0.25
    )
    alardDegGaussDeconv = pexConfig.Field(
        dtype=int,
        doc="""Degree of spatial modification of ALL gaussians in AL basis during deconvolution""",
        default=3,
        check=lambda x: x >= 1
    )
    alardMinSigDeconv = pexConfig.Field(
        dtype=float,
        doc="Minimum Sigma (pixels) for Gaussians during deconvolution; "
            "make smaller than alardMinSig as this is only indirectly used",
        default=0.4,
        check=lambda x: x >= 0.25
    )
    alardNGaussDeconv = pexConfig.Field(
        dtype=int,
        doc="Number of Gaussians in AL basis during deconvolution",
        default=3,
        check=lambda x: x >= 1
    )


class PsfMatchConfigDF(PsfMatchConfig):
    """The parameters specific to the delta-function (one basis per-pixel) Psf-matching basis"""

    def setDefaults(self):
        PsfMatchConfig.setDefaults(self)
        self.kernelBasisSet = "delta-function"
        self.maxConditionNumber = 5.0e6
        self.usePcaForSpatialKernel = True
        self.subtractMeanForPca = True
        self.useBicForKernelBasis = False

    useRegularization = pexConfig.Field(
        dtype=bool,
        doc="Use regularization to smooth the delta function kernels",
        default=True,
    )
    regularizationType = pexConfig.ChoiceField(
        dtype=str,
        doc="Type of regularization.",
        default="centralDifference",
        allowed={
            "centralDifference": "Penalize second derivative using 2-D stencil of central finite difference",
            "forwardDifference": "Penalize first, second, third derivatives using forward finite differeces"
        }
    )
    centralRegularizationStencil = pexConfig.ChoiceField(
        dtype=int,
        doc="Type of stencil to approximate central derivative (for centralDifference only)",
        default=9,
        allowed={
            5: "5-point stencil including only adjacent-in-x,y elements",
            9: "9-point stencil including diagonal elements"
        }
    )
    forwardRegularizationOrders = pexConfig.ListField(
        dtype=int,
        doc="Array showing which order derivatives to penalize (for forwardDifference only)",
        default=(1, 2),
        itemCheck=lambda x: (x > 0) and (x < 4)
    )
    regularizationBorderPenalty = pexConfig.Field(
        dtype=float,
        doc="Value of the penalty for kernel border pixels",
        default=3.0,
        check=lambda x: x >= 0.0
    )
    lambdaType = pexConfig.ChoiceField(
        dtype=str,
        doc="How to choose the value of the regularization strength",
        default="absolute",
        allowed={
            "absolute": "Use lambdaValue as the value of regularization strength",
            "relative": "Use lambdaValue as fraction of the default regularization strength (N.R. 18.5.8)",
            "minimizeBiasedRisk": "Minimize biased risk estimate",
            "minimizeUnbiasedRisk": "Minimize unbiased risk estimate",
        }
    )
    lambdaValue = pexConfig.Field(
        dtype=float,
        doc="Value used for absolute determinations of regularization strength",
        default=0.2,
    )
    lambdaScaling = pexConfig.Field(
        dtype=float,
        doc="Fraction of the default lambda strength (N.R. 18.5.8) to use.  1e-4 or 1e-5",
        default=1e-4,
    )
    lambdaStepType = pexConfig.ChoiceField(
        dtype=str,
        doc="""If a scan through lambda is needed (minimizeBiasedRisk, minimizeUnbiasedRisk),
                 use log or linear steps""",
        default="log",
        allowed={
            "log": "Step in log intervals; e.g. lambdaMin, lambdaMax, lambdaStep = -1.0, 2.0, 0.1",
            "linear": "Step in linear intervals; e.g. lambdaMin, lambdaMax, lambdaStep = 0.1, 100, 0.1",
        }
    )
    lambdaMin = pexConfig.Field(
        dtype=float,
        doc="""If scan through lambda needed (minimizeBiasedRisk, minimizeUnbiasedRisk),
                 start at this value.  If lambdaStepType = log:linear, suggest -1:0.1""",
        default=-1.0,
    )
    lambdaMax = pexConfig.Field(
        dtype=float,
        doc="""If scan through lambda needed (minimizeBiasedRisk, minimizeUnbiasedRisk),
                 stop at this value.  If lambdaStepType = log:linear, suggest 2:100""",
        default=2.0,
    )
    lambdaStep = pexConfig.Field(
        dtype=float,
        doc="""If scan through lambda needed (minimizeBiasedRisk, minimizeUnbiasedRisk),
                 step in these increments.  If lambdaStepType = log:linear, suggest 0.1:0.1""",
        default=0.1,
    )


class PsfMatchTask(pipeBase.Task):
    """Base class for Psf Matching; should not be called directly

    Notes
    -----
    PsfMatchTask is a base class that implements the core functionality for matching the
    Psfs of two images using a spatially varying Psf-matching lsst.afw.math.LinearCombinationKernel.
    The Task requires the user to provide an instance of an lsst.afw.math.SpatialCellSet,
    filled with lsst.ip.diffim.KernelCandidate instances, and a list of lsst.afw.math.Kernels
    of basis shapes that will be used for the decomposition.  If requested, the Task
    also performs background matching and returns the differential background model as an
    lsst.afw.math.Kernel.SpatialFunction.

    Invoking the Task

    As a base class, this Task is not directly invoked.  However, run() methods that are
    implemented on derived classes will make use of the core _solve() functionality,
    which defines a sequence of lsst.afw.math.CandidateVisitor classes that iterate
    through the KernelCandidates, first building up a per-candidate solution and then
    building up a spatial model from the ensemble of candidates.  Sigma clipping is
    performed using the mean and standard deviation of all kernel sums (to reject
    variable objects), on the per-candidate substamp diffim residuals
    (to indicate a bad choice of kernel basis shapes for that particular object),
    and on the substamp diffim residuals using the spatial kernel fit (to indicate a bad
    choice of spatial kernel order, or poor constraints on the spatial model).  The
    _diagnostic() method logs information on the quality of the spatial fit, and also
    modifies the Task metadata.

    .. list-table:: Quantities set in Metadata
       :header-rows: 1

       * - Parameter
         - Description
       * - `spatialConditionNum`
         - Condition number of the spatial kernel fit
       * - `spatialKernelSum`
         - Kernel sum (10^{-0.4 * ``Delta``; zeropoint}) of the spatial Psf-matching kernel
       * - `ALBasisNGauss`
         - If using sum-of-Gaussian basis, the number of gaussians used
       * - `ALBasisDegGauss`
         - If using sum-of-Gaussian basis, the deg of spatial variation of the Gaussians
       * - `ALBasisSigGauss`
         - If using sum-of-Gaussian basis, the widths (sigma) of the Gaussians
       * - `ALKernelSize`
         - If using sum-of-Gaussian basis, the kernel size
       * - `NFalsePositivesTotal`
         - Total number of diaSources
       * - `NFalsePositivesRefAssociated`
         - Number of diaSources that associate with the reference catalog
       * - `NFalsePositivesRefAssociated`
         - Number of diaSources that associate with the source catalog
       * - `NFalsePositivesUnassociated`
         - Number of diaSources that are orphans
       * - `metric_MEAN`
         - Mean value of substamp diffim quality metrics across all KernelCandidates,
           for both the per-candidate (LOCAL) and SPATIAL residuals
       * - `metric_MEDIAN`
         - Median value of substamp diffim quality metrics across all KernelCandidates,
           for both the per-candidate (LOCAL) and SPATIAL residuals
       * - `metric_STDEV`
         - Standard deviation of substamp diffim quality metrics across all KernelCandidates,
           for both the per-candidate (LOCAL) and SPATIAL residuals

    Debug variables

    The lsst.pipe.base.cmdLineTask.CmdLineTask command line task interface supports a
    flag -d/--debug to import @b debug.py from your PYTHONPATH.  The relevant contents of debug.py
    for this Task include:

    .. code-block:: py

        import sys
        import lsstDebug
        def DebugInfo(name):
            di = lsstDebug.getInfo(name)
            if name == "lsst.ip.diffim.psfMatch":
                # enable debug output
                di.display = True
                # display mask transparency
                di.maskTransparency = 80
                # show all the candidates and residuals
                di.displayCandidates = True
                # show kernel basis functions
                di.displayKernelBasis = False
                # show kernel realized across the image
                di.displayKernelMosaic = True
                # show coefficients of spatial model
                di.plotKernelSpatialModel = False
                # show the bad candidates (red) along with good (green)
                di.showBadCandidates = True
            return di
        lsstDebug.Info = DebugInfo
        lsstDebug.frame = 1

    Note that if you want addional logging info, you may add to your scripts:

    .. code-block:: py

        import lsst.log.utils as logUtils
        logUtils.traceSetAt("ip.diffim", 4)
    """
    ConfigClass = PsfMatchConfig
    _DefaultName = "psfMatch"

    def __init__(self, *args, **kwargs):
        """Create the psf-matching Task

        Parameters
        ----------
        *args
            Arguments to be passed to ``lsst.pipe.base.task.Task.__init__``
        **kwargs
            Keyword arguments to be passed to ``lsst.pipe.base.task.Task.__init__``

        Notes
        -----
        The initialization sets the Psf-matching kernel configuration using the value of
        self.config.kernel.active.  If the kernel is requested with regularization to moderate
        the bias/variance tradeoff, currently only used when a delta function kernel basis
        is provided, it creates a regularization matrix stored as member variable
        self.hMat.
        """
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.kConfig = self.config.kernel.active

        if 'useRegularization' in self.kConfig:
            self.useRegularization = self.kConfig.useRegularization
        else:
            self.useRegularization = False

        if self.useRegularization:
            self.hMat = diffimLib.makeRegularizationMatrix(pexConfig.makePolicy(self.kConfig))

    def _diagnostic(self, kernelCellSet, spatialSolution, spatialKernel, spatialBg):
        """Provide logging diagnostics on quality of spatial kernel fit

        Parameters
        ----------
        kernelCellSet : TYPE
            Cellset that contains the KernelCandidates used in the fitting
        spatialSolution : TYPE
            KernelSolution of best-fit
        spatialKernel : TYPE
            Best-fit spatial Kernel model
        spatialBg : TYPE
            Best-fit spatial background model
        """
        # What is the final kernel sum
        kImage = afwImage.ImageD(spatialKernel.getDimensions())
        kSum = spatialKernel.computeImage(kImage, False)
        self.log.info("Final spatial kernel sum %.3f" % (kSum))

        # Look at how well conditioned the matrix is
        conditionNum = spatialSolution.getConditionNumber(
            getattr(diffimLib.KernelSolution, self.kConfig.conditionNumberType))
        self.log.info("Spatial model condition number %.3e" % (conditionNum))

        if conditionNum < 0.0:
            self.log.warn("Condition number is negative (%.3e)" % (conditionNum))
        if conditionNum > self.kConfig.maxSpatialConditionNumber:
            self.log.warn("Spatial solution exceeds max condition number (%.3e > %.3e)" % (
                conditionNum, self.kConfig.maxSpatialConditionNumber))

        self.metadata.set("spatialConditionNum", conditionNum)
        self.metadata.set("spatialKernelSum", kSum)

        # Look at how well the solution is constrained
        nBasisKernels = spatialKernel.getNBasisKernels()
        nKernelTerms = spatialKernel.getNSpatialParameters()
        if nKernelTerms == 0:  # order 0
            nKernelTerms = 1

        # Not fit for
        nBgTerms = spatialBg.getNParameters()
        if nBgTerms == 1:
            if spatialBg.getParameters()[0] == 0.0:
                nBgTerms = 0

        nGood = 0
        nBad = 0
        nTot = 0
        for cell in kernelCellSet.getCellList():
            for cand in cell.begin(False):  # False = include bad candidates
                nTot += 1
                if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                    nGood += 1
                if cand.getStatus() == afwMath.SpatialCellCandidate.BAD:
                    nBad += 1

        self.log.info("Doing stats of kernel candidates used in the spatial fit.")

        # Counting statistics
        if nBad > 2*nGood:
            self.log.warn("Many more candidates rejected than accepted; %d total, %d rejected, %d used" % (
                          nTot, nBad, nGood))
        else:
            self.log.info("%d candidates total, %d rejected, %d used" % (nTot, nBad, nGood))

        # Some judgements on the quality of the spatial models
        if nGood < nKernelTerms:
            self.log.warn("Spatial kernel model underconstrained; %d candidates, %d terms, %d bases" % (
                nGood, nKernelTerms, nBasisKernels))
            self.log.warn("Consider lowering the spatial order")
        elif nGood <= 2*nKernelTerms:
            self.log.warn("Spatial kernel model poorly constrained; %d candidates, %d terms, %d bases" % (
                nGood, nKernelTerms, nBasisKernels))
            self.log.warn("Consider lowering the spatial order")
        else:
            self.log.info("Spatial kernel model well constrained; %d candidates, %d terms, %d bases" % (
                nGood, nKernelTerms, nBasisKernels))

        if nGood < nBgTerms:
            self.log.warn("Spatial background model underconstrained; %d candidates, %d terms" % (
                nGood, nBgTerms))
            self.log.warn("Consider lowering the spatial order")
        elif nGood <= 2*nBgTerms:
            self.log.warn("Spatial background model poorly constrained; %d candidates, %d terms" % (
                nGood, nBgTerms))
            self.log.warn("Consider lowering the spatial order")
        else:
            self.log.info("Spatial background model appears well constrained; %d candidates, %d terms" % (
                nGood, nBgTerms))

    def _displayDebug(self, kernelCellSet, spatialKernel, spatialBackground):
        """Provide visualization of the inputs and ouputs to the Psf-matching code

        Parameters
        ----------
        kernelCellSet : TYPE
            The SpatialCellSet used in determining the matching kernel and background
        spatialKernel : TYPE
            Spatially varying Psf-matching kernel
        spatialBackground : TYPE
            Spatially varying background-matching function
        """
        import lsstDebug
        displayCandidates = lsstDebug.Info(__name__).displayCandidates
        displayKernelBasis = lsstDebug.Info(__name__).displayKernelBasis
        displayKernelMosaic = lsstDebug.Info(__name__).displayKernelMosaic
        plotKernelSpatialModel = lsstDebug.Info(__name__).plotKernelSpatialModel
        showBadCandidates = lsstDebug.Info(__name__).showBadCandidates
        maskTransparency = lsstDebug.Info(__name__).maskTransparency
        if not maskTransparency:
            maskTransparency = 0
        afwDisplay.setDefaultMaskTransparency(maskTransparency)

        if displayCandidates:
            diutils.showKernelCandidates(kernelCellSet, kernel=spatialKernel, background=spatialBackground,
                                         frame=lsstDebug.frame,
                                         showBadCandidates=showBadCandidates)
            lsstDebug.frame += 1
            diutils.showKernelCandidates(kernelCellSet, kernel=spatialKernel, background=spatialBackground,
                                         frame=lsstDebug.frame,
                                         showBadCandidates=showBadCandidates,
                                         kernels=True)
            lsstDebug.frame += 1
            diutils.showKernelCandidates(kernelCellSet, kernel=spatialKernel, background=spatialBackground,
                                         frame=lsstDebug.frame,
                                         showBadCandidates=showBadCandidates,
                                         resids=True)
            lsstDebug.frame += 1

        if displayKernelBasis:
            diutils.showKernelBasis(spatialKernel, frame=lsstDebug.frame)
            lsstDebug.frame += 1

        if displayKernelMosaic:
            diutils.showKernelMosaic(kernelCellSet.getBBox(), spatialKernel, frame=lsstDebug.frame)
            lsstDebug.frame += 1

        if plotKernelSpatialModel:
            diutils.plotKernelSpatialModel(spatialKernel, kernelCellSet, showBadCandidates=showBadCandidates)

    def _createPcaBasis(self, kernelCellSet, nStarPerCell, policy):
        """Create Principal Component basis

        If a principal component analysis is requested, typically when using a delta function basis,
        perform the PCA here and return a new basis list containing the new principal components.

        Parameters
        ----------
        kernelCellSet : TYPE
            a SpatialCellSet containing KernelCandidates, from which components are derived
        nStarPerCell : TYPE
            the number of stars per cell to visit when doing the PCA
        policy : TYPE
            input policy controlling the single kernel visitor

        Returns
        -------
        nRejectedPca : TYPE
            number of KernelCandidates rejected during PCA loop
        spatialBasisList : TYPE
            basis list containing the principal shapes as Kernels

        Raises
        ------
        RuntimeError
            If the Eigenvalues sum to zero.
        """
        nComponents = self.kConfig.numPrincipalComponents
        imagePca = diffimLib.KernelPcaD()
        importStarVisitor = diffimLib.KernelPcaVisitorF(imagePca)
        kernelCellSet.visitCandidates(importStarVisitor, nStarPerCell)
        if self.kConfig.subtractMeanForPca:
            importStarVisitor.subtractMean()
        imagePca.analyze()

        eigenValues = imagePca.getEigenValues()
        pcaBasisList = importStarVisitor.getEigenKernels()

        eSum = np.sum(eigenValues)
        if eSum == 0.0:
            raise RuntimeError("Eigenvalues sum to zero")
        for j in range(len(eigenValues)):
            log.log("TRACE5." + self.log.getName() + "._solve", log.DEBUG,
                    "Eigenvalue %d : %f (%f)", j, eigenValues[j], eigenValues[j]/eSum)

        nToUse = min(nComponents, len(eigenValues))
        trimBasisList = []
        for j in range(nToUse):
            # Check for NaNs?
            kimage = afwImage.ImageD(pcaBasisList[j].getDimensions())
            pcaBasisList[j].computeImage(kimage, False)
            if not (True in np.isnan(kimage.getArray())):
                trimBasisList.append(pcaBasisList[j])

        # Put all the power in the first kernel, which will not vary spatially
        spatialBasisList = diffimLib.renormalizeKernelList(trimBasisList)

        # New Kernel visitor for this new basis list (no regularization explicitly)
        singlekvPca = diffimLib.BuildSingleKernelVisitorF(spatialBasisList, policy)
        singlekvPca.setSkipBuilt(False)
        kernelCellSet.visitCandidates(singlekvPca, nStarPerCell)
        singlekvPca.setSkipBuilt(True)
        nRejectedPca = singlekvPca.getNRejected()

        return nRejectedPca, spatialBasisList

    def _buildCellSet(self, *args):
        """Fill a SpatialCellSet with KernelCandidates for the Psf-matching process;
        override in derived classes"""
        return

    @pipeBase.timeMethod
    def _solve(self, kernelCellSet, basisList, returnOnExcept=False):
        """Solve for the PSF matching kernel

        Parameters
        ----------
        kernelCellSet : `lsst.afw.math.SpatialCellSet`
            a SpatialCellSet to use in determining the matching kernel
             (typically as provided by _buildCellSet)
        basisList : `list` of `lsst.afw.math.kernel.FixedKernel`
            list of Kernels to be used in the decomposition of the spatially varying kernel
            (typically as provided by makeKernelBasisList)
        returnOnExcept : `bool`, optional
            if True then return (None, None) if an error occurs, else raise the exception

        Returns
        -------
        psfMatchingKernel : TYPE
            PSF matching kernel
        backgroundModel : TYPE
            differential background model

        Raises
        ------
        RuntimeError :
            If unable to determine PSF matching kernel and `returnOnExcept==False`
        """

        import lsstDebug
        display = lsstDebug.Info(__name__).display

        maxSpatialIterations = self.kConfig.maxSpatialIterations
        nStarPerCell = self.kConfig.nStarPerCell
        usePcaForSpatialKernel = self.kConfig.usePcaForSpatialKernel

        # Visitor for the single kernel fit
        policy = pexConfig.makePolicy(self.kConfig)
        if self.useRegularization:
            singlekv = diffimLib.BuildSingleKernelVisitorF(basisList, policy, self.hMat)
        else:
            singlekv = diffimLib.BuildSingleKernelVisitorF(basisList, policy)

        # Visitor for the kernel sum rejection
        ksv = diffimLib.KernelSumVisitorF(policy)

        # Main loop
        t0 = time.time()
        try:
            totalIterations = 0
            thisIteration = 0
            while (thisIteration < maxSpatialIterations):

                # Make sure there are no uninitialized candidates as active occupants of Cell
                nRejectedSkf = -1
                while (nRejectedSkf != 0):
                    log.log("TRACE1." + self.log.getName() + "._solve", log.DEBUG,
                            "Building single kernels...")
                    kernelCellSet.visitCandidates(singlekv, nStarPerCell)
                    nRejectedSkf = singlekv.getNRejected()
                    log.log("TRACE1." + self.log.getName() + "._solve", log.DEBUG,
                            "Iteration %d, rejected %d candidates due to initial kernel fit",
                            thisIteration, nRejectedSkf)

                # Reject outliers in kernel sum
                ksv.resetKernelSum()
                ksv.setMode(diffimLib.KernelSumVisitorF.AGGREGATE)
                kernelCellSet.visitCandidates(ksv, nStarPerCell)
                ksv.processKsumDistribution()
                ksv.setMode(diffimLib.KernelSumVisitorF.REJECT)
                kernelCellSet.visitCandidates(ksv, nStarPerCell)

                nRejectedKsum = ksv.getNRejected()
                log.log("TRACE1." + self.log.getName() + "._solve", log.DEBUG,
                        "Iteration %d, rejected %d candidates due to kernel sum",
                        thisIteration, nRejectedKsum)

                # Do we jump back to the top without incrementing thisIteration?
                if nRejectedKsum > 0:
                    totalIterations += 1
                    continue

                # At this stage we can either apply the spatial fit to
                # the kernels, or we run a PCA, use these as a *new*
                # basis set with lower dimensionality, and then apply
                # the spatial fit to these kernels

                if (usePcaForSpatialKernel):
                    log.log("TRACE0." + self.log.getName() + "._solve", log.DEBUG,
                            "Building Pca basis")

                    nRejectedPca, spatialBasisList = self._createPcaBasis(kernelCellSet, nStarPerCell, policy)
                    log.log("TRACE1." + self.log.getName() + "._solve", log.DEBUG,
                            "Iteration %d, rejected %d candidates due to Pca kernel fit",
                            thisIteration, nRejectedPca)

                    # We don't want to continue on (yet) with the
                    # spatial modeling, because we have bad objects
                    # contributing to the Pca basis.  We basically
                    # need to restart from the beginning of this loop,
                    # since the cell-mates of those objects that were
                    # rejected need their original Kernels built by
                    # singleKernelFitter.

                    # Don't count against thisIteration
                    if (nRejectedPca > 0):
                        totalIterations += 1
                        continue
                else:
                    spatialBasisList = basisList

                # We have gotten on to the spatial modeling part
                regionBBox = kernelCellSet.getBBox()
                spatialkv = diffimLib.BuildSpatialKernelVisitorF(spatialBasisList, regionBBox, policy)
                kernelCellSet.visitCandidates(spatialkv, nStarPerCell)
                spatialkv.solveLinearEquation()
                log.log("TRACE2." + self.log.getName() + "._solve", log.DEBUG,
                        "Spatial kernel built with %d candidates", spatialkv.getNCandidates())
                spatialKernel, spatialBackground = spatialkv.getSolutionPair()

                # Check the quality of the spatial fit (look at residuals)
                assesskv = diffimLib.AssessSpatialKernelVisitorF(spatialKernel, spatialBackground, policy)
                kernelCellSet.visitCandidates(assesskv, nStarPerCell)
                nRejectedSpatial = assesskv.getNRejected()
                nGoodSpatial = assesskv.getNGood()
                log.log("TRACE1." + self.log.getName() + "._solve", log.DEBUG,
                        "Iteration %d, rejected %d candidates due to spatial kernel fit",
                        thisIteration, nRejectedSpatial)
                log.log("TRACE1." + self.log.getName() + "._solve", log.DEBUG,
                        "%d candidates used in fit", nGoodSpatial)

                # If only nGoodSpatial == 0, might be other candidates in the cells
                if nGoodSpatial == 0 and nRejectedSpatial == 0:
                    raise RuntimeError("No kernel candidates for spatial fit")

                if nRejectedSpatial == 0:
                    # Nothing rejected, finished with spatial fit
                    break

                # Otherwise, iterate on...
                thisIteration += 1

            # Final fit if above did not converge
            if (nRejectedSpatial > 0) and (thisIteration == maxSpatialIterations):
                log.log("TRACE1." + self.log.getName() + "._solve", log.DEBUG, "Final spatial fit")
                if (usePcaForSpatialKernel):
                    nRejectedPca, spatialBasisList = self._createPcaBasis(kernelCellSet, nStarPerCell, policy)
                regionBBox = kernelCellSet.getBBox()
                spatialkv = diffimLib.BuildSpatialKernelVisitorF(spatialBasisList, regionBBox, policy)
                kernelCellSet.visitCandidates(spatialkv, nStarPerCell)
                spatialkv.solveLinearEquation()
                log.log("TRACE2." + self.log.getName() + "._solve", log.DEBUG,
                        "Spatial kernel built with %d candidates", spatialkv.getNCandidates())
                spatialKernel, spatialBackground = spatialkv.getSolutionPair()

            spatialSolution = spatialkv.getKernelSolution()

        except Exception as e:
            self.log.error("ERROR: Unable to calculate psf matching kernel")

            log.log("TRACE1." + self.log.getName() + "._solve", log.DEBUG, str(e))
            raise e

        t1 = time.time()
        log.log("TRACE0." + self.log.getName() + "._solve", log.DEBUG,
                "Total time to compute the spatial kernel : %.2f s", (t1 - t0))

        if display:
            self._displayDebug(kernelCellSet, spatialKernel, spatialBackground)

        self._diagnostic(kernelCellSet, spatialSolution, spatialKernel, spatialBackground)

        return spatialSolution, spatialKernel, spatialBackground


PsfMatch = PsfMatchTask
