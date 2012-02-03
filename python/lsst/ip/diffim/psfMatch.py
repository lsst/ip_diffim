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
import time
import diffimLib
import lsst.pex.logging as pexLog
import lsst.pex.exceptions as pexExcept
import lsst.pex.config as pexConfig
import lsst.afw.math as afwMath

class DetectionConfig(pexConfig.Config):
    detThreshold = pexConfig.Field(
        dtype = float,
        doc = "Value of footprint detection threshold",
        default = 10.0
    )
    detThresholdType = pexConfig.ChoiceField(
        dtype = str,
        doc = "Type of detection threshold",
        default = "stdev",
        allowed = {
           "value"    : "Use counts as the detection threshold type",
           "stdev"    : "Use standard deviation as the detection threshold type",
           "variance" : "Use variance as the detection threshold type"
        }
    )
    detOnTemplate = pexConfig.Field(
        dtype = bool,
        doc = """If true run detection on the template (imageToConvolve); 
                 if false run detection on the science image (imageToNotConvolve)""",
        default = True
    )
    detBadMaskPlanes = pexConfig.ListField(
        dtype = str,
        doc = """Mask planes that lead to an invalid detection.
                 Options: EDGE SAT BAD CR INTRP
                 E.g. : EDGE SAT BAD allows CR-masked and interpolated pixels""",
        default = ("EDGE", "SAT", "BAD")
    )
    fpNpixMin = pexConfig.Field(
        dtype = int,
        doc = "Minimum number of pixels in an acceptable Footprint",
        default = 5
    )
    fpNpixMax = pexConfig.Field(
        dtype = int,
        doc = """Maximum number of pixels in an acceptable Footprint;
                 too big and the subsequent convolutions become unwieldy""",
        default = 500
    )
    fpGrowFwhmScaling = pexConfig.Field(
        dtype = float,
        doc = """Grow the footprint based on the Psf Fwhm;
                 should be larger than kernelRadiusFwhmScaling""",
        default = 10.
    )
    fpGrowPix = pexConfig.Field(
        dtype = int,
        doc = """Grow each raw detection footprint by
                 this many pixels.  The smaller the faster; however the kernel sum does
                 not converge if the stamp is too small; and the kernel is not
                 constrained at all if the stamp is the size of the kernel.  Rule of
                 thumb is at least 1.5 times the kernel size.  The grown stamp is
                 ~2*fpGrowPix pixels larger in each dimension.""",
        default = 30,
        check = lambda x: (x >= 20) and (x <= 40)
    )
    
class AfwBackgroundConfig(pexConfig.Config):
    algorithmName = pexConfig.ChoiceField(
        dtype = str,
        doc = "How to interpolate the background values",
        default = "NATURAL_SPLINE",
        allowed = {
            "CONSTANT": "CONSTANT",
            "LINEAR": "LINEAR",
            "NATURAL_SPLINE": "NATURAL_SPLINE",
            "CUBIC_SPLINE": "CUBIC_SPLINE",
            "CUBIC_SPLINE_PERIODIC": "CUBIC_SPLINE_PERIODIC",
            "AKIMA_SPLINE": "AKIMA_SPLINE",
            "AKIMA_SPLINE_PERIODIC": "AKIMA_SPLINE_PERIODIC"
        }
    )
    binsize = pexConfig.Field(
        dtype = int,
        doc = "How large of regions should be used for each background point",
        default = 2048
    )
    undersample = pexConfig.ChoiceField(
        dtype = str,
        doc = "What to do if there are not enough regions for the interpolation",
        default = "REDUCE_INTERP_ORDER",
        allowed = {
            "THROW_EXCEPTION" : "THROW_EXCEPTION",
            "REDUCE_INTERP_ORDER" : "REDUCE_INTERP_ORDER",
            "INCREASE_NXNYSAMPLE" : "INCREASE_NXNYSAMPLE"
        }
    )

class PsfMatchConfig(pexConfig.Config):
    warpingConfig = pexConfig.ConfigField("Config for warping exposures to a common alignment", afwMath.warper.WarperConfig)
    detectionConfig = pexConfig.ConfigField("Controlling the detection of sources for kernel building", DetectionConfig)
    afwBackgroundConfig = pexConfig.ConfigField("Controlling the Afw background fitting", AfwBackgroundConfig)

    ####
    # Background fitting
    useAfwBackground = pexConfig.Field(
        dtype = bool,
        doc = "Use afw background subtraction instead of ip_diffim",
        default = False,
    )
    fitForBackground = pexConfig.Field(
        dtype = bool,
        doc = "Include terms (including kernel cross terms) for background in ip_diffim",
        default = False,
    )

    ####
    # Basis set selection
    kernelBasisSet = pexConfig.ChoiceField(
        dtype = str,
        doc = "Type of basis set for PSF matching kernel.",
        default = "alard-lupton",
        allowed = {
            "alard-lupton" : """Alard-Lupton sum-of-gaussians basis set,
                           * The first term has no spatial variation
                           * The kernel sum is conserved
                           * You may want to turn off 'usePcaForSpatialKernel'""",
            "delta-function" : """Delta-function kernel basis set,
                           * You may enable the option useRegularization
                           * You should seriously consider usePcaForSpatialKernel, which will also
                             enable kernel sum conservation for the delta function kernels"""
        }
    )

    ####
    # Kernel size
    kernelSize = pexConfig.Field(
        dtype = int,
        doc = """Number of rows/columns in the convolution kernel; should be odd-valued.
                 Modified by kernelSizeFwhmScaling if scaleByFwhm = true""",
        default = 19,
        check = lambda x: (x >= 7) and (x <= 31)
    )
    scaleByFwhm = pexConfig.Field(
        dtype = bool,
        doc = "Scale kernelSize, alardSigGauss, and fpGrowPix by input Fwhm",
        default = False,
    )
    kernelSizeFwhmScaling = pexConfig.Field(
        dtype = float,
        doc = """How much to scale the kernel size based on the Psf Fwhm;,
                 should be smaller than fpGrowFwhmScaling.  Sets kernelSize.""",
        default = 4.0,
    )

    #####
    # Spatial modeling
    spatialModelType = pexConfig.ChoiceField(
        dtype = str,
        doc = "Type of spatial functions for kernel and background",
        default = "chebyshev1",
        allowed = {
            "chebyshev1" : "Chebyshev polynomial of the first kind",
            "polynomial" : "Standard x,y polynomial",
        }
    )
    spatialKernelOrder = pexConfig.Field(
        dtype = int,
        doc = "Spatial order of convolution kernel variation",
        default = 2,
        check = lambda x : x >= 0
    )
    spatialBgOrder = pexConfig.Field(
        dtype = int,
        doc = "Spatial order of differential background variation",
        default = 1,
        check = lambda x : x >= 0
    )
    sizeCellX = pexConfig.Field(
        dtype = int,
        doc = "Size (rows) in pixels of each SpatialCell for spatial modeling",
        default = 128,
    )
    sizeCellY = pexConfig.Field(
        dtype = int,
        doc = "Size (columns) in pixels of each SpatialCell for spatial modeling",
        default = 128,
    )
    nStarPerCell = pexConfig.Field(
        dtype = int,
        doc = "Number of KernelCandidates in each SpatialCell to use in the spatial fitting",
        default = 3,
    )
    maxSpatialIterations = pexConfig.Field(
        dtype = int,
        doc = "Maximum number of iterations for rejecting bad KernelCandidates in spatial fitting",
        default = 3,
    )

    ####
    # Spatial modeling; Pca
    usePcaForSpatialKernel = pexConfig.Field(
        dtype = bool,
        doc = """Use Pca to reduce the dimensionality of the kernel basis sets.
                 This is particularly useful for delta-function kernels.
                 Functionally, after all Cells have their raw kernels determined, we run 
                 a Pca on these Kernels, re-fit the Cells using the eigenKernels and then 
                 fit those for spatial variation using the same technique as for Alard-Lupton kernels.
                 If this option is used, the first term will have no spatial variation and the 
                 kernel sum will be conserved.""",
        default = False,
    )
    subtractMeanForPca = pexConfig.Field(
        dtype = bool,
        doc = "Subtract off the mean feature before doing the Pca",
        default = True,
    )
    numPrincipalComponents = pexConfig.Field(
        dtype = int,
        doc = """Number of principal components to use for Pca basis, including the
                 mean kernel if requested.""",
        default = 5,
    )
    fracEigenVal = pexConfig.Field(
        dtype = float,
        doc = "At what fraction of the eigenvalues do you cut off the expansion. Warning: not yet implemented",
        default = 0.99,
    )

    ####
    # What types of clipping of KernelCandidates to enable
    singleKernelClipping = pexConfig.Field(
        dtype = bool,
        doc = "Do sigma clipping on each raw kernel candidate",
        default = True,
    )
    kernelSumClipping = pexConfig.Field(
        dtype = bool,
        doc = "Do sigma clipping on the ensemble of kernel sums",
        default = True,
    )
    spatialKernelClipping = pexConfig.Field(
        dtype = bool,
        doc = "Do sigma clipping after building the spatial model",
        default = False,
    )
    checkConditionNumber = pexConfig.Field(
        dtype = bool,
        doc = """Test for maximum condition number when inverting a kernel matrix.
                 Anything above maxConditionNumber is not used and the candidate is set as BAD.        
                 Also used to truncate inverse matrix in estimateBiasedRisk.  However,
                 if you are doing any deconvolution you will want to turn this off, or use
                 a large maxConditionNumber""",
        default = False,
    )

    ####
    # Clipping of KernelCandidates based on diffim residuals; used with singleKernelClipping and spatialKernelClipping
    candidateResidualMeanMax = pexConfig.Field(
        dtype = float,
        doc = """Rejects KernelCandidates yielding bad difference image quality.
                 Used by BuildSingleKernelVisitor, AssessSpatialKernelVisitor.
                 Represents average over pixels of (image/sqrt(variance)).""",
        default = 0.25,
    )
    candidateResidualStdMax = pexConfig.Field(
        dtype = float,
        doc = """Rejects KernelCandidates yielding bad difference image quality.
                 Used by BuildSingleKernelVisitor, AssessSpatialKernelVisitor.
                 Represents stddev over pixels of (image/sqrt(variance)).""",
        default = 1.50,
    )
    useCoreStats = pexConfig.Field(
        dtype = bool,
        doc = """Use the core of the footprint for the quality statistics, instead of the entire footprint.
                 WARNING: if there is deconvolution we probably will need to turn this off""",
        default = True,
    )
    candidateCoreRadius = pexConfig.Field(
        dtype = int,
        doc = """Radius for calculation of stats in 'core' of KernelCandidate diffim.
                 Total number of pixels used will be (2*radius)**2. 
                 This is used both for 'core' diffim quality as well as ranking of
                 KernelCandidates by their total flux in this core""",
        default = 3,
    )

    ####
    # Clipping of KernelCandidates based on kernel sum distribution
    maxKsumSigma = pexConfig.Field(
        dtype = float,
        doc = """Maximum allowed sigma for outliers from kernel sum distribution.
                 Used to reject variable objects from the kernel model""",
        default = 3.0,
    )


    ####
    # Clipping of KernelCandidates based on condition number
    maxConditionNumber = pexConfig.Field(
        dtype = float,
        doc = "Maximum condition number for a well conditioned matrix",
        default = 5.0e7,
    )
    conditionNumberType = pexConfig.ChoiceField(
        dtype = str,
        doc = "Use singular values (SVD) or eigen values (EIGENVALUE) to determine condition number",
        default = "EIGENVALUE",
        allowed = {
            "SVD" : "Use singular values",
            "EIGENVALUE" : "Use eigen values (faster)",
        }
    )
    
    ####
    # Fitting of single kernel to object pair in KernelCandidate
    iterateSingleKernel = pexConfig.Field(
        dtype = bool,
        doc = """Remake KernelCandidate using better variance estimate after first pass?
                 Primarily useful when convolving a single-depth image, otherwise not necessary.""",
        default = False,
    )
    constantVarianceWeighting = pexConfig.Field(
        dtype = bool,
        doc = """Use constant variance weighting in single kernel fitting?
                 In some cases this is better for bright star residuals.""",
        default = False,
    )
    calculateKernelUncertainty = pexConfig.Field(
        dtype = bool,
        doc = """Calculate kernel and background uncertainties for each kernel candidate?
                 This comes from the inverse of the covariance matrix.
                 Warning: regularization can cause problems for this step.""",
        default = False,
    )


class PsfMatchConfigAL(PsfMatchConfig):
    def __init__(self):
        PsfMatchConfig.__init__(self)
        self.kernelBasisSet = "alard-lupton"
        self.maxConditionNumber = 5.0e7

    #####
    # Alard-Lupton Basis Parameters
    alardNGauss = pexConfig.Field(
        dtype = int,
        doc = "Number of gaussians in alard-lupton basis",
        default = 3,
    )
    alardDegGauss = pexConfig.ListField(
        dtype = int,
        doc = "Degree of spatial modification of gaussians in alard-lupton basis",
        default = (4, 3, 2),
    )
    alardSigGauss = pexConfig.ListField(
        dtype = float,
        doc = """Sigma in pixels of gaussians in alard-lupton basis (note: FWHM = 2.35 sigma).""",
        default = (0.7, 1.5, 3.0),
    )

    alardNGaussDeconv = pexConfig.Field(
        dtype = int,
        doc = "Number of gaussians in deconvolving alard-lupton basis",
        default = 4,
    )
    alardDegGaussDeconv = pexConfig.ListField(
        dtype = int,
        doc = """Degree of spatial modification of gaussians in alard-lupton basis during deconvolution""",
        default = (5, 5, 3, 3,),
    )


    
class PsfMatchConfigDF(PsfMatchConfig):
    def __init__(self):
        PsfMatchConfig.__init__(self)
        self.kernelBasisSet = "delta-function"
        self.maxConditionNumber = 5.0e6

    #####
    # Delta Function Basis Parameters
    useRegularization = pexConfig.Field(
        dtype = bool,
        doc = "Use regularization to smooth the delta function kernels",
        default = True,
    )

    #####
    # Regularization shape
    regularizationType = pexConfig.ChoiceField(
        dtype = str,
        doc = "Type of regularization.",
        default = "centralDifference",
        allowed = {
            "centralDifference": "Penalize second derivative using 2-D stencil of central finite difference",
            "forwardDifference": "Penalize first, second, third or combination of derivatives using forward finite differeces"
        }
    )
    centralRegularizationStencil = pexConfig.ChoiceField(
        dtype = int,
        doc = "Type of stencil to approximate central derivative (for centralDifference only)",
        default = 9,
        allowed = { 
            5 : "5-point stencil including only adjacent-in-x,y elements",
            9 : "9-point stencil including diagonal elements"
        }
    )
    forwardRegularizationOrders = pexConfig.ListField(
        dtype = int,
        doc = "Array showing which order derivatives to penalize (for forwardDifference only)",
        default = (1, 2),
        itemCheck = lambda x: (x > 0) and (x < 4)
    )
    regularizationBorderPenalty = pexConfig.Field(
        dtype = float,
        doc = "Value of the penalty for kernel border pixels",
        default = 3.0,
    )

    #####
    # Regularization strength
    lambdaType = pexConfig.ChoiceField(
        dtype = str,
        doc = "How to choose the value of the regularization strength",
        default = "absolute",
        allowed = {
            "absolute" : "Use lambdaValue as the value of regularization strength",
            "relative" : "Use lambdaValue as fraction of the default regularization strength (N.R. 18.5.8)",
            "minimizeBiasedRisk" : "Minimize biased risk estimate",
            "minimizeUnbiasedRisk" : "Minimize unbiased risk estimate",
        }
    )
    lambdaValue = pexConfig.Field(
        dtype = float,
        doc = "Value used for absolute determinations of regularization strength",
        default = 0.2,
    )
    lambdaScaling = pexConfig.Field(
        dtype = float,
        doc = "Fraction of the default lambda strength (N.R. 18.5.8) to use.  1e-4 or 1e-5",
        default = 1e-4,
    )
    lambdaStepType = pexConfig.ChoiceField(
        dtype = str,
        doc = """If a scan through lambda is needed (minimizeBiasedRisk, minimizeUnbiasedRisk),
                 use log or linear steps""",
        default = "log",
        allowed = {
            "log" : "Step in log intervals; e.g. lambdaMin, lambdaMax, lambdaStep = -1.0, 2.0, 0.1",
            "linear" : "Step in linear intervals; e.g. lambdaMin, lambdaMax, lambdaStep = 0.1, 100, 0.1",
        }
    )
    lambdaMin = pexConfig.Field(
        dtype = float,
        doc = """If scan through lambda needed (minimizeBiasedRisk, minimizeUnbiasedRisk),
                 start at this value.  If lambdaStepType = log:linear, suggest -1:0.1""",
        default = -1.0,
    )
    lambdaMax = pexConfig.Field(
        dtype = float,
        doc = """If scan through lambda needed (minimizeBiasedRisk, minimizeUnbiasedRisk),
                 stop at this value.  If lambdaStepType = log:linear, suggest 2:100""",
        default = 2.0,
    )
    lambdaStep = pexConfig.Field(
        dtype = float,
        doc = """If scan through lambda needed (minimizeBiasedRisk, minimizeUnbiasedRisk),
                 step in these increments.  If lambdaStepType = log:linear, suggest 0.1:0.1""",
        default = 0.1,
    )


class PsfMatch(object):
    """Base class for PSF matching
    """
    def __init__(self, config, logName="lsst.ip.diffim.PsfMatch"):
        """Create a PsfMatchToImage
        
        @param config: configuration for ip_diffim
        @param logName: name by which messages are logged
        """
        self._config = config
        self._log = pexLog.Log(pexLog.Log.getDefaultLog(), logName)

    def _diagnostic(self, kernelCellSet, spatialKernel, spatialBg):
        # Look at how well the solution is constrained

        nKernelTerms = spatialKernel.getNSpatialParameters()
        if nKernelTerms == 0: # order 0
            nKernelTerms = 1
        nBgTerms     = len(spatialBg.getParameters())
        
        nGood = 0
        nBad  = 0
        nTot  = 0
        for cell in kernelCellSet.getCellList():
            for cand in cell.begin(False): # False = include bad candidates
                cand = diffimLib.cast_KernelCandidateF(cand)
                nTot += 1
                if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                    nGood += 1
                if cand.getStatus() == afwMath.SpatialCellCandidate.BAD:
                    nBad += 1

        # Counting statistics
        if nBad > 2*nGood:
            pexLog.Trace(self._log.getName()+"._diagnostic", 1,
                         "WARNING: many more candidates rejected than accepted; %d total, %d rejected, %d used" % (
                          nTot, nBad, nGood) )
        else:
            pexLog.Trace(self._log.getName()+"._diagnostic", 1,
                         "NOTE: %d candidates total, %d rejected, %d used" % (nTot, nBad, nGood))
            
        # Some judgements on the quality of the spatial models
        if nGood < nKernelTerms:
            pexLog.Trace(self._log.getName()+"._diagnostic", 1,
                         "WARNING: spatial kernel model underconstrained; %d candidates, %d terms" % (nGood, nKernelTerms))
            pexLog.Trace(self._log.getName()+"._diagnostic", 2,
                         "Consider lowering the spatial order")
        elif nGood <= 2*nKernelTerms:
            pexLog.Trace(self._log.getName()+"._diagnostic", 1,
                         "WARNING: spatial kernel model poorly constrained; %d candidates, %d terms" % (nGood, nKernelTerms))
            pexLog.Trace(self._log.getName()+"._diagnostic", 2,
                         "Consider lowering the spatial order")
        else:
            pexLog.Trace(self._log.getName()+"._diagnostic", 1,
                         "NOTE: spatial kernel model appears well constrained; %d candidates, %d terms" % (nGood, nKernelTerms))

        if nGood < nBgTerms:
            pexLog.Trace(self._log.getName()+"._diagnostic", 1,
                         "WARNING: spatial background model underconstrained; %d candidates, %d terms" % (nGood, nBgTerms))
            pexLog.Trace(self._log.getName()+"._diagnostic", 2,
                         "Consider lowering the spatial order")
        elif nGood <= 2*nBgTerms:
            pexLog.Trace(self._log.getName()+"._diagnostic", 1,
                         "WARNING: spatial background model poorly constrained; %d candidates, %d terms" % (nGood, nBgTerms))
            pexLog.Trace(self._log.getName()+"._diagnostic", 2,
                         "Consider lowering the spatial order")
        else:
            pexLog.Trace(self._log.getName()+"._diagnostic", 1,
                         "NOTE: spatial background model appears well constrained; %d candidates, %d terms" % (nGood, nBgTerms))
        
    
    def _solve(self, kernelCellSet, basisList, hMat = None, returnOnExcept = False):
        """Determine the PSF matching kernel
        
        @param kernelCellSet: a SpatialCellSet to use in determining the PSF matching kernel
            as provided by self._buildCellSet
        @param returnOnExcept: if True then return (None, None) if an error occurs, else raise the exception
        
        @return
        - psfMatchingKernel: PSF matching kernel
        - backgroundModel: differential background model
        
        @raise Exception if unable to determine PSF matching kernel and returnOnExcept False
        """


        maxSpatialIterations   = self._config.maxSpatialIterations
        nStarPerCell           = self._config.nStarPerCell
        usePcaForSpatialKernel = self._config.usePcaForSpatialKernel
        subtractMeanForPca     = self._config.subtractMeanForPca
        useRegularization      = self._config.useRegularization

        # Visitor for the single kernel fit
        if useRegularization:
            singlekv = ipDiffim.BuildSingleKernelVisitorF(basisList, self._config)
        else:
            singlekv = ipDiffim.BuildSingleKernelVisitorF(basisList, self._config, hMat)

        # Visitor for the kernel sum rejection
        ksv = ipDiffim.kernelSumVisitorF(self._config)
        
        # Main loop
        t0 = time.time()
        try:
            totalIterations = 0
            thisIteration   = 0
            while (thisIteration < maxSpatialIterations):

                # Make sure there are no uninitialized candidates as active occupants of Cell
                nRejectedSkf = -1
                while (nRejectedSkf != 0):
                    pexLog.Trace(self._log.getName()+"._solve", 2, 
                                 "Building single kernels...")
                    kernelCells.visitCandidates(singlekv, nStarPerCell)
                    nRejectedSkf = singleKernelFitter.getNRejected()
                    pexLog.Trace(self._log.getName()+"._solve", 2, 
                                 "Iteration %d, rejected %d candidates due to initial kernel fit" % (thisIteration, nRejectedKsum))

                # Reject outliers in kernel sum 
                ksv.resetKernelSum()
                ksv.setMode(ipDiffimDetail.AGGREGATE)
                kernelCells.visitCandidates(ksv, nStarPerCell)
                ksv.processKsumDistribution()
                ksv.setMode(ipDiffimDetail.REJECT)
                kernelCells.visitCandidates(ksv, nStarPerCell)      

                nRejectedKsum = ksv.getNRejected()
                pexLog.Trace(self._log.getName()+"._solve", 2, 
                             "Iteration %d, rejected %d candidates due to kernel sum" % (thisIteration, nRejectedKsum))


                # Do we jump back to the top without incrementing thisIteration?
                if nRejectedKsum > 0:
                    totalIterations += 1
                    continue

                # At this stage we can either apply the spatial fit to
                # the kernels, or we run a PCA, use these as a *new*
                # basis set with lower dimensionality, and then apply
                # the spatial fit to these kernels
                
                if (usePcaForSpatialKernel):
                    pexLog.Trace(self._log.getName()+"._solve", 5, 
                                 "Building Pca basis")

                    nComponents       = self._config.numPrincipalComponents
                    imagePca          = afwImage.ImagePcaF()
                    importStarVisitor = ipDiffim.KernelPcaVisitor(imagePca)
                    kernelCells.visitCandidates(importStarVisitor, nStarPerCell)
                    if self._config.subtractMeanForPca:
                        importStarVisitor.subtractMean()
                    imagePca.analyze()

                    # eigenImages  = imagePca.getEigenImages()
                    eigenValues  = imagePca.getEigenValues()
                    pcaBasisList = importStarVisitor.getEigenKernels()

                    eSum = num.sum(eigenValues)
                    for j in range(len(eSum)):
                        pexLog.Trace(self._log.getName()+"._solve", 6, 
                                     "Eigenvalue %d : %f (%f \%)" % (j, eigenValues[j], eigenValues[j]/eSum))
                                     
                    nToUse = min(nComponents, len(eigenValues))
                    trimBasisList = afwMath.KernelList()
                    for j in range(nToUse):
                        # Check for NaNs!
                        if not (True in num.isnan(pcaBasisList[j])[0]):
                            trimBasisList.push_back(pcaBasisList[j])
                            
                    # Put all the power in the first kernel, which will not vary spatially
                    spatialBasisList = ipDiffim.renormalizeKernelList(trimBasisList)

                    # New Kernel visitor for this new basis list (no regularization explicitly)
                    singlekvPca = ipDiffim.BuildSingleKernelVisitorF(spatialBasisList, self._config)
                    singlekvPca.setSkipBuilt(False)
                    kernelCells.visitCandidates(singlekvPca, nStarPerCell)
                    singlekvPca.setSkipBuilt(True)
                    nRejectedPca = singlekvPca.getNRejected()
                    pexLog.Trace(self._log.getName()+"._solve", 2, 
                                 "Iteration %d, rejected %d candidates due to Pca kernel fit" % (thisIteration, nRejectedPca))
                    
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
                regionBBox = kernelCells.getBBox()
                spatialkv  = ipDiffim.BuildSpatialKernelVisitor(spatialBasisList, regionBBox, self._config)
                kernelCells.visitCandidates(spatialkv, nStarPerCell)
                spatialkv.solveLinearEquation()
                pexLog.Trace(self._log.getName()+"._solve", 3, 
                             "Spatial kernel built with %d candidates" % (spatialkv.getNCandidates()))
                spatialKernel, spatialBackground = spatialkv.getSolutionPair()

                # Check the quality of the spatial fit (look at residuals)
                assesskv   = ipDiffim.AssessSpatialKernelVisitor(spatialKernel, spatialBackground, self._config)
                kernelCells.visitCandidates(assesskv, nStarPerCell)
                nRejectedSpatial = assesskv.getNRejected()
                nGoodSpatial     = assesskv.getNGood()
                pexLog.Trace(self._log.getName()+"._solve", 2, 
                             "Iteration %d, rejected %d candidates due to spatial kernel fit" % (thisIteration, nRejectedPca))
                pexLog.Trace(self._log.getName()+"._solve", 2, 
                             "%d candidates used in fit" % (nGoodSpatial))

                if nGoodSpatial == 0:
                    raise RuntimeError("No kernel candidates for spatial fit")

                if nRejectedSpatial == 0:
                    # Nothing rejected, finished with spatial fit
                    break
                
                # Otherwise, iterate on...
                thisIteration   += 1

        except Exception, e:
            pexLog.Trace(self._log.getName()+"._solve", 1, "ERROR: Unable to calculate psf matching kernel")
            pexLog.Trace(self._log.getName()+"._solve", 2, e.args[0].what())
            raise e
            
      
        t1 = time.time()
        pexLog.Trace(self._log.getName()+"._solve", 1, 
                     "Total time to compute the spatial kernel : %.2f s" % (t1 - t0))
        pexLog.Trace(self._log.getName()+"._solve", 2, "")

        self._diagnostic(spatialKernel, spatialBackground)
        
        return spatialKernel, spatialBackground
