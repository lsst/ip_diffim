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
import sys
import diffimLib
import lsst.pex.logging as pexLog
import lsst.pex.exceptions as pexExcept
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

class WarpConfig(pexConfig.Config):
    warpingKernelName = pexConfig.ChoiceField(
        dtype = str,
        doc = "Warping kernel",
        default = "lanczos4",
        allowed = {
            "bilinear": "bilinear interpolation",
            "lanczos3": "Lanczos kernel of order 3",
            "lanczos4": "Lanczos kernel of order 4",
            "lanczos5": "Lanczos kernel of order 5",
        }
    )
    interpLength = pexConfig.Field(
        dtype = int,
        doc = "interpLength argument to lsst.afw.math.warpExposure",
        default = 10
    )
    cacheSize = pexConfig.Field(
        dtype = int,
        doc = "cacheSize argument to lsst.afw.math.SeparableKernel.computeCache",
        default = 0
    )

class DetectionConfig(pexConfig.Config):
    detThreshold = pexConfig.Field(
        dtype = float,
        doc = "Value of footprint detection threshold"
        default = 10.0
    )
    detThresholdType = pexConfig.ChoiceField(
        dtype = str,
        doc = "Type of detection threshold"
        default = "stdev"
        allowed = {
           "value"    : "Use counts as the detection threshold type",
           "stdev"    : "Use standard deviation as the detection threshold type"
           "variance" : "Use variance as the detection threshold type"
        }
    )
    detOnTemplate = pexConfig.Field(
        dtype = bool,
        doc = "If true run detection on the template (imageToConvolve); 
               if false run detection on the science image (imageToNotConvolve)"
        default = True
    )
    detBadMaskPlanes = pexConfig.ListField(
        dtype = str,
        doc = "Mask planes that lead to an invalid detection.
               Options: EDGE SAT BAD CR INTRP
               E.g. : EDGE SAT BAD allows CR-masked and interpolated pixels",
        default = ("EDGE", "SAT", "BAD")
    )
    fpNpixMin = pexConfig.Field(
        dtype = int,
        doc = "Minimum number of pixels in an acceptable Footprint",
        default = 5
    )
    fpNpixMax = pexConfig.Field(
        dtype = int
        doc = "Maximum number of pixels in an acceptable Footprint;
               too big and the subsequent convolutions become unwieldy",
        default = 500
    )
    fpGrowFwhmScaling = pexConfig.Field(
        dtype = float,
        doc = "Grow the footprint based on the Psf Fwhm;
               should be larger than kernelRadiusFwhmScaling",
        default = 10.
    }
    fpGrowPix = pexConfig.Field(
        dtype = int,
        doc = "Grow each raw detection footprint by
               this many pixels.  The smaller the faster; however the kernel sum does
               not converge if the stamp is too small; and the kernel is not
               constrained at all if the stamp is the size of the kernel.  Rule of
               thumb is at least 1.5 times the kernel size.  The grown stamp is
               ~2*fpGrowPix pixels larger in each dimension.",
        default = 30,
        check = lambda x: (x >= 20) && (x <= 40)
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

class PsfMatchConfigAL(pexConfig.Config):
    pass
class PsfMatchConfigDF(pexConfig.Config):
    pass

class PsfMatchConfig(pexConfig.Config):
    warpConfig = pexConfig.ConfigField("Controlling the astrometric WCS remapping", WarpConfig)
    detectionConfig = pexConfig.ConfigField("Controlling the detection of sources for kernel building", DetectionConfig)
    afwBackgroundConfig = pexConfig.ConfigField("Controlling the Afw background fitting", AfwBackgroundConfig)

    ###
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
            "alard-lupton" : "Alard-Lupton sum-of-gaussians basis set,
                           * The first term has no spatial variation
                           * The kernel sum is conserved
                           * You may want to turn off 'usePcaForSpatialKernel'",
            "delta-function" : "Delta-function kernel basis set,
                           * You may enable the option useRegularization
                           * You should seriously consider usePcaForSpatialKernel, which will also
                             enable kernel sum conservation for the delta function kernels"
        }
    )













    ######
    #
    # Do we use createDefaultPolicy to modify terms based on FWHM?
    #
    scaleByFwhm = pexConfig.Field(
        dtype = bool,
        doc = "Scale kernelSize, alardSigGauss, and fpGrowPix by input Fwhm",
        default = false,
    )



    ######
    #
    # Kernel size
    #
    kernelSize = pexConfig.Field(
        dtype = int,
        doc = "Number of rows/columns in the convolution kernel; odd-valued.,
                      Modified by kernelSizeFwhmScaling if scaleByFwhm = true"
        default = 19,
    )

    kernelSizeFwhmScaling = pexConfig.Field(
        dtype = double,
        doc = "How much to scale the kernel size based on the Psf Fwhm;,
                      should be smaller than fpGrowFwhmScaling.  Sets kernelSize."
        default = 4.0,
    )

    kernelSizeMin = pexConfig.Field(
        dtype = int,
        doc = "Minimum kernel dimensions",
        default = 7,
    )

    kernelSizeMax = pexConfig.Field(
        dtype = int,
        doc = "Maximum kernel dimensions",
        default = 31,
    )

    ######
    #
    # Alard-Lupton Basis Parameters
    #
    alardNGauss = pexConfig.Field(
        dtype = int,
        doc = "Number of gaussians in alard-lupton basis",
        default = 3,
    )

    alardDegGauss = pexConfig.Field(
        dtype = int,
        doc = "Degree of spatial modification of gaussians in alard-lupton basis",
        default = 4 3 2,
    )

    alardSigGauss = pexConfig.Field(
        dtype = double,
        doc = "Sigma in pixels of gaussians in alard-lupton basis (note: FWHM = 2.35 sigma). ,
                      Scaled by alardSigFwhmScaling if scaleByFwhm = true"
        default = 0.7 1.5 3.0,
    )

    alardSigFwhmScaling = pexConfig.Field(
        dtype = double,
        doc = "Scaling of the alard-lupton gaussian sigmas.  Sets alardSigGauss",
        default = 0.50 1.00 2.00,
    )

    ######
    #
    # Delta Function Basis Parameters
    #
    useRegularization = pexConfig.Field(
        dtype = bool,
        doc = "Use regularization to smooth the delta function kernels",
        default = true,
    )
    
    regularizationType = pexConfig.Field(
        dtype = string,
        doc = "Type of regularization.",
        default = "centralDifference",
        allowed = pexConfig.Field(
            value:        "centralDifference"
            doc =  "Penalize second derivative using 2-D stencil",
        )
        allowed = pexConfig.Field(
            value:        "forwardDifference"
            doc =  "Penalize first, second, third or combination of derivatives",
        )
    )

    centralRegularizationStencil = pexConfig.Field(
        dtype = int,
        doc = "Type of stencil to approximate central derivative (for centralDifference only)",
        default = 9,
        allowed = pexConfig.Field(
            value: 5
            doc = "5-point stencil including only adjacent-in-x,y elements",
        )
        allowed = pexConfig.Field(
            value: 9
            doc = "9-point stencil including diagonal elements",
        )
    )

    forwardRegularizationOrders = pexConfig.Field(
        dtype = int,
        doc = "Array showing which order derivatives to penalize (for forwardDifference only)",
        default = 1 2,
    )

    regularizationBorderPenalty = pexConfig.Field(
        dtype = double,
        doc = "Value of the penalty for kernel border pixels",
        default = 3.0,
    )

    regularizationScaling = pexConfig.Field(
        dtype = double,
        doc = "Fraction of the default lambda strength (N.R. 18.5.8) to use. ,
                      somewhere around 1e-4 to 1e-5 seems to work.
                      some kernels need high freq power"
        default = 1e-4,
    )

    lambdaType = pexConfig.Field(
        dtype = string,
        doc = "How to choose the value of the regularization strength",
        default = "absolute",
        allowed = pexConfig.Field(
            value:        "absolute"
            doc =  "Use lambdaValue as the value of regularization strength",
        )
        allowed = pexConfig.Field(
            value:        "relative"
            doc =  "Use lambdaValue as fraction of the default regularization strength (N.R. 18.5.8)",
        )
        allowed = pexConfig.Field(
            value:        "minimizeBiasedRisk"
            doc =  "Minimize biased risk estimate",
        )
        allowed = pexConfig.Field(
            value:        "minimizeUnbiasedRisk"
            doc =  "Minimize unbiased risk estimate",
        )
    )
    
    lambdaValue = pexConfig.Field(
        dtype = double,
        doc = "Value used for absolute or relative determinations of regularization strength",
        default = 0.2,
    )

    lambdaStepType = pexConfig.Field(
        dtype = string,
        doc = "If scan through lambda needed (minimizeBiasedRisk, minimizeUnbiasedRisk) use log,
                      or linear steps"
        default = "log",
        allowed = pexConfig.Field(
            value: "log"
            doc = "Step in log intervals; e.g. lambdaMin, lambdaMax, lambdaStep = -1.0, 2.0, 0.1",
        )
        allowed = pexConfig.Field(
            value: "linear"
            doc = "Step in linear intervals; e.g. lambdaMin, lambdaMax, lambdaStep = 0.1, 100, 0.1",
        )
    )
    
    lambdaMin = pexConfig.Field(
        dtype = double,
        doc = "If scan through lambda needed (minimizeBiasedRisk, minimizeUnbiasedRisk) ,
                      start at this value.  If lambdaStepType = log:linear, suggest -1:0.1"
        default = -1.0,
    )
    
    lambdaMax = pexConfig.Field(
        dtype = double,
        doc = "If scan through lambda needed (minimizeBiasedRisk, minimizeUnbiasedRisk) ,
                      stop at this value.  If lambdaStepType = log:linear, suggest 2:100"
        default = 2.0,
    )
    
    lambdaStep = pexConfig.Field(
        dtype = double,
        doc = "If scan through lambda needed (minimizeBiasedRisk, minimizeUnbiasedRisk) ,
                      step in these increments.  If lambdaStepType = log:linear, suggest 0.1:0.1"
        default = 0.1,
    )

    ######
    #
    # Spatial modeling
    #
    spatialKernelType = pexConfig.Field(
        dtype = string,
        doc = "Type of spatial function for kernel",
        default = "chebyshev1",
        allowed = pexConfig.Field(
            value:        "chebyshev1"
            doc =  "Chebyshev polynomial of the first kind",
        )
        allowed = pexConfig.Field(
            value:        "polynomial"
            doc =  "Standard x,y polynomial",
        )
    )

    spatialKernelOrder = pexConfig.Field(
        dtype = int,
        doc = "Spatial order of convolution kernel variation",
        default = 1,
    )

    spatialBgType = pexConfig.Field(
        dtype = string,
        doc = "Type of spatial function for kernel",
        default = "chebyshev1",
        allowed = pexConfig.Field(
            value:        "chebyshev1"
            doc =  "Chebyshev polynomial of the first kind",
        )
        allowed = pexConfig.Field(
            value:        "polynomial"
            doc =  "Standard x,y polynomial",
        )
    )

    spatialBgOrder = pexConfig.Field(
        dtype = int,
        doc = "Spatial order of differential background variation",
        default = 1,
    )

    sizeCellX = pexConfig.Field(
        dtype = int,
        doc = "Size (rows) in pixels of each SpatialCell for spatial modeling",
        default = 256,
    )

    sizeCellY = pexConfig.Field(
        dtype = int,
        doc = "Size (columns) in pixels of each SpatialCell for spatial modeling",
        default = 256,
    )

    nStarPerCell = pexConfig.Field(
        dtype = int,
        doc = "Number of KernelCandidates in each SpatialCell to use in the spatial fitting",
        default = 1,
    )

    maxSpatialIterations = pexConfig.Field(
        dtype = int,
        doc = "Maximum number of iterations for rejecting bad KernelCandidates in spatial fitting",
        default = 3,
    )

    ######
    #
    # Spatial modeling; Pca
    #
    usePcaForSpatialKernel = pexConfig.Field(
        dtype = bool,
        doc = "Use Pca to reduce the dimensionality of the kernel basis sets.,
                      This is particularly useful for delta-function kernels.
                      Functionally, after all Cells have their raw kernels determined, we run 
                      a Pca on these Kernels, re-fit the Cells using the eigenKernels and then 
                      fit those for spatial variation using the same technique as for Alard-Lupton kernels.
                      If this option is used, the first term will have no spatial variation and the 
                      kernel sum will be conserved."
        default = false,
    )

    subtractMeanForPca = pexConfig.Field(
        dtype = bool,
        doc = "Subtract off the mean feature before doing the Pca",
        default = true,
    )

    numPrincipalComponents = pexConfig.Field(
        dtype = int,
        doc = "Number of principal components to use for Pca basis, including the ,
                      mean kernel if requested."
        default = 50,
    )

    fracEigenVal = pexConfig.Field(
        dtype = double,
        doc = "At what fraction of the eigenvalues do you cut off the expansion.,
                      Warning: not yet implemented"
        default = 0.99,
    )

    ######
    # 
    # What types of clipping of KernelCandidates to enable
    #
    singleKernelClipping = pexConfig.Field(
        dtype = bool,
        doc = "Do sigma clipping on each raw kernel candidate",
        default = true,
    )

    kernelSumClipping = pexConfig.Field(
        dtype = bool,
        doc = "Do sigma clipping on the ensemble of kernel sums",
        default = true,
    )

    spatialKernelClipping = pexConfig.Field(
        dtype = bool,
        doc = "Do sigma clipping after building the spatial model",
        default = true,
    )

    ######
    # 
    # Clipping of KernelCandidates based on diffim residuals
    #
    candidateResidualMeanMax = pexConfig.Field(
        dtype = double,
        doc = "Rejects KernelCandidates yielding bad difference image quality.,
                      Represents average over pixels of (image/sqrt(variance))."
        default = 0.25,
    )

    candidateResidualStdMax = pexConfig.Field(
        dtype = double,
        doc = "Rejects KernelCandidates yielding bad difference image quality.,
                      Represents stddev over pixels of (image/sqrt(variance))."
        default = 1.50,
    )

    useCoreStats = pexConfig.Field(
        dtype = bool,
        doc = "Use the core of the stamp for the quality statistics, instead of the entire footprint",
        default = true,
    )

    candidateCoreRadius = pexConfig.Field(
        dtype = int,
        doc = "Radius for calculation of stats in 'core' of KernelCandidate diffim.,
                      Total number of pixels used will be (2*radius)**2. 
                      This is used both for 'core' diffim quality as well as ranking of
                      KernelCandidates by their total flux in this core"
        default = 3,
    )

    ######
    # 
    # Clipping of KernelCandidates based on kernel sum distribution
    #
    maxKsumSigma = pexConfig.Field(
        dtype = double,
        doc = "Maximum allowed sigma for outliers from kernel sum distribution.,
                      Used to reject variable objects from the kernel model"
        default = 3.0,
    )

    ######
    # 
    # Clipping of KernelCandidates based on their matrices
    #
    checkConditionNumber = pexConfig.Field(
        dtype = bool,
        doc = "Test for maximum condition number when inverting a kernel matrix?,
                      Anything above the value is not used and the candidate is set as BAD.        
                      Also used to truncate inverse matrix in estimateBiasedRisk.  However,
                      if you are doing any deconvolution you will want to turn this off, or use
                      a large maxConditionNumber"
        default = false,
    )

    maxConditionNumber = pexConfig.Field(
        dtype = double,
        doc = "Maximum condition number for a well conditioned matrix.,
                      Suggested values:
                      * 5.0e6 for 'delta-function' basis
                      * 5.0e7 for 'alard-lupton' basis"
        default = 5.0e7,
    )

    conditionNumberType = pexConfig.Field(
        dtype = string,
        doc = "Use singular values (SVD) or eigen values (EIGENVALUE) to determine condition number",
        allowed = pexConfig.Field(
            value: "SVD"
            doc = "Use singular values",
        )
        allowed = pexConfig.Field(
            value: "EIGENVALUE"
            doc = "Use eigen values (faster)",
        )
        default = "EIGENVALUE",
    )
    
    ######
    # 
    # Fitting of single kernel to object pair
    #
    iterateSingleKernel = pexConfig.Field(
        dtype = bool,
        doc = "Remake single kernel using better variance estimate after first pass?,
                      Primarily useful when convolving a single-depth image, otherwise not necessary."
        default = false,
    )

    constantVarianceWeighting = pexConfig.Field(
        dtype = bool,
        doc = "Use constant variance weighting in single kernel fitting?,
                      In some cases this is better for bright star residuals."
        default = false,
    )

    calculateKernelUncertainty = pexConfig.Field(
        dtype = bool,
        doc = "Calculate kernel and background uncertainties for each kernel candidate?,
                      This comes from the inverse of the covariance matrix.
                      Warning: regularization can cause problems for this step."
        default = false,
    )

    ######
    # 
    # Any modifications by subpolicies
    #
    modifiedForImagePsfMatch = pexConfig.Field(
        dtype = bool,
        doc = "Policy modified for ImagePsfMatch class",
        default = false,
    )

    modifiedForDeconvolution = pexConfig.Field(
        dtype = bool,
        doc = "Policy modified for deconvolution",
        default = false,
    )

    useOuterForDeconv = pexConfig.Field(
        dtype = bool,
        doc = "Use outer fifth gaussian",
        default = false,
    )

    modifiedForModelPsfMatch = pexConfig.Field(
        dtype = bool,
        doc = "Policy modified for ModelPsfMatch class",
        default = false,
    )

    modifiedForSnapSubtraction = pexConfig.Field(
        dtype = bool,
        doc = "Policy modified for subtraction of back-to-back snaps",
        default = false,
    )

    

class PsfMatch(object):
    """Base class for PSF matching
    """
    def __init__(self, policy, logName="lsst.ip.diffim.PsfMatch"):
        """Create a PsfMatchToImage
        
        @param policy: see lsst/ip/diffim/policy/PsfMatchingDictionary.paf
        @param logName: name by which messages are logged
        """
        self._policy = policy
        self._log = pexLog.Log(pexLog.Log.getDefaultLog(), logName)

    def _solve(self, kernelCellSet, returnOnExcept = False):
        """Determine the PSF matching kernel
        
        @param kernelCellSet: a SpatialCellSet to use in determining the PSF matching kernel
            as provided by self._buildCellSet
        @param returnOnExcept: if True then return (None, None) if an error occurs, else raise an exception
        
        @return
        - psfMatchingKernel: PSF matching kernel
        - backgroundModel: differential background model
        
        @raise Exception if unable to determine PSF matching kernel and returnOnExcept False
        """
        try:
            kb = diffimLib.fitSpatialKernelFromCandidates(kernelCellSet, self._policy)
        except pexExcept.LsstCppException, e:
            pexLog.Trace(self._log.getName(), 1, "ERROR: Unable to calculate psf matching kernel")
            pexLog.Trace(self._log.getName(), 2, e.args[0].what())
    
            if returnOnExcept:
                return (None, None)
            else:
                raise

        psfMatchingKernel = kb[0]
        backgroundModel   = kb[1]
    
        # What is the status of the processing?
        nGood = 0
        for cell in kernelCellSet.getCellList():
            for cand in cell.begin(True):
                cand = diffimLib.cast_KernelCandidateF(cand)
                if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                    nGood += 1
        if nGood == 0:
            pexLog.Trace(self._log.getName(), 1, "WARNING")
        pexLog.Trace(self._log.getName(), 1, "Used %d kernels for spatial fit" % (nGood))
    
        return (psfMatchingKernel, backgroundModel)


