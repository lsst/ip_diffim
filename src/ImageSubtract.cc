// -*- lsst-c++ -*-
/**
 * @file ImageSubtract.cc
 *
 * @brief Implementation of image subtraction functions declared in ImageSubtract.h
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */
#include <iostream>
#include <numeric>
#include <limits>

#include "boost/timer.hpp" 

#include "Eigen/Core"

#include "lsst/afw/image.h"
#include "lsst/afw/math.h"

#include "lsst/pex/logging/Trace.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/exceptions/Runtime.h"

#include "lsst/ip/diffim.h"

namespace afwImage   = lsst::afw::image;
namespace afwMath    = lsst::afw::math;
namespace pexExcept  = lsst::pex::exceptions; 
namespace pexLog     = lsst::pex::logging; 
namespace pexPolicy  = lsst::pex::policy; 

namespace lsst { 
namespace ip { 
namespace diffim {

template<typename PixelT>
std::pair<lsst::afw::math::LinearCombinationKernel::Ptr, lsst::afw::math::Kernel::SpatialFunctionPtr>
fitSpatialKernelFromCandidates(
    lsst::afw::math::SpatialCellSet &kernelCells,       ///< A SpatialCellSet containing KernelCandidates
    lsst::pex::policy::Policy const& policy             ///< Policy to control the processing
    ) {
    typedef typename afwImage::Image<afwMath::Kernel::Pixel> ImageT;
    
    /* There are a variety of recipes for creating a spatial kernel which I will
     * outline here :
     *
     * 1a) Using unregularized delta function kernels, run a full spatial model
     where effectively each delta function basis varies spatially
     individually.  While this is the most general solution and may be
     fast due to the specialization of delta-function convolution, it has
     also been shown to lead to noisy kernels.  This is not recommended.
     
     * 1b) Using unregularized delta function kernels, do a PCA of the returned
     Kernels, and use these to create a new basis set.  This requires a
     first call to singleKernelFitter, then an instance of
     KernelPcaVisitor() to do the PCA, creation of a new kFunctor with
     the eigenBases, a new call to singleKernelFitter using these new
     bases then a call to spatialKernelFitter.  It appears that the
     kernels are not self-similar enough to make this a viable solution.
     
     * 2a) Using regularized delta function kernels, run a full spatial model
     where effectively each delta function basis varies spatially
     individually.  This merely requires repeated calls to
     singleKernelFitter and spatialKernelFitter with the supplied
     kFunctor, same as option 1a) and 3).  While this is general and may
     be fast due to the specialized delta-function convolution, we cannot
     enforce that the kernel sum does not vary spatially.
     
     * 2b) Using regularized delta function kernels, do a PCA of the returned
     Kernels, and use these to create a new basis set.  This requires a
     first call to singleKernelFitter, then an instance of
     KernelPcaVisitor() to do the PCA, creation of a new kFunctor with
     the eigenBases, a new call to singleKernelFitter using these new
     bases then a call to spatialKernelFitter.  While this seems somewhat
     circuitous, we should be able to use many fewer basis functions,
     making the final image convolution faster.  We can also enforce that
     the kernel sum does not vary spatially by modifying the eigenBases.
     
     * 3)  Use Alard Lupton basis set.  This merely requires repeated calls to
     singleKernelFitter and spatialKernelFitter with the supplied
     kFunctor.  With these we can enforce that the kernel sum does not
     vary spatially.
     * 
     */
    
    int const maxSpatialIterations    = policy.getInt("maxSpatialIterations");
    int const nStarPerCell            = policy.getInt("nStarPerCell");
    bool const usePcaForSpatialKernel = policy.getBool("usePcaForSpatialKernel");
    bool const subtractMeanForPca     = policy.getBool("subtractMeanForPca");
    
    /* Build basis list here */
    lsst::afw::math::KernelList basisList = makeKernelBasisList(policy);

    /* Build regularization matrix here */
    boost::shared_ptr<Eigen::MatrixXd> hMat;
    bool useRegularization     = policy.getBool("useRegularization");
    std::string kernelBasisSet = policy.getString("kernelBasisSet");

    if (useRegularization) {
        if (kernelBasisSet == "alard-lupton") {
            pexLog::TTrace<1>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                              "Regularization not enabled for Alard-Lupton kernels");
            useRegularization = false;
        }
        else {
            hMat = makeRegularizationMatrix(policy);
        }
    }
        

    boost::timer t;
    t.restart();
    
    afwMath::LinearCombinationKernel::Ptr spatialKernel;
    afwMath::Kernel::SpatialFunctionPtr spatialBackground;
    
    /* Visitor for the single kernel fit */
    detail::BuildSingleKernelVisitor<PixelT> singleKernelFitter(basisList, policy);
    
    /* Visitor for the kernel sum rejection */
    detail::KernelSumVisitor<PixelT> kernelSumVisitor(policy);
    
    /* Main loop; catch any exceptions */
    try {
        int totalIterations = 0;
        for (int i = 0; i < maxSpatialIterations; i++, totalIterations++) {

            /* Make sure there are no uninitialized candidates as current occupant of Cell */
            int nRejectedSkf = -1;
            while (nRejectedSkf != 0) {
                pexLog::TTrace<2>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                                  "Building single kernels...");
                kernelCells.visitCandidates(&singleKernelFitter, nStarPerCell);
                nRejectedSkf = singleKernelFitter.getNRejected();
            }
            
            /* Reject outliers in kernel sum */
            kernelSumVisitor.resetKernelSum();
            kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::AGGREGATE);
            kernelCells.visitCandidates(&kernelSumVisitor, nStarPerCell);
            kernelSumVisitor.processKsumDistribution();
            kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::REJECT);
            kernelCells.visitCandidates(&kernelSumVisitor, nStarPerCell);

            int nRejectedKsum = kernelSumVisitor.getNRejected();
            pexLog::TTrace<2>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                              "Iteration %d, Spatial Iteration %d, rejected %d Kernels", 
                              totalIterations, i, nRejectedKsum);

            if (nRejectedKsum > 0) {
                /* Jump back to the top; don't count against index i */
                i -= 1;
                continue;
            }
                
            /* 
               At this stage we can either apply the spatial fit to the kernels, or
               we run a PCA, use these as a *new* basis set with lower
               dimensionality, and then apply the spatial fit to these kernels 
            */
            afwMath::KernelList spatialBasisList;
            if (usePcaForSpatialKernel) {
                int const nComponents = policy.getInt("numPrincipalComponents");
                
                pexLog::TTrace<5>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                                  "Building Pca Basis");
                afwImage::ImagePca<ImageT> imagePca;
                detail::KernelPcaVisitor<PixelT> importStarVisitor(&imagePca);
                kernelCells.visitCandidates(&importStarVisitor, nStarPerCell);
                if (subtractMeanForPca) {
                    importStarVisitor.subtractMean();
                }
                imagePca.analyze();
                std::vector<typename ImageT::Ptr> eigenImages = imagePca.getEigenImages();
                std::vector<double> eigenValues = imagePca.getEigenValues();
                
                double eSum = std::accumulate(eigenValues.begin(), eigenValues.end(), 0.);
                for(unsigned int j=0; j < eigenValues.size(); j++) {
                    pexLog::TTrace<6>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                                      "Eigenvalue %d : %f (%f \%)", 
                                      j, eigenValues[j], eigenValues[j]/eSum);
                }
                int const nEigen = subtractMeanForPca ? 
                    static_cast<int>(eigenValues.size()) - 1 :
                    static_cast<int>(eigenValues.size());
                int const ncomp  = (nComponents <= 0 || nEigen < nComponents) ? nEigen : nComponents;

                afwMath::KernelList basisListRaw;
                if (subtractMeanForPca) {
                    basisListRaw.push_back(afwMath::Kernel::Ptr(
                                               new afwMath::FixedKernel(
                                                   afwImage::Image<afwMath::Kernel::Pixel>
                                                   (*(importStarVisitor.returnMean()), true))));
                }
                for (int j = 0; j != ncomp; ++j) {
                    /* Any NaN? */
                    afwImage::Image<afwMath::Kernel::Pixel> img = 
                        afwImage::Image<afwMath::Kernel::Pixel>(*eigenImages[j], true);
                    afwMath::Statistics stats = afwMath::makeStatistics(img, afwMath::SUM);
                    
                    if (std::isnan(stats.getValue(afwMath::SUM))) {
                        pexLog::TTrace<2>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                                          "WARNING : NaN in principal component %d; skipping", j);
                    } else {
                        basisListRaw.push_back(afwMath::Kernel::Ptr(
                                                   new afwMath::FixedKernel(img)
                                                   ));
                    }
                }
                /* Put all the power in the first kernel, which will not vary spatially */
                afwMath::KernelList basisListPca = renormalizeKernelList(basisListRaw);
                
                /* New PsfMatchingFunctor and Kernel visitor for this new basis list */
                detail::BuildSingleKernelVisitor<PixelT> singleKernelFitterPca(basisListPca, policy);
                
                pexLog::TTrace<2>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                                  "Rebuilding kernels using Pca basis");
                
                /* Only true for the first visit so we rebuild each good kernel with
                 * its PcaBasis representation */
                singleKernelFitterPca.setSkipBuilt(false);
                kernelCells.visitCandidates(&singleKernelFitterPca, nStarPerCell);
                /* Once they are built we don't have to revisit */
                singleKernelFitterPca.setSkipBuilt(true);

                int nRejectedPca = singleKernelFitterPca.getNRejected();
                pexLog::TTrace<2>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                                  "Iteration %d, Spatial Iteration %d, rejected %d Kernels", 
                                  totalIterations, i, nRejectedPca);
                
                if (nRejectedPca > 0) {
                    /* We don't want to continue on (yet) with the spatial modeling,
                     * because we have bad objects contributing to the Pca basis.
                     * We basically need to restart from the beginning of this loop,
                     * since the cell-mates of those objects that were rejected need
                     * their original Kernels built by singleKernelFitter.  
                     */

                    /* Jump back to the top; don't count against index i */
                    i -= 1;
                    continue;
                }
                spatialBasisList = basisListPca;
            }
            else {
                spatialBasisList = basisList;
            }

                
            /* We have gotten on to the spatial modeling part */
            detail::BuildSpatialKernelVisitor<PixelT> spatialKernelFitter(spatialBasisList, 
                                                                          policy);
            kernelCells.visitCandidates(&spatialKernelFitter, nStarPerCell);
            spatialKernelFitter.solveLinearEquation();
            pexLog::TTrace<3>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                              "Spatial kernel built with %d candidates", 
                              spatialKernelFitter.getNCandidates());

            std::pair<afwMath::LinearCombinationKernel::Ptr, afwMath::Kernel::SpatialFunctionPtr> kb =
                spatialKernelFitter.getSolutionPair();

            spatialKernel     = kb.first;
            spatialBackground = kb.second;
            
            /* Visitor for the spatial kernel result */
            detail::AssessSpatialKernelVisitor<PixelT> spatialKernelAssessor(spatialKernel, 
                                                                             spatialBackground, 
                                                                             policy);
            kernelCells.visitCandidates(&spatialKernelAssessor, nStarPerCell);
            int nRejectedSpatial = spatialKernelAssessor.getNRejected();
            pexLog::TTrace<1>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                              "Spatial Kernel iteration %d, rejected %d Kernels, using %d Kernels", 
                              i+1, nRejectedSpatial, spatialKernelAssessor.getNGood());

            if (nRejectedSpatial == 0) {
                /* Nothing rejected, finished with spatial fit */
                break;
            }
        }
    } catch (pexExcept::Exception &e) {
        LSST_EXCEPT_ADD(e, "Unable to calculate spatial kernel model");
        throw e; 
    }
    double time = t.elapsed();
    pexLog::TTrace<1>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                      "Total time to compute the spatial kernel : %.2f s", time);
    pexLog::TTrace<2>("lsst.ip.diffim.fitSpatialKernelFromCandidates", "");
    
    return std::make_pair(spatialKernel, spatialBackground);
}

/**
 * @brief Turns Image into a 2-D Eigen Matrix
 */
template <typename PixelT>
Eigen::MatrixXd imageToEigenMatrix(
    lsst::afw::image::Image<PixelT> const &img
    ) {
    unsigned int rows = img.getHeight();
    unsigned int cols = img.getWidth();
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(rows, cols);
    for (int y = 0; y != img.getHeight(); ++y) {
        int x = 0;
        for (typename afwImage::Image<PixelT>::x_iterator ptr = img.row_begin(y); 
             ptr != img.row_end(y); ++ptr, ++x) {
            // M is addressed row, col
            M(y,x) = *ptr;
        }
    }
    return M;
}
    

/** 
 * @brief Implement fundamental difference imaging step of convolution and
 * subtraction : D = I - (K*T + bg) where * denotes convolution
 * 
 * @note If you convolve the science image, D = (K*I + bg) - T, set invert=False
 *
 * @note The template is taken to be an MaskedImage; this takes c 1.6 times as long
 * as using an Image.
 *
 * @note Instantiated such that background can be a double or Function2D
 *
 * @return Difference image
 *
 * @ingroup diffim
 */
template <typename PixelT, typename BackgroundT>
afwImage::MaskedImage<PixelT> convolveAndSubtract(
    lsst::afw::image::MaskedImage<PixelT> const &imageToConvolve,    ///< Image T to convolve with Kernel
    lsst::afw::image::MaskedImage<PixelT> const &imageToNotConvolve, ///< Image I to subtract T from
    lsst::afw::math::Kernel const &convolutionKernel,                ///< PSF-matching Kernel used
    BackgroundT background,                                  ///< Differential background 
    bool invert                                              ///< Invert the output difference image
    ) {

    boost::timer t;
    t.restart();

    afwImage::MaskedImage<PixelT> convolvedMaskedImage(imageToConvolve.getDimensions());
    convolvedMaskedImage.setXY0(imageToConvolve.getXY0());
    afwMath::convolve(convolvedMaskedImage, imageToConvolve, convolutionKernel, false);
    
    /* Add in background */
    *(convolvedMaskedImage.getImage()) += background;
    
    /* Do actual subtraction */
    convolvedMaskedImage -= imageToNotConvolve;

    /* Invert */
    if (invert) {
        convolvedMaskedImage *= -1.0;
    }

    double time = t.elapsed();
    pexLog::TTrace<5>("lsst.ip.diffim.convolveAndSubtract", 
                      "Total compute time to convolve and subtract : %.2f s", time);

    return convolvedMaskedImage;
}

/** 
 * @brief Implement fundamental difference imaging step of convolution and
 * subtraction : D = I - (K.x.T + bg)
 *
 * @note The template is taken to be an Image, not a MaskedImage; it therefore
 * has neither variance nor bad pixels
 *
 * @note If you convolve the science image, D = (K*I + bg) - T, set invert=False
 * 
 * @note Instantiated such that background can be a double or Function2D
 *
 * @return Difference image
 *
 * @ingroup diffim
 */
template <typename PixelT, typename BackgroundT>
afwImage::MaskedImage<PixelT> convolveAndSubtract(
    lsst::afw::image::Image<PixelT> const &imageToConvolve,          ///< Image T to convolve with Kernel
    lsst::afw::image::MaskedImage<PixelT> const &imageToNotConvolve, ///< Image I to subtract T from
    lsst::afw::math::Kernel const &convolutionKernel,                ///< PSF-matching Kernel used
    BackgroundT background,                                  ///< Differential background 
    bool invert                                              ///< Invert the output difference image
    ) {
    
    boost::timer t;
    t.restart();

    afwImage::MaskedImage<PixelT> convolvedMaskedImage(imageToConvolve.getDimensions());
    convolvedMaskedImage.setXY0(imageToConvolve.getXY0());
    afwMath::convolve(*convolvedMaskedImage.getImage(), imageToConvolve, convolutionKernel, false);
    
    /* Add in background */
    *(convolvedMaskedImage.getImage()) += background;
    
    /* Do actual subtraction */
    *convolvedMaskedImage.getImage() -= *imageToNotConvolve.getImage();

    /* Invert */
    if (invert) {
        *convolvedMaskedImage.getImage() *= -1.0;
    }
    *convolvedMaskedImage.getMask() <<= *imageToNotConvolve.getMask();
    *convolvedMaskedImage.getVariance() <<= *imageToNotConvolve.getVariance();
    
    double time = t.elapsed();
    pexLog::TTrace<5>("lsst.ip.diffim.convolveAndSubtract", 
                      "Total compute time to convolve and subtract : %.2f s", time);

    return convolvedMaskedImage;
}

/***********************************************************************************************************/
//
// Explicit instantiations
//

template 
Eigen::MatrixXd imageToEigenMatrix(lsst::afw::image::Image<float> const &);

template 
Eigen::MatrixXd imageToEigenMatrix(lsst::afw::image::Image<double> const &);

template class FindSetBits<lsst::afw::image::Mask<> >;
template class ImageStatistics<float>;
template class ImageStatistics<double>;

/* */

#define p_INSTANTIATE_convolveAndSubtract(TEMPLATE_IMAGE_T, TYPE)     \
    template \
    lsst::afw::image::MaskedImage<TYPE> convolveAndSubtract(            \
        lsst::afw::image::TEMPLATE_IMAGE_T<TYPE> const& imageToConvolve, \
        lsst::afw::image::MaskedImage<TYPE> const& imageToNotConvolve,  \
        lsst::afw::math::Kernel const& convolutionKernel,               \
        double background,                                              \
        bool invert);                                                   \
    \
    template \
    afwImage::MaskedImage<TYPE> convolveAndSubtract( \
        lsst::afw::image::TEMPLATE_IMAGE_T<TYPE> const& imageToConvolve, \
        lsst::afw::image::MaskedImage<TYPE> const& imageToNotConvolve, \
        lsst::afw::math::Kernel const& convolutionKernel, \
        lsst::afw::math::Function2<double> const& backgroundFunction, \
        bool invert); \

#define INSTANTIATE_convolveAndSubtract(TYPE) \
p_INSTANTIATE_convolveAndSubtract(Image, TYPE) \
p_INSTANTIATE_convolveAndSubtract(MaskedImage, TYPE)
/*
 * Here are the instantiations.
 *
 * Do we need double diffim code?  It isn't sufficient to remove it here; you'll have to also remove at
 * least SpatialModelKernel<double> and swig instantiations thereof
 */
INSTANTIATE_convolveAndSubtract(float);
INSTANTIATE_convolveAndSubtract(double);

template
std::pair<lsst::afw::math::LinearCombinationKernel::Ptr, lsst::afw::math::Kernel::SpatialFunctionPtr>
fitSpatialKernelFromCandidates<float>(lsst::afw::math::SpatialCellSet &,
                                      lsst::pex::policy::Policy const&);

}}} // end of namespace lsst::ip::diffim
