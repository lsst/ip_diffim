// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Implementation of SpatialModelKernel class
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include <numeric>

#include "boost/timer.hpp" 

#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/ImagePca.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/afw/detection/Footprint.h"

#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/ip/diffim/BasisSets.h"
#include "lsst/ip/diffim/PsfMatchingFunctor.h"
#include "lsst/ip/diffim/SpatialModelKernel.h"
#include "lsst/ip/diffim/SpatialModelVisitors.h"

namespace afwMath        = lsst::afw::math;
namespace afwImage       = lsst::afw::image;
namespace pexLogging     = lsst::pex::logging; 
namespace pexExcept      = lsst::pex::exceptions; 
namespace pexPolicy      = lsst::pex::policy; 

namespace lsst { 
namespace ip { 
namespace diffim {
            
template <typename PixelT>
KernelCandidate<PixelT>::ImageT::ConstPtr KernelCandidate<PixelT>::getImage() const {
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Kernel to make KernelCandidate Image from");
    }
    return _image;
}
    
template <typename PixelT>
KernelCandidate<PixelT>::ImageT::Ptr KernelCandidate<PixelT>::copyImage() const {
    return typename KernelCandidate<PixelT>::ImageT::Ptr(
        new typename KernelCandidate<PixelT>::ImageT(*getImage(), true)
        );
    /*
      typename KernelCandidate<PixelT>::ImageT::Ptr imcopy(
      new typename KernelCandidate<PixelT>::ImageT(*_image, true)
      );
      return imcopy;
    */
}


template <typename PixelT>
void KernelCandidate<PixelT>::setKernel(lsst::afw::math::Kernel::Ptr kernel) {
    _kernel     = kernel; 
    _haveKernel = true;
    
    setWidth(_kernel->getWidth());
    setHeight(_kernel->getHeight());
    
    typename KernelCandidate<PixelT>::ImageT::Ptr image (
        new typename KernelCandidate<PixelT>::ImageT(_kernel->getDimensions())
        );
    _kSum  = _kernel->computeImage(*image, false);                    
    _image = image;
}

template <typename PixelT>
lsst::afw::math::Kernel::Ptr KernelCandidate<PixelT>::getKernel() const {
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Kernel for KernelCandidate");
    }
    return _kernel;
}

template <typename PixelT>
double KernelCandidate<PixelT>::getBackground() const {
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Background for KernelCandidate");
    }
    return _background;
}

template <typename PixelT>
double KernelCandidate<PixelT>::getKsum() const {
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Ksum for KernelCandidate");
    }
    return _kSum;
}

template <typename PixelT>
lsst::afw::image::MaskedImage<PixelT> KernelCandidate<PixelT>::returnDifferenceImage() {
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Kernel for KernelCandidate");
    }
    return returnDifferenceImage(_kernel, _background);
}

template <typename PixelT>
lsst::afw::image::MaskedImage<PixelT> KernelCandidate<PixelT>::returnDifferenceImage(
    lsst::afw::math::Kernel::Ptr kernel,
    double background
    ) {
    
    /* Make diffim and set chi2 from result */
    afwImage::MaskedImage<PixelT> diffIm = convolveAndSubtract(*_miToConvolvePtr,
                                                               *_miToNotConvolvePtr,
                                                               *kernel,
                                                               background);
    return diffIm;
    
}
    
/***********************************************************************************************************/
    
template<typename PixelT>
std::pair<lsst::afw::math::LinearCombinationKernel::Ptr, lsst::afw::math::Kernel::SpatialFunctionPtr>
fitSpatialKernelFromCandidates(
    PsfMatchingFunctor<PixelT> &kFunctor,               ///< kFunctor used to build the kernels
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
     SetPcaImageVisitor() to do the PCA, creation of a new kFunctor with
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
     SetPcaImageVisitor() to do the PCA, creation of a new kFunctor with
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
    
    
    int const maxKsumIterations       = policy.getInt("maxKsumIterations");
    int const maxSpatialIterations    = policy.getInt("maxSpatialIterations");
    int const nStarPerCell            = policy.getInt("nStarPerCell");
    int const spatialKernelOrder      = policy.getInt("spatialKernelOrder");
    int const spatialBgOrder          = policy.getInt("spatialBgOrder");
    bool const usePcaForSpatialKernel = policy.getBool("usePcaForSpatialKernel");
    
    boost::timer t;
    t.restart();
    
    /* The basis used for the spatial fit may change if we run Pca */
    boost::shared_ptr<afwMath::KernelList> basisListToUse;
    
    afwMath::LinearCombinationKernel::Ptr spatialKernel;
    afwMath::Kernel::SpatialFunctionPtr spatialBackground;
    
    /* Visitor for the single kernel fit */
    detail::BuildSingleKernelVisitor<PixelT> singleKernelFitter(kFunctor, policy);
    
    /* Visitor for the kernel sum rejection */
    detail::KernelSumVisitor<PixelT> kernelSumVisitor(policy);
    
    /* Main loop; catch any exceptions */
    try {
        for (int i=0, nRejected=-1; i < maxSpatialIterations; i++) {
            /* Make sure there are no uninitialized candidates as current occupant of Cell */
            while (nRejected != 0) {
                pexLogging::TTrace<2>("lsst.ip.diffim.BuildSingleKernelVisitor", 
                                      "Building single kernels A...");
                kernelCells.visitCandidates(&singleKernelFitter, nStarPerCell);
                nRejected = singleKernelFitter.getNRejected();
            }
            
            /* Reject outliers in kernel sum */
            for (int j=0; j < maxKsumIterations; j++) {
                kernelSumVisitor.resetDerived();
                kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::AGGREGATE);
                kernelCells.visitCandidates(&kernelSumVisitor, nStarPerCell);
                kernelSumVisitor.processKsumDistribution();
                kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::REJECT);
                kernelCells.visitCandidates(&kernelSumVisitor, nStarPerCell);
                nRejected = kernelSumVisitor.getNRejected();
                pexLogging::TTrace<2>("lsst.ip.diffim.KernelSumVisitor", 
                                      "Ksum Iteration %d, rejected %d Kernels", j+1, nRejected);
                if (nRejected == 0) {
                    break;
                }
                else {
                    /* Make sure there are no uninitialized candidates as current occupant of Cell */
                    while (nRejected != 0) {
                        pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor", 
                                              "Building single kernels B...");
                        kernelCells.visitCandidates(&singleKernelFitter, nStarPerCell);
                        nRejected = singleKernelFitter.getNRejected();
                    }
                }
            }
            
            /* 
               At this stage we can either apply the spatial fit to the kernels, or
               we run a PCA, use these as a *new* basis set with lower
               dimensionality, and then apply the spatial fit to these kernels 
            */
            if (usePcaForSpatialKernel) {
                int const nComponents = policy.getInt("numPrincipalComponents");
                
                pexLogging::TTrace<5>("lsst.ip.diffim.SetPcaImageVisitor", 
                                      "Building Pca Basis");
                afwImage::ImagePca<ImageT> imagePca;
                detail::SetPcaImageVisitor<PixelT> importStarVisitor(&imagePca);
                kernelCells.visitCandidates(&importStarVisitor, nStarPerCell);
                importStarVisitor.subtractMean();
                imagePca.analyze();
                std::vector<typename ImageT::Ptr> eigenImages = imagePca.getEigenImages();
                std::vector<double> eigenValues = imagePca.getEigenValues();
                
                double eSum = std::accumulate(eigenValues.begin(), eigenValues.end(), 0.);
                for(unsigned int j=0; j < eigenValues.size(); j++) {
                    pexLogging::TTrace<6>("lsst.ip.diffim.SetPcaImageVisitor", 
                                          "Eigenvalue %d : %f (%f \%)", 
                                          j, eigenValues[j], eigenValues[j]/eSum);
                }
                
                int const nEigen = static_cast<int>(eigenValues.size()) - 1; /* -1 as we subtracted mean */
                int const ncomp  = (nComponents <= 0 || nEigen < nComponents) ? nEigen : nComponents;
                //
                afwMath::KernelList kernelListRaw;
                kernelListRaw.push_back(afwMath::Kernel::Ptr(
                                            new afwMath::FixedKernel(
                                                afwImage::Image<afwMath::Kernel::Pixel>
                                                (*(importStarVisitor.returnMean()), true))));
                for (int j = 0; j != ncomp; ++j) {
                    /* Any NaN? */
                    afwImage::Image<afwMath::Kernel::Pixel> img = 
                        afwImage::Image<afwMath::Kernel::Pixel>(*eigenImages[j], true);
                    afwMath::Statistics stats = afwMath::makeStatistics(img, afwMath::SUM);
                    
                    if (std::isnan(stats.getValue(afwMath::SUM))) {
                        pexLogging::TTrace<2>("lsst.ip.diffim.SetPcaImageVisitor", 
                                              "WARNING : NaN in principal component %d; skipping", j);
                    } else {
                        kernelListRaw.push_back(afwMath::Kernel::Ptr(
                                                    new afwMath::FixedKernel(img)
                                                    ));
                    }
                }
                /* Put all the power in the first kernel, which will not vary spatially */
                afwMath::KernelList kernelListPca = renormalizeKernelList(kernelListRaw);
                
                /* New PsfMatchingFunctor and Kernel visitor for this new basis list */
                PsfMatchingFunctor<PixelT> kFunctorPca(kernelListPca);
                detail::BuildSingleKernelVisitor<PixelT> singleKernelFitterPca(kFunctorPca, policy);
                
                /* Always true for Pca kernel; leave original kernel alone for future Pca iterations */
                singleKernelFitterPca.setCandidateKernel(false);
                
                pexLogging::TTrace<2>("lsst.ip.diffim.BuildSingleKernelVisitor", 
                                      "Rebuilding kernels using Pca basis");
                
                /* Only true for the first visit so we rebuild each good kernel with
                 * its PcaBasis representation */
                singleKernelFitterPca.setSkipBuilt(false);
                kernelCells.visitCandidates(&singleKernelFitterPca, nStarPerCell);
                
                /* Once they are built we don't have to revisit */
                singleKernelFitterPca.setSkipBuilt(true);
                
                /* Were any rejected?  If so revisit the bad Cells */
                nRejected = singleKernelFitterPca.getNRejected();
                while (nRejected != 0) {
                    pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor", 
                                          "Building single kernels C...");
                    kernelCells.visitCandidates(&singleKernelFitterPca, nStarPerCell);
                    nRejected = singleKernelFitterPca.getNRejected();
                }
                basisListToUse.reset(new afwMath::KernelList(kFunctorPca.getBasisList()));
            }
            else {
                basisListToUse.reset(new afwMath::KernelList(kFunctor.getBasisList()));
            }
            
            detail::BuildSpatialKernelVisitor<PixelT> spatialKernelFitter(*basisListToUse, 
                                                                          spatialKernelOrder, 
                                                                          spatialBgOrder, 
                                                                          policy);
            kernelCells.visitCandidates(&spatialKernelFitter, nStarPerCell);
            spatialKernelFitter.solveLinearEquation();
            pexLogging::TTrace<3>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                                  "Spatial kernel built with %d candidates", 
                                  spatialKernelFitter.getNCandidates());

            std::pair<afwMath::LinearCombinationKernel::Ptr, afwMath::Kernel::SpatialFunctionPtr> KB =
                spatialKernelFitter.getSpatialModel();

            spatialKernel     = KB.first;
            spatialBackground = KB.second;
            
            /* Visitor for the spatial kernel result */
            detail::AssessSpatialKernelVisitor<PixelT> spatialKernelAssessor(spatialKernel, 
                                                                             spatialBackground, 
                                                                             policy);
            kernelCells.visitCandidates(&spatialKernelAssessor, nStarPerCell);
            nRejected = spatialKernelAssessor.getNRejected();
            pexLogging::TTrace<1>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                                  "Spatial Kernel iteration %d, rejected %d Kernels, using %d Kernels", 
                                  i+1, nRejected, spatialKernelAssessor.getNGood());
            if (nRejected == 0) {
                break;
            }
            /* 
             * Don't need to call kernelCells.visitCandidates() here since its done
             * at the begnning of the loop
             */
            
        }
    } catch (pexExcept::Exception &e) {
        LSST_EXCEPT_ADD(e, "Unable to calculate spatial kernel model");
        throw e; 
    }
    double time = t.elapsed();
    pexLogging::TTrace<1>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                          "Total time to compute the spatial kernel : %.2f s", time);
    pexLogging::TTrace<2>("lsst.ip.diffim.fitSpatialKernelFromCandidates", "");
    
    return std::make_pair(spatialKernel, spatialBackground);
}
    
/***********************************************************************************************************/
//
// Explicit instantiations
//

typedef float PixelT;
template class KernelCandidate<PixelT>;

template
std::pair<lsst::afw::math::LinearCombinationKernel::Ptr, lsst::afw::math::Kernel::SpatialFunctionPtr>
fitSpatialKernelFromCandidates<PixelT>(PsfMatchingFunctor<PixelT> &,
                                       lsst::afw::math::SpatialCellSet &,
                                       lsst::pex::policy::Policy const&);

}}} // end of namespace lsst::ip::diffim

