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

#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/Image.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/detection/Footprint.h>
#include <lsst/daf/base.h>

#include <lsst/pex/policy/Policy.h>
#include <lsst/pex/logging/Trace.h>

#include <lsst/ip/diffim/SpatialModelKernel.h>

namespace exceptions = lsst::pex::exceptions; 
namespace policy     = lsst::pex::policy; 
namespace logging    = lsst::pex::logging; 
namespace image      = lsst::afw::image;
namespace math       = lsst::afw::math;
namespace detection  = lsst::afw::detection;

namespace lsst {
namespace ip {
namespace diffim {
    
template <typename ImageT>
SpatialModelKernel<ImageT>::SpatialModelKernel(
    lsst::afw::detection::Footprint::Ptr const &fpPtr,
    MaskedImagePtr const &miToConvolvePtr,
    MaskedImagePtr const &miToNotConvolvePtr,
    boost::shared_ptr<PsfMatchingFunctor<ImageT> > const &kFunctor,
    lsst::pex::policy::Policy const &policy,
    bool build
    ) :
    _fpPtr(fpPtr),
    _miToConvolvePtr(miToConvolvePtr),
    _miToNotConvolvePtr(miToNotConvolvePtr),
    _kFunctor(kFunctor),
    _policy(policy),
    _colc(0.),
    _rowc(0.),
    _kPtr(),
    _kErrPtr(),
    _kSum(0.),
    _bg(0.),
    _bgErr(0.),
    _kStats(),
    _isBuilt(false),
    _isGood(false)
{
    if (build == true) {
        buildModel();
    }
}

template <typename ImageT>
bool SpatialModelKernel<ImageT>::buildModel() {

    if (isBuilt() == true) {
        return false;
    }

    // fill in information on position in the image
    image::BBox fpBBox = _fpPtr->getBBox();

    // NOTE : since we can't remap pixel range to go from -1 to 1 in convolve(),
    // we have to use the actual pixel value here.  Not optimal.
    // setColc(float(fpBBox.getX0() + fpBBox.getX1()) / _miToConvolveParentPtr->getWidth() - 1.0);
    // setRowc(float(fpBBox.getY0() + fpBBox.getY1()) / _miToConvolveParentPtr->getHeight() - 1.0);
    setColc(0.5 * float(fpBBox.getX0() + fpBBox.getX1()));
    setRowc(0.5 * float(fpBBox.getY0() + fpBBox.getY1()));

    logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                       "Footprint = %d,%d -> %d,%d (center %d,%d)",
                       fpBBox.getX0(), fpBBox.getY0(),
                       fpBBox.getX1(), fpBBox.getY1(),
		       int(getColc()), int(getRowc()));

    // Estimate of the variance for first kernel pass
    // True argument is for a deep copy, so -= does not modify the original pixels
    image::MaskedImage<ImageT> varEstimate = image::MaskedImage<ImageT>(*_miToNotConvolvePtr, true);
    varEstimate -= *_miToConvolvePtr;

    try {
        _kFunctor->apply(*_miToConvolvePtr->getImage(), *_miToNotConvolvePtr->getImage(),
                         *varEstimate.getVariance(), _policy);
    } catch (exceptions::Exception& e) {
        setStatus(false);
        logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                           "Exception caught from computePsfMatchingKernelForFootprint"); 
        logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                           e.what());
        return false;
    }
    math::Kernel::PtrT kernelPtr       = _kFunctor->getKernel();
    math::Kernel::PtrT kernelErrorPtr  = _kFunctor->getKernelError();
    double background =      _kFunctor->getBackground();
    double backgroundError = _kFunctor->getBackgroundError();

    // Compute kernel sum
    double kSum = 0.;
    unsigned int kCols = _policy.getInt("kernelCols");
    unsigned int kRows = _policy.getInt("kernelRows");
    image::Image<double> kImage(kCols, kRows);
    kSum = kernelPtr->computeImage(kImage, false);

    // Create difference image and calculate associated statistics
    image::MaskedImage<ImageT> diffIm = convolveAndSubtract(*_miToConvolvePtr->getImage(), *_miToNotConvolvePtr,
                                                            *kernelPtr, background);

    // Find statistics of difference image 
    typename diffim::ImageStatistics<ImageT>::Ptr kStats(new diffim::ImageStatistics<ImageT>());
    (*kStats).apply(diffIm);
    logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                       "Kernel pass 1 : Kernel Sum = %.3f; Background = %.3f +/- %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                       kSum, background, backgroundError,
                       (*kStats).getMean(),
                       (*kStats).getRms());

    // A second pass with a better variance estimate from first difference image
    bool iterateKernel = _policy.getBool("iterateKernel");
    if (iterateKernel) {
        try {
            try {
                _kFunctor->apply(*_miToConvolvePtr->getImage(), *_miToNotConvolvePtr->getImage(),
                                 *diffIm.getVariance(), _policy);
            } catch (exceptions::Exception& e) {
                throw;
            }
            kernelPtr       = _kFunctor->getKernel();
            kernelErrorPtr  = _kFunctor->getKernelError();
            background      = _kFunctor->getBackground();
            backgroundError = _kFunctor->getBackgroundError();
            
            kSum    = 0.;
            kSum    = kernelPtr->computeImage(kImage, false);
            diffIm  = convolveAndSubtract(*_miToConvolvePtr->getImage(), *_miToNotConvolvePtr, *kernelPtr, background);

            // Reset the image its looking at
            kStats->apply(diffIm);

            logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                               "Kernel pass 2 : Kernel Sum = %.3f; Background = %.3f +/- %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                               kSum, background, backgroundError, kStats->getMean(), kStats->getRms());
            
        } catch (exceptions::Exception& e) {
            logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                               "Exception caught from computePsfMatchingKernelForFootprint, reverting to first solution");
            logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                               e.what());
        }
    }
        
    // Updates for derived class
    _kPtr    = kernelPtr;
    _kErrPtr = kernelErrorPtr;
    _kSum    = kSum;
    _bg      = background;
    _bgErr   = backgroundError;
    _kStats  = kStats;
    // Updates for base class
    setStatus(kStats->evaluateQuality(_policy));
    setBuildStatus(true);

    logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                       "Kernel : Kernel Sum = %.3f; Background = %.3f +/- %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                       _kSum, background, backgroundError, _kStats->getMean(), _kStats->getRms());

    // Return quality of the kernel
    return getStatus();
}

template <typename ImageT>
double SpatialModelKernel<ImageT>::returnCellRating() {
    // Currently, just check the total flux in the template image
    FindCounts<ImageT> counter;
    counter.apply(*_miToConvolvePtr);
    return counter.getCounts();
}

// Explicit instantiations
template class SpatialModelKernel<float>;
template class SpatialModelKernel<double>;

}}} // end of namespace lsst::ip::diffim

