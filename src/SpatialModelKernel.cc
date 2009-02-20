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

#include <lsst/ip/diffim/SpatialModelBase.h>
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
    
/*
template <typename ImageT>
SpatialModelKernel<ImageT>::SpatialModelKernel() :
    lsst::ip::diffim::SpatialModelBase<ImageT>(),
    _fpPtr(lsst::afw::detection::Footprint::Ptr()),
    _miToConvolvePtr(MaskedImagePtr()),
    _miToNotConvolvePtr(MaskedImagePtr()),
    _kFunctor(boost::shared_ptr<PsfMatchingFunctor<ImageT> >()),
    _policy(),
    _kPtr(boost::shared_ptr<lsst::afw::math::Kernel>()),
    _kErrPtr(boost::shared_ptr<lsst::afw::math::Kernel>()),
    _kSum(0.),
    _bg(0.),
    _bgErr(0.),
    _kStats(boost::shared_ptr<ImageStatistics<lsst::afw::image::MaskedImage<ImageT> > >())

{;}
*/

template <typename ImageT>
SpatialModelKernel<ImageT>::SpatialModelKernel(
    detection::Footprint::Ptr const &fpPtr,
    MaskedImagePtr const &miToConvolvePtr,
    MaskedImagePtr const &miToNotConvolvePtr,
    boost::shared_ptr<PsfMatchingFunctor<ImageT> > const &kFunctor,
    policy::Policy const &policy,
    bool build
    ) :
    diffim::SpatialModelBase<ImageT>(),
    _fpPtr(fpPtr),
    _miToConvolvePtr(miToConvolvePtr),
    _miToNotConvolvePtr(miToNotConvolvePtr),
    _kFunctor(kFunctor),
    _policy(policy),
    _kPtr(boost::shared_ptr<math::Kernel>()),
    _kErrPtr(boost::shared_ptr<math::Kernel>()),
    _kSum(0.),
    _bg(0.),
    _bgErr(0.),
    _kStats(boost::shared_ptr<ImageStatistics<image::MaskedImage<ImageT> > >())
{
    if (build == true) {
        this->buildModel();
    }
}

template <typename ImageT>
bool SpatialModelKernel<ImageT>::buildModel() {

    if (this->getBuildStatus() == true) {
        return false;
    }

    // fill in information on position in the image
    image::BBox fpBBox = this->_fpPtr->getBBox();

    // NOTE : since we can't remap pixel range to go from -1 to 1 in convolve(),
    // we have to use the actual pixel value here.  Not optimal.
    // this->setColc(float(fpBBox.getX0() + fpBBox.getX1()) / this->_miToConvolveParentPtr->getWidth() - 1.0);
    // this->setRowc(float(fpBBox.getY0() + fpBBox.getY1()) / this->_miToConvolveParentPtr->getHeight() - 1.0);
    this->setColc(0.5 * float(fpBBox.getX0() + fpBBox.getX1()));
    this->setRowc(0.5 * float(fpBBox.getY0() + fpBBox.getY1()));

    logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                       "Footprint = %d,%d -> %d,%d",
                       fpBBox.getX0(), fpBBox.getY0(),
                       fpBBox.getX1(), fpBBox.getY1());

    // Estimate of the variance for first kernel pass
    // True argument is for a deep copy, so -= does not modify the original pixels
    image::MaskedImage<ImageT> varEstimate = 
        image::MaskedImage<ImageT>(*(this->_miToNotConvolvePtr), true);
    varEstimate -= *(this->_miToConvolvePtr);
    

    try {
        this->_kFunctor->apply(*(this->_miToConvolvePtr), 
                               *(this->_miToNotConvolvePtr),
                               *(varEstimate.getVariance()),
                               this->_policy);
    } catch (exceptions::Exception& e) {
        this->setSdqaStatus(false);
        logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                           "Exception caught from computePsfMatchingKernelForFootprint"); 
        logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                           e.what());
        return false;
    }
    boost::shared_ptr<math::Kernel> kernelPtr       = this->_kFunctor->getKernel();
    boost::shared_ptr<math::Kernel> kernelErrorPtr  = this->_kFunctor->getKernelError();
    double                          background      = this->_kFunctor->getBackground();
    double                          backgroundError = this->_kFunctor->getBackgroundError();

    // Compute kernel sum
    double kSum = 0.;
    unsigned int kCols = this->_policy.getInt("kernelCols");
    unsigned int kRows = this->_policy.getInt("kernelRows");
    image::Image<double> kImage(kCols, kRows);
    kSum = kernelPtr->computeImage(kImage, false);

    // Create difference image and calculate associated statistics
    image::MaskedImage<ImageT> diffIm = convolveAndSubtract( *(this->_miToConvolvePtr),
                                                             *(this->_miToNotConvolvePtr),
                                                             *(kernelPtr), 
                                                             background);

    // Find statistics of difference image 
    boost::shared_ptr<diffim::ImageStatistics<image::MaskedImage<ImageT> > > kStats = 
        boost::shared_ptr<diffim::ImageStatistics<image::MaskedImage<ImageT> > > (
            new diffim::ImageStatistics<image::MaskedImage<ImageT> >(diffIm)
            );
    detection::Footprint fp(
        image::BBox(image::PointI(0,0), 
                    diffIm.getWidth(),
                    diffIm.getHeight()
            )
        );
    (*kStats).apply(fp);
    logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                       "Kernel pass 1 : Kernel Sum = %.3f; Background = %.3f +/- %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                       kSum, background, backgroundError,
                       (*kStats).getMean(),
                       sqrt((*kStats).getVariance()));

    // A second pass with a better variance estimate from first difference image
    bool iterateKernel = this->_policy.getBool("iterateKernel");
    if (iterateKernel) {
        try {
            try {
                this->_kFunctor->apply(*(this->_miToConvolvePtr), 
                                       *(this->_miToNotConvolvePtr),
                                       *(diffIm.getVariance()),
                                       this->_policy);
            } catch (exceptions::Exception& e) {
                throw;
            }
            kernelPtr       = this->_kFunctor->getKernel();
            kernelErrorPtr  = this->_kFunctor->getKernelError();
            background      = this->_kFunctor->getBackground();
            backgroundError = this->_kFunctor->getBackgroundError();
            
            kSum    = 0.;
            kSum    = kernelPtr->computeImage(kImage, false);
            diffIm  = convolveAndSubtract( *(this->_miToConvolvePtr),
                                           *(this->_miToNotConvolvePtr),
                                           *(kernelPtr), 
                                           background);

            // Reset image
            kStats = 
                boost::shared_ptr<diffim::ImageStatistics<image::MaskedImage<ImageT> > > (
                    new diffim::ImageStatistics<image::MaskedImage<ImageT> >(diffIm)
                    );
            (*kStats).apply(fp);

            logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                               "Kernel pass 2 : Kernel Sum = %.3f; Background = %.3f +/- %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                               kSum, background, backgroundError,
                               (*kStats).getMean(),
                               sqrt((*kStats).getVariance()));
            
        } catch (exceptions::Exception& e) {
            logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                               "Exception caught from computePsfMatchingKernelForFootprint, reverting to first solution");
            logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                               e.what());
        }
    }
        
    // Updates for derived class
    this->_kPtr    = kernelPtr;
    this->_kErrPtr = kernelErrorPtr;
    this->_kSum    = kSum;
    this->_bg      = background;
    this->_bgErr   = backgroundError;
    this->_kStats  = kStats;
    // Updates for base class
    this->setSdqaStatus((*kStats).evaluateQuality(this->_policy));
    this->setBuildStatus(true);

    logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                       "Kernel : Kernel Sum = %.3f; Background = %.3f +/- %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                       this->_kSum, background, backgroundError,
                       this->_kStats->getMean(),
                       sqrt(this->_kStats->getVariance()));


    // Return quality of the kernel
    return this->getSdqaStatus();
}

template <typename ImageT>
double SpatialModelKernel<ImageT>::returnSdqaRating(lsst::pex::policy::Policy &policy) {
    return this->_kStats->evaluateQuality(policy);
}

// Explicit instantiations
template class SpatialModelKernel<float>;
template class SpatialModelKernel<double>;

}}} // end of namespace lsst::ip::diffim

