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

namespace pexExcept = lsst::pex::exceptions; 

namespace lsst {
namespace ip {
namespace diffim {
    
template <typename ImageT>
SpatialModelKernel<ImageT>::SpatialModelKernel() :
    lsst::ip::diffim::SpatialModelBase<ImageT>()
{;}

template <typename ImageT>
SpatialModelKernel<ImageT>::SpatialModelKernel(
    lsst::afw::detection::Footprint::Ptr fpPtr,
    MaskedImagePtr miToConvolveParentPtr,
    MaskedImagePtr miToNotConvolveParentPtr,
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> kBasisList,
    lsst::pex::policy::Policy &policy,
    bool build
    ) :
    lsst::ip::diffim::SpatialModelBase<ImageT>(),
    _miToConvolveParentPtr(miToConvolveParentPtr),
    _miToNotConvolveParentPtr(miToNotConvolveParentPtr),
    _kBasisList(kBasisList),
    _policy(policy),
    _fpPtr(fpPtr)
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
    lsst::afw::image::BBox fpBBox = this->_fpPtr->getBBox();

    // NOTE : since we can't remap pixel range to go from -1 to 1 in convolve(),
    // we have to use the actual pixel value here.  Not optimal.

    // this->setColc(float(fpBBox.getX0() + fpBBox.getX1()) / this->_miToConvolveParentPtr->getWidth() - 1.0);
    // this->setRowc(float(fpBBox.getY0() + fpBBox.getY1()) / this->_miToConvolveParentPtr->getHeight() - 1.0);

    this->setColc(0.5 * float(fpBBox.getX0() + fpBBox.getX1()));
    this->setRowc(0.5 * float(fpBBox.getY0() + fpBBox.getY1()));

    lsst::pex::logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                                  "Footprint = %d,%d -> %d,%d",
                                  fpBBox.getX0(), fpBBox.getY0(),
                                  fpBBox.getX1(), fpBBox.getY1());

    // Fill in information on the actual pixels used
    MaskedImagePtr miToConvolvePtr    = MaskedImagePtr ( 
        new lsst::afw::image::MaskedImage<ImageT>(*(this->_miToConvolveParentPtr), fpBBox)
        );

    MaskedImagePtr miToNotConvolvePtr = MaskedImagePtr ( 
        new lsst::afw::image::MaskedImage<ImageT>(*(this->_miToNotConvolveParentPtr), fpBBox)
        );
    this->_miToConvolvePtr = miToConvolvePtr;
    this->_miToNotConvolvePtr = miToNotConvolvePtr;

    // Estimate of the variance for first kernel pass
    lsst::afw::image::MaskedImage<ImageT> varEstimate = 
        lsst::afw::image::MaskedImage<ImageT>(this->_miToConvolvePtr->getDimensions());

    varEstimate += *(this->_miToNotConvolvePtr);
    varEstimate -= *(this->_miToConvolvePtr);
    
    boost::shared_ptr<lsst::afw::math::Kernel> kernelPtr;
    boost::shared_ptr<lsst::afw::math::Kernel> kernelErrorPtr;
    double                                     background;
    double                                     backgroundError;
    try {
        computePsfMatchingKernelForFootprint(*(this->_miToConvolvePtr), 
                                             *(this->_miToNotConvolvePtr), 
                                             varEstimate, 
                                             this->_kBasisList, 
                                             this->_policy, 
                                             kernelPtr, kernelErrorPtr,
                                             background, backgroundError);
    } catch (pexExcept::Exception& e) {
        this->setSdqaStatus(false);
        lsst::pex::logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                                      "Exception caught from computePsfMatchingKernelForFootprint"); 
        lsst::pex::logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                                      e.what());
        return false;
    }

    // Compute kernel sum
    double kSum = 0.;
    unsigned int kCols = this->_policy.getInt("kernelCols");
    unsigned int kRows = this->_policy.getInt("kernelRows");
    lsst::afw::image::Image<double> kImage(kCols, kRows);
    kSum = kernelPtr->computeImage(kImage, false);

    // Create difference image and calculate associated statistics
    lsst::afw::image::MaskedImage<ImageT> diffIm = convolveAndSubtract( *(this->_miToConvolvePtr),
                                                                               *(this->_miToNotConvolvePtr),
                                                                               *(kernelPtr), 
                                                                               background);
    lsst::ip::diffim::DifferenceImageStatistics<ImageT> kStats = 
        lsst::ip::diffim::DifferenceImageStatistics<ImageT>(diffIm);

    lsst::pex::logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                                  "Kernel pass 1 : Kernel Sum = %.3f; Background = %.3f +/- %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                                  kSum, background, backgroundError,
                                  kStats.getResidualMean(),
                                  kStats.getResidualStd());

    // A second pass with a better variance estimate from first difference image
    bool iterateKernel = this->_policy.getBool("iterateKernel");
    if (iterateKernel) {
        try {
            try {
                computePsfMatchingKernelForFootprint(*(this->_miToConvolvePtr), 
                                                     *(this->_miToNotConvolvePtr), 
                                                     diffIm,
                                                     this->_kBasisList, 
                                                     this->_policy, 
                                                     kernelPtr, kernelErrorPtr,
                                                     background, backgroundError);
            } catch (pexExcept::Exception& e) {
                throw;
            }
            
            kSum    = 0.;
            kSum = kernelPtr->computeImage(kImage, false);
            diffIm  = convolveAndSubtract( *(this->_miToConvolvePtr),
                                           *(this->_miToNotConvolvePtr),
                                           *(kernelPtr), 
                                           background);

            kStats = lsst::ip::diffim::DifferenceImageStatistics<ImageT>(diffIm);
            lsst::pex::logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                                          "Kernel pass 2 : Kernel Sum = %.3f; Background = %.3f +/- %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                                          kSum, background, backgroundError,
                                          kStats.getResidualMean(),
                                          kStats.getResidualStd());
            
        } catch (pexExcept::Exception& e) {
            lsst::pex::logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                                          "Exception caught from computePsfMatchingKernelForFootprint, reverting to first solution");
            lsst::pex::logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
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
    this->setSdqaStatus(kStats.evaluateQuality(this->_policy));
    this->setBuildStatus(true);

    lsst::pex::logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                                  "Kernel : Kernel Sum = %.3f; Background = %.3f +/- %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                                  this->_kSum, background, backgroundError,
                                  this->_kStats.getResidualMean(),
                                  this->_kStats.getResidualStd());


    // Return quality of the kernel
    return this->getSdqaStatus();
}

template <typename ImageT>
double SpatialModelKernel<ImageT>::returnSdqaRating() {
    return this->_kStats.getResidualMean();
}

// Explicit instantiations
template class SpatialModelKernel<float>;
template class SpatialModelKernel<double>;

}}} // end of namespace lsst::ip::diffim

