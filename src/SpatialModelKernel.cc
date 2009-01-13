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

#include <lsst/pex/policy/Policy.h>
#include <lsst/pex/logging/Trace.h>

#include <lsst/ip/diffim/SpatialModelBase.h>
#include <lsst/ip/diffim/SpatialModelKernel.h>

namespace lsst {
namespace ip {
namespace diffim {
    
template <typename ImageT, typename MaskT>
SpatialModelKernel<ImageT, MaskT>::SpatialModelKernel() :
    lsst::ip::diffim::SpatialModelBase()
{;}

template <typename ImageT, typename MaskT>
SpatialModelKernel<ImageT, MaskT>::SpatialModelKernel(
    lsst::detection::Footprint::PtrType fpPtr,
    MaskedImagePtr miToConvolveParentPtr,
    MaskedImagePtr miToNotConvolveParentPtr,
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> kBasisList,
    lsst::pex::policy::Policy &policy,
    bool build
    ) :
    lsst::ip::diffim::SpatialModelBase(),
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

template <typename ImageT, typename MaskT>
bool SpatialModelKernel<ImageT, MaskT>::buildModel() {

    typedef lsst::afw::image::Image<double>::pixel_accessor pixelAccessor;

    if (this->getBuildStatus() == true) {
        return false;
    }

    // fill in information on position in the image
    vw::BBox2i   fpBBox = this->_fpPtr->getBBox();
    vw::Vector2i fpMin  = fpBBox.min();
    vw::Vector2i fpMax  = fpBBox.max();
    this->setColcNorm(float(fpMin.x() + fpMax.x()) / this->_miToConvolveParentPtr->getCols() - 1.0);
    this->setRowcNorm(float(fpMin.y() + fpMax.y()) / this->_miToConvolveParentPtr->getRows() - 1.0);

    lsst::pex::logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                                  "Footprint = %d,%d -> %d,%d",
                                  fpBBox.min().x(), fpBBox.min().y(),
                                  fpBBox.max().x(), fpBBox.max().y());

    // Fill in information on the actual pixels used
    MaskedImagePtr miToConvolvePtr    = this->_miToConvolveParentPtr->getSubImage(fpBBox);
    MaskedImagePtr miToNotConvolvePtr = this->_miToNotConvolveParentPtr->getSubImage(fpBBox);
    this->_miToConvolvePtr = miToConvolvePtr;
    this->_miToNotConvolvePtr = miToNotConvolvePtr;

    // Estimate of the variance for first kernel pass
    lsst::afw::image::MaskedImage<ImageT, MaskT> varEstimate = 
        lsst::afw::image::MaskedImage<ImageT, MaskT>(this->_miToConvolvePtr->getCols(), 
                                                     this->_miToConvolvePtr->getRows());
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
    } catch (lsst::pex::exceptions::ExceptionStack &e) {
        this->setSdqaStatus(false);
        return false;
    }

    // Compute kernel sum
    double kSum = 0.;
    unsigned int kCols = this->_policy.getInt("kernelCols");
    unsigned int kRows = this->_policy.getInt("kernelRows");
    lsst::afw::image::Image<double> kImage = 
        lsst::afw::image::Image<double>(kCols, kRows);
    kImage = kernelPtr->computeNewImage(kSum, false);

    // Create difference image and calculate associated statistics
    lsst::afw::image::MaskedImage<ImageT, MaskT> diffIm = convolveAndSubtract( *(this->_miToConvolvePtr),
                                                                               *(this->_miToNotConvolvePtr),
                                                                               kernelPtr, 
                                                                               background);
    lsst::ip::diffim::DifferenceImageStatistics<ImageT, MaskT> kStats = 
        lsst::ip::diffim::DifferenceImageStatistics<ImageT, MaskT>(diffIm);

    lsst::pex::logging::TTrace<6>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                                  "Kernel pass 1 : Kernel Sum = %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                                  kSum, 
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
            } catch (lsst::pex::exceptions::ExceptionStack &e) {
                // Exit to first solution
                throw;
            }
            
            kSum    = 0.;
            kImage  = kernelPtr->computeNewImage(kSum, false);
            diffIm  = convolveAndSubtract( *(this->_miToConvolvePtr),
                                           *(this->_miToNotConvolvePtr),
                                           kernelPtr, 
                                           background);

            kStats = lsst::ip::diffim::DifferenceImageStatistics<ImageT, MaskT>(diffIm);
            lsst::pex::logging::TTrace<6>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                                          "Kernel pass 2 : Kernel Sum = %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                                          kSum, 
                                          kStats.getResidualMean(),
                                          kStats.getResidualStd());
            
        } catch (lsst::pex::exceptions::ExceptionStack &e) {
            // Use the first solution
            ;
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
                                  "Kernel : Kernel Sum = %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                                  this->_kSum, 
                                  this->_kStats.getResidualMean(),
                                  this->_kStats.getResidualStd());


    // Return quality of the kernel
    return this->getSdqaStatus();
}

template <typename ImageT, typename MaskT>
double SpatialModelKernel<ImageT, MaskT>::returnSdqaRating() {
    return this->_kStats.getResidualMean();
}

// Explicit instantiations
template class SpatialModelKernel<float, lsst::afw::image::maskPixelType>;
template class SpatialModelKernel<double, lsst::afw::image::maskPixelType>;

}}} // end of namespace lsst::ip::diffim

