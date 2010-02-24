// -*- lsst-c++ -*-
/**
 * @file KernelCandidate.cc
 *
 * @brief Implementation of KernelCandidate class
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/pex/exceptions/Runtime.h"

#include "lsst/ip/diffim/ImageSubtract.h"
#include "lsst/ip/diffim/KernelCandidate.h"

namespace afwMath        = lsst::afw::math;
namespace afwImage       = lsst::afw::image;
namespace pexLogging     = lsst::pex::logging; 
namespace pexExcept      = lsst::pex::exceptions; 

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
//
// Explicit instantiations
//
typedef float PixelT;
template class KernelCandidate<PixelT>;

}}} // end of namespace lsst::ip::diffim
