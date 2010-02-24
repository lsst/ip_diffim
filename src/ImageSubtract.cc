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
#include <limits>

#include "boost/timer.hpp" 

#include "Eigen/Core"

// NOTE -  trace statements >= 6 can ENTIRELY kill the run time
// #define LSST_MAX_TRACE 5

#include "lsst/ip/diffim/ImageSubtract.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/pex/exceptions/Exception.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/pex/logging/Log.h"

namespace pexExcept  = lsst::pex::exceptions; 
namespace pexLog     = lsst::pex::logging; 
namespace pexPolicy  = lsst::pex::policy; 
namespace afwImage   = lsst::afw::image;
namespace afwMath    = lsst::afw::math;
namespace afwDetect  = lsst::afw::detection;

namespace lsst { 
namespace ip { 
namespace diffim {

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

// Explicit instantiations
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

}}} // end of namespace lsst::ip::diffim
