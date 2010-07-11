// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
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

/** 
 * @brief Runs Detection on a single image for significant peaks, and checks
 * returned Footprints for Masked pixels.
 *
 * @note Accepts two MaskedImages, one of which is to be convolved to match the
 * other.  The Detection package is run on the image to be convolved
 * (assumed to be higher S/N than the other image).  The subimages
 * associated with each returned Footprint in both images are checked for
 * Masked pixels; Footprints containing Masked pixels are rejected.  The
 * Footprints are grown by an amount specified in the Policy.  The
 * acceptible Footprints are returned in a vector.
 *
 * @return Vector of "clean" Footprints around which Image Subtraction
 * Kernels will be built.
 *
 */
template <typename PixelT>
std::vector<afwDetect::Footprint::Ptr> getCollectionOfFootprintsForPsfMatching(
    lsst::afw::image::MaskedImage<PixelT> const &imageToConvolve,    
    lsst::afw::image::MaskedImage<PixelT> const &imageToNotConvolve, 
    lsst::pex::policy::Policy             const &policy                                       
    ) {
    
    // Parse the Policy
    unsigned int fpNpixMin      = policy.getInt("fpNpixMin");
    unsigned int fpNpixMax      = policy.getInt("fpNpixMax");
    int fpGrowPix               = policy.getInt("fpGrowPix");

    double detThreshold         = policy.getDouble("detThreshold");
    std::string detThresholdType = policy.getString("detThresholdType");

    // New mask plane that tells us which pixels are already in sources
    // Add to both images so mask planes are aligned
    int diffimMaskPlane = imageToConvolve.getMask()->addMaskPlane(diffimStampCandidateStr);
    (void)imageToNotConvolve.getMask()->addMaskPlane(diffimStampCandidateStr);
    afwImage::MaskPixel const diffimBitMask = 
        imageToConvolve.getMask()->getPlaneBitMask(diffimStampCandidateStr);

    // Add in new plane that will tell us which ones are used
    (void)imageToConvolve.getMask()->addMaskPlane(diffimStampUsedStr);
    (void)imageToNotConvolve.getMask()->addMaskPlane(diffimStampUsedStr);

    // List of Footprints
    std::vector<afwDetect::Footprint::Ptr> footprintListIn;
    std::vector<afwDetect::Footprint::Ptr> footprintListOut;

    // Functors to search through the images for masked pixels within candidate footprints
    FindSetBits<afwImage::Mask<afwImage::MaskPixel> > fsb;
 
    int nCleanFp = 0;
    imageToConvolve.getMask()->clearMaskPlane(diffimMaskPlane);
    imageToNotConvolve.getMask()->clearMaskPlane(diffimMaskPlane);
    
    footprintListIn.clear();
    footprintListOut.clear();
    
    // Find detections
    afwDetect::Threshold threshold = 
        afwDetect::createThreshold(detThreshold, detThresholdType);
    afwDetect::FootprintSet<PixelT> footprintSet(
        imageToConvolve, 
        threshold,
        "",
        fpNpixMin);
    
    // Get the associated footprints
    footprintListIn = footprintSet.getFootprints();
    pexLog::TTrace<4>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                      "Found %d total footprints above threshold %.3f",
                      footprintListIn.size(), detThreshold);
    
    // Iterate over footprints, look for "good" ones
    nCleanFp = 0;
    for (std::vector<afwDetect::Footprint::Ptr>::iterator i = footprintListIn.begin(); 
         i != footprintListIn.end(); ++i) {
        // footprint has too many pixels
        if (static_cast<unsigned int>((*i)->getNpix()) > fpNpixMax) {
            pexLog::TTrace<6>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                              "Footprint has too many pix: %d (max =%d)", 
                              (*i)->getNpix(), fpNpixMax);
            continue;
        } 
        
        pexLog::TTrace<8>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                          "Footprint in : %d,%d -> %d,%d",
                          (*i)->getBBox().getX0(), (*i)->getBBox().getY0(), 
                          (*i)->getBBox().getX1(), (*i)->getBBox().getY1());
        
        pexLog::TTrace<8>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                          "Grow by : %d pixels", fpGrowPix);
        
        /* Grow the footprint
           flag true  = isotropic grow   = slow
           flag false = 'manhattan grow' = fast
           
           The manhattan masks are rotated 45 degree w.r.t. the coordinate
           system.  They intersect the vertices of the rectangle that would
           connect pixels (X0,Y0) (X1,Y0), (X0,Y1), (X1,Y1).
           
           The isotropic masks do take considerably longer to grow and are
           basically elliptical.  X0, X1, Y0, Y1 delimit the extent of the
           ellipse.
           
           In both cases, since the masks aren't rectangles oriented with
           the image coordinate system, when we DO extract such rectangles
           as subimages for kernel fitting, some corner pixels can be found
           in multiple subimages.
           
        */
        afwDetect::Footprint::Ptr fpGrow = 
            afwDetect::growFootprint(*i, fpGrowPix, false);

        /* Operate on the BBox when grabbing sub images from the main image */
        afwImage::BBox fpBBox = fpGrow->getBBox();
        
        pexLog::TTrace<7>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                          "Footprint out : %d,%d -> %d,%d (center %d,%d)",
                          fpBBox.getX0(), fpBBox.getY0(),
                          fpBBox.getX1(), fpBBox.getY1(),
                          int(0.5 * (fpBBox.getX0() + fpBBox.getX1())),
                          int(0.5 * (fpBBox.getY0() + fpBBox.getY1())));
        
        
        /* 
           Note we need to translate to pixel coordinates here since thats how
           BBoxes address (sub)images.  This does not affect the Footprint.  
        */
        fpBBox.shift(-imageToConvolve.getX0(), -imageToConvolve.getY0());

        pexLog::TTrace<7>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                          "Footprint shifted : %d,%d -> %d,%d (center %d,%d)",
                          fpBBox.getX0(), fpBBox.getY0(),
                          fpBBox.getX1(), fpBBox.getY1(),
                          int(0.5 * (fpBBox.getX0() + fpBBox.getX1())),
                          int(0.5 * (fpBBox.getY0() + fpBBox.getY1())));

        // Ignore if its too close to the edge of the image 
        bool belowOriginX = fpBBox.getX0() < 0;
        bool belowOriginY = fpBBox.getY0() < 0;
        bool offImageX    = fpBBox.getX1() > imageToConvolve.getWidth();
        bool offImageY    = fpBBox.getY1() > imageToConvolve.getHeight();
        if (belowOriginX || belowOriginY || offImageX || offImageY) {
            pexLog::TTrace<6>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                              "Footprint grown off image"); 
            continue;
        }
        
        // Grab a subimage; report any exception
        try {
            afwImage::MaskedImage<PixelT> subImageToConvolve(imageToConvolve, fpBBox);
            afwImage::MaskedImage<PixelT> subImageToNotConvolve(imageToNotConvolve, fpBBox);
            
            // Search for any masked pixels within the footprint
            fsb.apply(*(subImageToConvolve.getMask()));
            if (fsb.getBits() > 0) {
                pexLog::TTrace<6>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                                  "Footprint has masked pix (val=%d) in image to convolve", 
                                  fsb.getBits()); 
                continue;
            }
            
            fsb.apply(*(subImageToNotConvolve.getMask()));
            if (fsb.getBits() > 0) {
                pexLog::TTrace<6>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                                  "Footprint has masked pix (val=%d) in image not to convolve", 
                                  fsb.getBits());
                continue;
            }
            
        } catch (pexExcept::Exception& e) {
            pexLog::TTrace<6>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching",
                              "Exception caught extracting Footprint");
            pexLog::TTrace<7>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching",
                              e.what());
            continue;
        }

        // If we get this far, we have a clean footprint
        footprintListOut.push_back(fpGrow);
        (void)afwDetect::setMaskFromFootprint(&(*imageToConvolve.getMask()), *fpGrow, diffimBitMask);
        (void)afwDetect::setMaskFromFootprint(&(*imageToNotConvolve.getMask()), *fpGrow, diffimBitMask);
        nCleanFp += 1;
    }
    
    imageToConvolve.getMask()->clearMaskPlane(diffimMaskPlane);
    imageToNotConvolve.getMask()->clearMaskPlane(diffimMaskPlane);
    
    if (footprintListOut.size() == 0) {
      throw LSST_EXCEPT(pexExcept::Exception, 
                        "Unable to find any footprints for Psf matching");
    }

    pexLog::TTrace<1>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                      "Found %d clean footprints above threshold %.3f",
                      footprintListOut.size(), detThreshold);
    
    return footprintListOut;
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

/* */


template
std::vector<lsst::afw::detection::Footprint::Ptr> getCollectionOfFootprintsForPsfMatching(
    lsst::afw::image::MaskedImage<float> const &,
    lsst::afw::image::MaskedImage<float> const &,
    lsst::pex::policy::Policy const &);

template
std::vector<lsst::afw::detection::Footprint::Ptr> getCollectionOfFootprintsForPsfMatching(
    lsst::afw::image::MaskedImage<double> const &,
    lsst::afw::image::MaskedImage<double> const &,
    pexPolicy::Policy  const &);

}}} // end of namespace lsst::ip::diffim
