// -*- lsst-c++ -*-
/**
 * @file ImageSubtract.h
 *
 * @brief Image Subtraction helper functions
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_IMAGESUBTRACT_H
#define LSST_IP_DIFFIM_IMAGESUBTRACT_H

#include <vector>
#include <string>

#include "Eigen/Core"

#include "boost/shared_ptr.hpp"

#include "lsst/pex/policy/Policy.h"
#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Footprint.h"

namespace lsst { 
namespace ip { 
namespace diffim {

    
    /** Mask plane definitions */
    std::string const diffimStampCandidateStr = "DIFFIM_STAMP_CANDIDATE";
    std::string const diffimStampUsedStr      = "DIFFIM_STAMP_USED";
    
    /**
     * @brief Uses a functor to accumulate Mask bits
     *
     * @note Search through a footprint for any set mask fits.
     * 
     * @note May need to modify this as our mask planes evolve to include
     * non-bad mask information
     *
     * @code
        FindSetBits<image::Mask<image::MaskPixel> > count(mask); 
        count.reset(); 
        count.apply(footprint); 
        nSet = count.getBits();
     * @endcode
     * 
     * @ingroup ip_diffim
     *
     */
    template <typename MaskT>
    class FindSetBits : public lsst::afw::detection::FootprintFunctor<MaskT> {
    public:
        FindSetBits(MaskT const& mask) : 
            lsst::afw::detection::FootprintFunctor<MaskT>(mask), _bits(0) {;}
        
        void operator()(typename MaskT::xy_locator loc, ///< locator pointing at the pixel
                        int x,                          ///< column-position of pixel
                        int y                           ///< row-position of pixel
            ) {
            _bits |= *loc;
        }
        
        // Return the bits set
        typename MaskT::Pixel getBits() const { return _bits; }
        
        // Clear the accumulator
        void reset() { _bits = 0; }
        
    private:
        typename MaskT::Pixel _bits;
    };

    /**
     * @brief Class to calculate difference image statistics
     *
     * @note Find mean and unbiased variance of pixel residuals in units of
     * sqrt(variance)
     * 
     * @ingroup ip_diffim
     */
    template <typename PixelT>
    class ImageStatistics {
    public:
        typedef boost::shared_ptr<ImageStatistics> Ptr;
        typedef typename lsst::afw::image::MaskedImage<PixelT>::x_iterator x_iterator;

        ImageStatistics() : 
            _xsum(0.), _x2sum(0.), _npix(0) {} ;
        virtual ~ImageStatistics() {} ;

        // Clear the accumulators
        void reset() { _xsum = _x2sum = 0.; _npix = 0;}

        // Work your magic
        void apply(lsst::afw::image::MaskedImage<PixelT> const& image) {
            reset();
            for (int y = 0; y != image.getHeight(); ++y) {
                for (x_iterator ptr = image.row_begin(y), end = image.row_end(y); ptr != end; ++ptr) {
		    if ((*ptr).mask() == 0) {
                        double const ivar = 1. / (*ptr).variance();
                        _xsum  += (*ptr).image() * sqrt(ivar);
                        _x2sum += (*ptr).image() * (*ptr).image() * ivar;
                        _npix  += 1;
                    }
                }
            }
        }
        
        // Mean of distribution
        double getMean() const { 
            return (_npix > 0) ? _xsum/_npix : std::numeric_limits<double>::quiet_NaN(); 
        }
        // Variance of distribution 
        double getVariance() const { 
            return (_npix > 1) ? (_x2sum/_npix - _xsum/_npix * _xsum/_npix) * _npix/(_npix-1.) : std::numeric_limits<double>::quiet_NaN(); 
        }
        // RMS
        double getRms() const { 
            return sqrt(getVariance());
        }
        // Return the number of good pixels
        int getNpix() const { return _npix; }

        // Return Sdqa rating
        bool evaluateQuality(lsst::pex::policy::Policy const& policy) {
            if ( fabs(getMean())     > policy.getDouble("maximumFootprintResidualMean") ) return false;
            if ( getRms()            > policy.getDouble("maximumFootprintResidualStd")  ) return false;
            return true;
        }           
        
    private:
        double _xsum;
        double _x2sum;
        int    _npix;
    };


    /**
     * @brief Execute fundamental task of convolving template and subtracting it from science image
     * 
     * @note This version accepts a MaskedImage for the template
     * 
     * @param imageToConvolve  MaskedImage to apply convolutionKernel to
     * @param imageToNotConvolve  MaskedImage from which convolved imageToConvolve is subtracted 
     * @param convolutionKernel  Kernel to apply to imageToConvolve
     * @param background  Background scalar or function to subtract after convolution
     * @param invert  Invert the output difference image
     * 
     * @ingroup ip_diffim
     */
    template <typename PixelT, typename BackgroundT>
    lsst::afw::image::MaskedImage<PixelT> convolveAndSubtract(
        lsst::afw::image::MaskedImage<PixelT> const& imageToConvolve,
        lsst::afw::image::MaskedImage<PixelT> const& imageToNotConvolve,
        lsst::afw::math::Kernel const& convolutionKernel,
        BackgroundT background,
        bool invert=true
        );

    /**
     * @brief Execute fundamental task of convolving template and subtracting it from science image
     * 
     * @note This version accepts an Image for the template, and is thus faster during convolution
     * 
     * @param imageToConvolve  Image to apply convolutionKernel to
     * @param imageToNotConvolve  MaskedImage from which convolved imageToConvolve is subtracted 
     * @param convolutionKernel  Kernel to apply to imageToConvolve
     * @param background  Background scalar or function to subtract after convolution
     * @param invert  Invert the output difference image
     * 
     * @ingroup ip_diffim
     */
    template <typename PixelT, typename BackgroundT>
    lsst::afw::image::MaskedImage<PixelT> convolveAndSubtract(
        lsst::afw::image::Image<PixelT> const& imageToConvolve,
        lsst::afw::image::MaskedImage<PixelT> const& imageToNotConvolve,
        lsst::afw::math::Kernel const& convolutionKernel,
        BackgroundT background,
        bool invert=true
        );

    /**
     * @brief Search through images for Footprints with no masked pixels
     *
     * @note Runs detection on the template; searches through both images for masked pixels
     *
     * @param imageToConvolve  MaskedImage that will be convolved with kernel; detection is run on this image
     * @param imageToNotConvolve  MaskedImage to subtract convolved template from
     * @param policy  Policy for operations; in particular object detection
     *
     * @ingroup ip_diffim
     */    
    template <typename PixelT>
    std::vector<lsst::afw::detection::Footprint::Ptr> getCollectionOfFootprintsForPsfMatching(
        lsst::afw::image::MaskedImage<PixelT> const& imageToConvolve,
        lsst::afw::image::MaskedImage<PixelT> const& imageToNotConvolve,
        lsst::pex::policy::Policy             const& policy
        );

    /**
     * @brief Turns a 2-d Image into a 2-d Eigen Matrix
     *
     * @param img  Image whose pixel values are read into an Eigen::MatrixXd
     *
     * @ingroup ip_diffim
     */
    template <typename PixelT>
    Eigen::MatrixXd imageToEigenMatrix(
        lsst::afw::image::Image<PixelT> const& img
        );
    
    /**
     * @brief Helper method to add a Function to an Image
     *
     * @param image  Image to be modified
     * @param function  Funtion that is evaluated at all pixel values and added to image
     *
     * @ingroup ip_diffim
     */
    template <typename PixelT>
    void addToImage(lsst::afw::image::Image<PixelT> &image,
                    lsst::afw::math::Function2<double> const &function);

    /**
     * @brief Helper method to add a double to an Image
     *
     * @param image  Image to be modified
     * @param value  Value to be added to image
     *
     * @ingroup ip_diffim
     */
    template <typename PixelT>
    void addToImage(lsst::afw::image::Image<PixelT> &image,
                    double value);

}}} // end of namespace lsst::ip::diffim

#endif



