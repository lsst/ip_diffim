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

    
    /**
     * @brief Class to accumulate Mask bits
     *
     * @note Search through a Mask for any set bits.
     * 
     * @ingroup ip_diffim
     */
    template <typename MaskT>
    class FindSetBits {
    public:
        typedef typename MaskT::x_iterator x_iterator;

        FindSetBits() : 
            _bits(0) {;}
        virtual ~FindSetBits() {} ;

        // Clear the accumulators
        void reset() { _bits = 0;}

        // Return the bits set
        typename MaskT::Pixel getBits() const { return _bits; }

        // Work your magic
        void apply(MaskT const& mask) {
            reset();
            for (int y = 0; y != mask.getHeight(); ++y) {
                for (x_iterator ptr = mask.row_begin(y), end = mask.row_end(y); ptr != end; ++ptr) {
                    _bits |= (*ptr);
                }
            }
        }

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
            _xsum(0.), _x2sum(0.), _npix(0) {
            lsst::afw::image::MaskPixel const edgeMask   = 
                lsst::afw::image::Mask<lsst::afw::image::MaskPixel>::getPlaneBitMask("EDGE");
            lsst::afw::image::MaskPixel const crMask     = 
                lsst::afw::image::Mask<lsst::afw::image::MaskPixel>::getPlaneBitMask("CR");
            lsst::afw::image::MaskPixel const satMask    = 
                lsst::afw::image::Mask<lsst::afw::image::MaskPixel>::getPlaneBitMask("SAT");
            lsst::afw::image::MaskPixel const badMask    = 
                lsst::afw::image::Mask<lsst::afw::image::MaskPixel>::getPlaneBitMask("BAD");
            lsst::afw::image::MaskPixel const interpMask = 
                lsst::afw::image::Mask<lsst::afw::image::MaskPixel>::getPlaneBitMask("INTRP");
            _bpMask = edgeMask | crMask | satMask | badMask | interpMask;
        } ;
        virtual ~ImageStatistics() {} ;

        // Clear the accumulators
        void reset() { _xsum = _x2sum = 0.; _npix = 0;}

        // Work your magic
        void apply(lsst::afw::image::MaskedImage<PixelT> const& image) {
            reset();
            for (int y = 0; y != image.getHeight(); ++y) {
                for (x_iterator ptr = image.row_begin(y), end = image.row_end(y); ptr != end; ++ptr) {
                    if (!((*ptr).mask() & _bpMask)) {
                        double const ivar = 1. / (*ptr).variance();
                        _xsum  += (*ptr).image() * sqrt(ivar);
                        _x2sum += (*ptr).image() * (*ptr).image() * ivar;
                        _npix  += 1;
                    }
                }
            }
        }

        void apply(lsst::afw::image::MaskedImage<PixelT> const& image, int core) {
            reset();
            int y0 = std::max(0, image.getHeight()/2 - core);
            int y1 = std::min(image.getHeight(), image.getHeight()/2 + core);
            int x0 = std::max(0, image.getWidth()/2 - core);
            int x1 = std::min(image.getWidth(), image.getWidth()/2 + core);

            for (int y = y0; y != y1; ++y) {
                for (x_iterator ptr = image.x_at(x0, y), end = image.x_at(x1, y); 
                     ptr != end; ++ptr) {
                    if (!((*ptr).mask() & _bpMask)) {
                        double const ivar = 1. / (*ptr).variance();
                        _xsum  += (*ptr).image() * sqrt(ivar);
                        _x2sum += (*ptr).image() * (*ptr).image() * ivar;
                        _npix  += 1;
                    }
                }
            }
        }

        void setBpMask(lsst::afw::image::MaskPixel bpMask) {_bpMask = bpMask;}
        lsst::afw::image::MaskPixel getBpMask() {return _bpMask;}

        // Mean of distribution
        double getMean() const { 
            return (_npix > 0) ? _xsum/_npix : std::numeric_limits<double>::quiet_NaN(); 
        }
        // Variance of distribution 
        double getVariance() const { 
            return (_npix > 1) ? (_x2sum/_npix - _xsum/_npix * _xsum/_npix) * _npix/(_npix-1.) : 
                std::numeric_limits<double>::quiet_NaN(); 
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
        lsst::afw::image::MaskPixel _bpMask;
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
    
}}} // end of namespace lsst::ip::diffim

#endif



