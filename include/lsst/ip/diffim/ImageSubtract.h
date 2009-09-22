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

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>

#include <lsst/pex/policy/Policy.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/math/KernelFunctions.h>
#include <lsst/afw/math/Function.h>
#include <lsst/afw/math/SpatialCell.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/detection/Footprint.h>

namespace lsst {
namespace ip {
namespace diffim {

    
    template<typename ImageT>
    class SpatialCellKernelCandidate : public lsst::afw::math::SpatialCellImageCandidate<ImageT> {
    public:
        typedef boost::shared_ptr<SpatialCellKernelCandidate<ImageT> > Ptr;
        typedef boost::shared_ptr<const SpatialCellKernelCandidate<ImageT> > ConstPtr;

        SpatialCellKernelCandidate(float const xCenter, ///< The object's column-centre
                                   float const yCenter  ///< The object's row-centre
            ) : lsst::afw::math::SpatialCellImageCandidate<ImageT>(xCenter, yCenter),
                _kernel(typename lsst::afw::math::Kernel::PtrT()){
        }

        /// Return the Candidate's Kernel
        virtual typename lsst::afw::math::Kernel::PtrT getKernel() const = 0;
        
    protected:
        typename lsst::afw::math::Kernel::PtrT mutable _kernel; 
    private:
    };
        

    /** Mask plane definitions */
    std::string const diffimStampCandidateStr = "DIFFIM_STAMP_CANDIDATE";
    std::string const diffimStampUsedStr      = "DIFFIM_STAMP_USED";
    
    /** Uses a functor to accumulate Mask bits
     *
     * @ingroup diffim
     *
     * @note Search through a footprint for any set mask fits.
     * 
     * @note May need to modify this as our mask planes evolve to include
     * non-bad mask information
     *
     * Example usage : 
     *  FindSetBits<image::Mask<image::MaskPixel> > count(mask); 
     *  count.reset(); 
     *  count.apply(footprint); 
     *  nSet = count.getBits();
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

    /** Uses a functor to sum over the MaskedImage pixels
     *
     * @ingroup diffim
     *
     * @note Count the total flux within the image, excluding masked pixels
     * 
     * @note Still needs a background model to correct for
     *
     */
    template <typename ImageT>
    class FindCounts {
    public:
        typedef typename lsst::afw::image::MaskedImage<ImageT>::x_iterator x_iterator;
        FindCounts() : 
            _counts(0.) {} ;
        virtual ~FindCounts() {};

        // Clear the accumulator
        void reset() { _counts = 0.; }

        // Count pixels
        void apply(lsst::afw::image::MaskedImage<ImageT> const& image) {
            reset();
            for (int y = 0; y != image.getHeight(); ++y) {
                for (x_iterator ptr = image.row_begin(y), end = image.row_end(y); ptr != end; ++ptr) {
                    if ((*ptr).mask() == 0) {
                        _counts += (*ptr).image();
                    }
                }
            }
        }

        // Return the total counts
        double getCounts() const { return _counts; }
        
    private:
        double _counts;
    };

    /** Uses a functor to calculate difference image statistics
     *
     * @ingroup diffim
     *
     * @note Looks like this is (almost) implemented in lsst/afw/math/Statistics.h
     * 
     * @note Find mean and unbiased variance of pixel residuals in units of
     * sqrt(variance)
     * 
     */
    template <typename ImageT>
    class ImageStatistics {
    public:
        typedef boost::shared_ptr<ImageStatistics> Ptr;
        typedef typename lsst::afw::image::MaskedImage<ImageT>::x_iterator x_iterator;

        ImageStatistics() : 
            _xsum(0.), _x2sum(0.), _npix(0) {} ;
        virtual ~ImageStatistics() {} ;

        // Clear the accumulators
        void reset() { _xsum = _x2sum = 0.; _npix = 0;}

        // Work your magic
        void apply(lsst::afw::image::MaskedImage<ImageT> const& image) {
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


    /* Build a set of Delta Function basis kernels
     * 
     * @note Total number of basis functions is width*height
     * 
     * @param width  Width of basis set
     * @param height Height of basis set
     */    
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> generateDeltaFunctionKernelSet(
        unsigned int width,
        unsigned int height
        );

    /* Build a regularization matrix for Delta function kernels
     * 
     * @param width            Width of basis set you want to regularize
     * @param height           Height of basis set you want to regularize
     * @param order            Which derivative you expect to be smooth (derivative order+1 is penalized) 
     * @param boundary_style   0 = unwrapped, 1 = wrapped, 2 = order-tappered ('order' is highest used) 
     * @param difference_style 0 = forward, 1 = central
     * @param printB           debugging
     */    
    Eigen::MatrixXd generateFiniteDifferenceRegularization(
        unsigned int width,
        unsigned int height,
        unsigned int order,
	unsigned int boundary_style = 1, 
	unsigned int difference_style = 0,
	bool printB=false
        );

    /** Build a set of Alard/Lupton basis kernels
     *
     * @note Should consider implementing as SeparableKernels for additional speed,
     * but this will make the normalization a bit more complicated
     * 
     * @param halfWidth  size is 2*N + 1
     * @param nGauss     number of gaussians
     * @param sigGauss   Widths of the Gaussian Kernels
     * @param degGauss   Local spatial variation of bases
     */    
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> generateAlardLuptonKernelSet(
        unsigned int halfWidth,                ///< size is 2*N + 1
        unsigned int nGauss,                   ///< number of gaussians
        std::vector<double> const& sigGauss,   ///< width of the gaussians
        std::vector<int>    const& degGauss    ///< local spatial variation of gaussians
        );

    /*
     * Execute fundamental task of convolving template and subtracting it from science image
     */
    template <typename ImageT, typename BackgroundT>
    lsst::afw::image::MaskedImage<ImageT> convolveAndSubtract(
        lsst::afw::image::MaskedImage<ImageT> const& imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT> const& imageToNotConvolve,
        lsst::afw::math::Kernel const& convolutionKernel,
        BackgroundT background,
        bool invert=true
        );

    template <typename ImageT, typename BackgroundT>
    lsst::afw::image::MaskedImage<ImageT> convolveAndSubtract(
        lsst::afw::image::Image<ImageT> const& imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT> const& imageToNotConvolve,
        lsst::afw::math::Kernel const& convolutionKernel,
        BackgroundT background,
        bool invert=true
        );

    /** Search through images for Footprints with no masked pixels
     *
     * @note Uses Eigen math backend
     *
     * @param imageToConvolve  MaskedImage to convolve with Kernel
     * @param imageToNotConvolve  MaskedImage to subtract convolved template from
     * @param policy  Policy for operations; in particular object detection
     */    
    template <typename ImageT>
    std::vector<lsst::afw::detection::Footprint::Ptr> getCollectionOfFootprintsForPsfMatching(
        lsst::afw::image::MaskedImage<ImageT> const& imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT> const& imageToNotConvolve,
        lsst::pex::policy::Policy const& policy
        );

    template <typename ImageT>
    Eigen::MatrixXd imageToEigenMatrix(
        lsst::afw::image::Image<ImageT> const& img
        );
    
    /** Functor to create PSF Matching Kernel
     *
     * @ingroup diffim
     *
     * @note If constructed with a regularization matrix, will use it by default
     * 
     */
    template <typename ImageT, typename VarT=lsst::afw::image::VariancePixel>
    class PsfMatchingFunctor {
    public:
        typedef boost::shared_ptr<PsfMatchingFunctor> Ptr;
        typedef typename lsst::afw::image::MaskedImage<ImageT>::xy_locator xy_locator;
        typedef typename lsst::afw::image::Image<VarT>::xy_locator         xyi_locator;

        PsfMatchingFunctor(
            lsst::afw::math::KernelList<lsst::afw::math::Kernel> const& basisList
            );
        PsfMatchingFunctor(
            lsst::afw::math::KernelList<lsst::afw::math::Kernel> const& basisList,
            Eigen::MatrixXd const H
            );
        virtual ~PsfMatchingFunctor() {};

        /** Return background value
         */
        double getBackground()                   const { return _background; }

        /** Return uncertainty on background value
         */
        double getBackgroundError()              const { return _backgroundError; }

        /** Return PSF matching kernel
         */
        boost::shared_ptr<lsst::afw::math::Kernel> getKernel()      const { return _kernel; }

        /** Return uncertainty on matching kernel, as kernel itself
         */
        boost::shared_ptr<lsst::afw::math::Kernel> getKernelError() const { return _kernelError; }

        /* Create PSF matching kernel */
        void apply(lsst::afw::image::Image<ImageT> const& imageToConvolve,
                   lsst::afw::image::Image<ImageT> const& imageToNotConvolve,
                   lsst::afw::image::Image<VarT>   const& varianceEstimate,
                   lsst::pex::policy::Policy       const& policy
            );

        /* Slow and outdated version */
        void apply2(lsst::afw::image::Image<ImageT> const& imageToConvolve,
                    lsst::afw::image::Image<ImageT> const& imageToNotConvolve,
                    lsst::afw::image::Image<VarT>   const& varianceEstimate,
                    lsst::pex::policy::Policy       const& policy
            );

    protected:
        lsst::afw::math::KernelList<lsst::afw::math::Kernel> _basisList;        ///< List of Kernel basis functions
        Eigen::MatrixXd _H;                                                     ///< Regularization matrix
        double _background;                                                     ///< Differenaitl background estimate
        double _backgroundError;                                                ///< Uncertainty on background
        boost::shared_ptr<lsst::afw::math::Kernel> _kernel;                     ///< PSF matching kernel
        boost::shared_ptr<lsst::afw::math::Kernel> _kernelError;                ///< Uncertainty on kernel
    };
    
}}}

#endif



