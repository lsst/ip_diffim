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

#include <boost/shared_ptr.hpp>

#include <lsst/pex/policy/Policy.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/math/KernelFunctions.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/math/Function.h>
#include <lsst/detection/Footprint.h>

namespace lsst {
namespace ip {
namespace diffim {
    
    /** 
     * 
     * @brief Class to store the summary statistics of a difference MaskedImage
     * 
     * @ingroup ip_diffim
     *
     * @note This class will be deprecated as soon as Sdqa classes are implemented
     */
    template <typename ImageT, typename MaskT>
    class DifferenceImageStatistics {
    public:
        DifferenceImageStatistics();
        DifferenceImageStatistics(const lsst::afw::image::MaskedImage<ImageT, MaskT> differenceMaskedImage);

        virtual ~DifferenceImageStatistics() {};
        void setResidualMean(double mean) {_residualMean = mean;}
        void setResidualStd(double std) {_residualStd = std;}
        double getResidualMean() {return _residualMean;}
        double getResidualStd() {return _residualStd;}

        bool evaluateQuality(lsst::pex::policy::Policy &policy);
    private:
        double _residualMean;
        double _residualStd;
    };

    /** Uses a functor to iterate over the Mask image pixels
     *
     * @ingroup diffim
     *
     * @note Will replace maskOK subroutine.  Use this as an example to sort
     * through footprints for the "best" ones for difference imaging.
     *
     * @note Need up update detection first
     * 
     * Example usage : 
     *  FindSetBits<image::Mask<image::MaskPixel> > count(mask); 
     *  count.reset(); 
     *  count.apply(footprint); 
     *  nSet = count.getBits();
     * 
     */

    /*
    template <typename MaskT>
    class FindSetBits : public lsst::detection::FootprintFunctor<MaskT> {
    public:
        FindSetBits(MaskT const& mask) : 
            lsst::detection::FootprintFunctor<MaskT>(mask), _bits(0) {;}
        
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
    */


    /** Build a set of Delta Function basis kernels
     *
     * @param nCols  Number of rows in the set
     * @param nRows  Number of columns in the set
     */    
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> generateDeltaFunctionKernelSet(
        unsigned int nCols,
        unsigned int nRows
        );

    /** Build a set of Alard/Lupton basis kernels
     *
     * @note NOT IMPLEMENTED
     *
     * @param nCols  Number of rows in the set
     * @param nRows  Number of columns in the set
     * @param sigGauss  Widths of the Gaussian Kernels
     * @param degGauss  Local spatial variation of bases
     */    
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> generateAlardLuptonKernelSet(
        unsigned int nCols,
        unsigned int nRows,
        std::vector<double> const &sigGauss,
        std::vector<double> const &degGauss
        );

    /** Execute fundamental task of convolving template and subtracting it from science image
     *
     * @note D = I - (K.x.T + bg)
     * 
     * @param imageToConvolve  MaskedImage to convolve with Kernel
     * @param imageToNotConvolve  MaskedImage to subtract convolved template from
     * @param convolutionKernelPtr  PSF-matching Kernel used for convolution
     * @param background  Differential background value
     */    
    template <typename ImageT, typename MaskT>
    lsst::afw::image::MaskedImage<ImageT, MaskT> convolveAndSubtract(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        boost::shared_ptr<lsst::afw::math::Kernel> const &convolutionKernelPtr,
        double background
        );

    /** Execute fundamental task of convolving template and subtracting it from science image
     *
     * @note D = I - (K.x.T + bg)
     * 
     * @param imageToConvolve  MaskedImage to convolve with Kernel
     * @param imageToNotConvolve  MaskedImage to subtract convolved template from
     * @param convolutionKernelPtr  PSF-matching LinearCombinationKernelKernel used for convolution
     * @param background  Differential background value
     */    
    template <typename ImageT, typename MaskT>
    lsst::afw::image::MaskedImage<ImageT, MaskT> convolveAndSubtract(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> const &convolutionKernelPtr,
        double background
        );

    /** Execute fundamental task of convolving template and subtracting it from science image
     *
     * @note D = I - (K.x.T + bg)
     * 
     * @param imageToConvolve  MaskedImage to convolve with Kernel
     * @param imageToNotConvolve  MaskedImage to subtract convolved template from
     * @param convolutionKernelPtr  PSF-matching Kernel used for convolution
     * @param background  Differential background function
     */    
    template <typename ImageT, typename MaskT, typename FunctionT>
    lsst::afw::image::MaskedImage<ImageT, MaskT> convolveAndSubtract(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        boost::shared_ptr<lsst::afw::math::Kernel> const &convolutionKernelPtr,
        boost::shared_ptr<lsst::afw::math::Function2<FunctionT> > backgroundFunction
        );

    /** Execute fundamental task of convolving template and subtracting it from science image
     *
     * @note D = I - (K.x.T + bg)
     * 
     * @param imageToConvolve  MaskedImage to convolve with Kernel
     * @param imageToNotConvolve  MaskedImage to subtract convolved template from
     * @param convolutionKernelPtr  PSF-matching LinearCombinationKernel used for convolution
     * @param background  Differential background function
     */    
    template <typename ImageT, typename MaskT, typename FunctionT>
    lsst::afw::image::MaskedImage<ImageT, MaskT> convolveAndSubtract(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> const &convolutionKernelPtr,
        boost::shared_ptr<lsst::afw::math::Function2<FunctionT> > backgroundFunction
        );

    /** Search through images for Footprints with no masked pixels
     *
     * @param imageToConvolve  MaskedImage to convolve with Kernel
     * @param imageToNotConvolve  MaskedImage to subtract convolved template from
     * @param policy  Policy for operations; in particular object detection
     */    
    template <typename ImageT, typename MaskT>
    std::vector<lsst::detection::Footprint::PtrType> getCollectionOfFootprintsForPsfMatching(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::pex::policy::Policy &policy
        );
    
    /** Build a single PSF-matching Kernel for a Footprint; core of ip_diffim processing
     *
     * @note This version uses an input MaskedImage as an estimate of the variance
     *
     * @param imageToConvolve  MaskedImage to convolve with Kernel
     * @param imageToNotConvolve  MaskedImage to subtract convolved template from
     * @param varianceImage  Estimate of diffim variance
     * @param kernelInBasisList  Input kernel basis set
     * @param policy  Policy for operations; in particular object detection
     *
     * @param kernelPtr  Pointer to resulting PSF matching kernel
     * @param kernelErrorPtr  Uncertainty on PSF matching kernel
     * @param background  Differential background
     * @param backgroundError  Uncertainty on differential background
     */    
    template <typename ImageT, typename MaskT>
    void computePsfMatchingKernelForFootprint(
        lsst::afw::image::MaskedImage<ImageT, MaskT>         const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT>         const &imageToNotConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT>         const &varianceImage,
        lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelInBasisList,
        lsst::pex::policy::Policy                  &policy,
        boost::shared_ptr<lsst::afw::math::Kernel> &kernelPtr,
        boost::shared_ptr<lsst::afw::math::Kernel> &kernelErrorPtr,
        double                                     &background,
        double                                     &backgroundError
        );

    /** Build a single PSF-matching Kernel for a Footprint; core of ip_diffim processing
     *
     * @note This version uses an input MaskedImage as an estimate of the
     * variance, as well as GSL for the matrix math.
     *
     * @param imageToConvolve  MaskedImage to convolve with Kernel
     * @param imageToNotConvolve  MaskedImage to subtract convolved template from
     * @param varianceImage  Estimate of diffim variance
     * @param kernelInBasisList  Input kernel basis set
     * @param policy  Policy for operations; in particular object detection
     */    
    template <typename ImageT, typename MaskT>
    std::vector<std::pair<double,double> > computePsfMatchingKernelForFootprintGSL(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &varianceImage,
        lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelInBasisList,
        lsst::pex::policy::Policy &policy
        );

    /** Build a single PSF-matching Kernel for a Footprint; core of ip_diffim processing
     *
     * @param background  Differential background value
     * @param imageToConvolve  MaskedImage to convolve with Kernel
     * @param imageToNotConvolve  MaskedImage to subtract convolved template from
     * @param kernelInBasisList  Input kernel basis set
     * @param policy  Policy for operations; in particular object detection
     */    
    template <typename ImageT, typename MaskT>
    std::vector<double> computePsfMatchingKernelForFootprint_Legacy(
        double &background,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelInBasisList,
        lsst::pex::policy::Policy &policy
        );

    /** Searches through Mask for bit badPixelMask
     *
     * @param inputMask  Input Mask to search through
     * @param badPixelMask  Mask bit to search for
     */
    template <typename MaskT>
    bool maskOk(
        lsst::afw::image::Mask<MaskT> const &inputMask,
        MaskT const badPixelMask
        );

    /** Calculate pixel statistics of a MaskedImage
     *
     * @note Typically run on a difference image
     *
     * @note This should eventually be replaced by afw::math functions
     *
     * @param nGoodPixels  Returned number of pixels in the calculation
     * @param mean  Returned mean of pixel values
     * @param variance  Returned variance of pixel values
     * @param inputImage  MaskedImage to calculate statistics for
     * @param badPixelMask  Mask bit to ignore
     */
    template <typename ImageT, typename MaskT>
    void calculateMaskedImageStatistics(
        int &nGoodPixels,
        double &mean,
        double &variance,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &inputImage,
        MaskT const badPixelMask
        );

    /** Calculate pixel statistics of a MaskedImage
     *
     * @note Typically run on a difference image; this version ignores Mask values
     *
     * @note This should eventually be replaced by afw::math functions
     *
     * @param nGoodPixels  Returned number of pixels in the calculation
     * @param mean  Returned mean of pixel values
     * @param variance  Returned variance of pixel values
     * @param inputImage  MaskedImage to calculate statistics for
     */
    template <typename ImageT, typename MaskT>
    void calculateMaskedImageStatistics(
        int &nGoodPixels,
        double &mean,
        double &variance,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &inputImage
        );

    /** Calculate pixel statistics of an Image
     *
     * @note This should eventually be replaced by afw::math functions
     *
     * @param nGoodPixels  Returned number of pixels in the calculation
     * @param mean  Returned mean of pixel values
     * @param variance  Returned variance of pixel values
     * @param inputImage  Image to calculate statistics for
     */
    template <typename ImageT>
    void calculateImageStatistics(
        int &nGoodPixels,
        double &mean,
        double &variance,
        lsst::afw::image::Image<ImageT> const &inputImage
        );

    /** Calculate statistics of a Vector
     *
     * @note This should eventually be replaced by afw::math functions
     *
     * @param inputVector  Input data
     * @param mean  Returned mean of values
     * @param variance  Returned variance of values
     */
    template <typename VectorT>
    void calculateVectorStatistics(
        vw::Vector<VectorT> const &inputVector,
        double &mean,
        double &variance
        );

    /** Add a spatially varying function to an Image
     *
     * @note Typically used to add a background Function to an Image
     *
     * @param image Image to add function to
     * @param function  Function added to image
     */
    template <typename PixelT, typename FunctionT>
    void addFunctionToImage(
        lsst::afw::image::Image<PixelT> &image,
        lsst::afw::math::Function2<FunctionT> const &function
        );


}}}

#endif



