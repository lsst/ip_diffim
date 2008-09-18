// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Implementation of Image Subtraction
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup diffim
 */

#ifndef LSST_IMAGEPROC_IMAGESUBTRACT_H
#define LSST_IMAGEPROC_IMAGESUBTRACT_H

#include <vector>

#include <boost/shared_ptr.hpp>

#include <lsst/pex/policy/Policy.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/math/KernelFunctions.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/math/Function.h>
#include <lsst/detection/Footprint.h>
#include <lsst/ip/diffim/Pca.h>

namespace lsst {
namespace ip {
namespace diffim {
    
    /**
     * @brief Class to store the summary statistics of a difference MaskedImage
     * 
     * @ingroup diffim
     */
    class DifferenceImageStatistics {
    public:
        explicit DifferenceImageStatistics() {};
        virtual ~DifferenceImageStatistics() {};
        void setFootprintResidualMean(double mean) {footprintResidualMean = mean;}
        void setFootprintResidualVariance(double variance) {footprintResidualVariance = variance;}
        double getFootprintResidualMean() {return footprintResidualMean;}
        double getFootprintResidualVariance() {return footprintResidualVariance;}
    private:
        double footprintResidualMean;
        double footprintResidualVariance;
    };

    /**
     * @brief Class containing Image Subtraction information for a given
     * Footprint
     * 
     * For each set of input template,science MaskedImages, the Detection
     * pipeline or a database query will return a list of acceptible coordinates
     * around which to build a difference imaging Kernel.  Each of these
     * positions will be assigned a DifferenceImageFootprintInformation instance
     * containing pointers to the Footprint generated at this position, the
     * row,col position itself, and associated subimages in the MaskedImages.
     * 
     * @ingroup diffim
     */
    template <typename ImageT, typename MaskT>
    class DifferenceImageFootprintInformation {
    public:
        typedef lsst::afw::image::MaskedImage<ImageT, MaskT> maskedImageType;
        typedef boost::shared_ptr<maskedImageType> maskedImagePtrType;

        explicit DifferenceImageFootprintInformation() {};
        virtual ~DifferenceImageFootprintInformation() {};

        void setColcNorm(double colc) {colcNorm = colc;};
        void setRowcNorm(double rowc) {rowcNorm = rowc;};
        double getColcNorm() {return colcNorm;};
        double getRowcNorm() {return rowcNorm;};

        void setFootprintPtr(lsst::detection::Footprint::PtrType ptr) {footprintPtr = ptr;};
        lsst::detection::Footprint::PtrType getFootprintPtr() {return footprintPtr;};

        void setImageToNotConvolvePtr(maskedImagePtrType ptr) {imageToNotConvolvePtr = ptr;};
        maskedImagePtrType getImageToNotConvolvePtr() {return imageToNotConvolvePtr;};

        void setImageToConvolvePtr(maskedImagePtrType ptr) {imageToConvolvePtr = ptr;};
        maskedImagePtrType getImageToConvolvePtr() {return imageToConvolvePtr;};

        void setSingleKernelPtr(boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> ptr) {singleKernelPtr = ptr;};
        boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> getSingleKernelPtr() {return singleKernelPtr;};

        void setSingleKernelSum(double kernelSum) {singleKernelSum = kernelSum;};
        double getSingleKernelSum() {return singleKernelSum;};

        void setSingleBackground(double background) {singleBackground = background;};
        double getSingleBackground() {return singleBackground;};

        void setSingleStats(DifferenceImageStatistics stats) {singleKernelStats = stats;};
        DifferenceImageStatistics getSingleStats() {return singleKernelStats;};

        DifferenceImageStatistics computeImageStatistics(boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> kernel,
                                                         double background);

    private:
        /* position of the Footprint in the image, -1 to 1 */
        double colcNorm;
        double rowcNorm;

        /* footprint assocated with the object we're building the kernels around */
        lsst::detection::Footprint::PtrType footprintPtr;

        /* subimages associated with the Footprint */
        maskedImagePtrType imageToNotConvolvePtr; /* Typically the science image */
        maskedImagePtrType imageToConvolvePtr;    /* Typically the template image */

        /* results from individual kernel fit */
        boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> singleKernelPtr;
        double singleKernelSum;
        double singleBackground;
        DifferenceImageStatistics singleKernelStats;
    };




    lsst::afw::math::KernelList<lsst::afw::math::Kernel> generateDeltaFunctionKernelSet(
        unsigned int nCols,
        unsigned int nRows
        );

    lsst::afw::math::KernelList<lsst::afw::math::Kernel> generateAlardLuptonKernelSet(
        unsigned int nCols,
        unsigned int nRows,
        std::vector<double> const &sigGauss,
        std::vector<double> const &degGauss
        );

    template <typename ImageT, typename MaskT>
    lsst::afw::image::MaskedImage<ImageT, MaskT> convolveAndSubtract(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> convolutionKernel,
        double background,
        lsst::pex::policy::Policy &policy
        );

    template <typename ImageT, typename MaskT>
    std::vector<lsst::detection::Footprint::PtrType> getCollectionOfFootprintsForPsfMatching(
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::pex::policy::Policy &policy
        );

    template <typename ImageT, typename MaskT>
    std::vector<double> computePsfMatchingKernelForFootprint(
        double &background,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelInBasisList, ///< Input kernel basis set
        lsst::pex::policy::Policy &policy
        );

    template <typename MaskT>
    bool maskOk(
        lsst::afw::image::Mask<MaskT> const &inputMask,
        MaskT const badPixelMask
        );

    template <typename ImageT, typename MaskT>
    void calculateMaskedImageStatistics(
        int &nGoodPixels,
        double &mean,
        double &variance,
        lsst::afw::image::MaskedImage<ImageT, MaskT> const &inputImage,
        MaskT const badPixelMask
        );

    template <typename ImageT>
    void calculateImageStatistics(
        int &nGoodPixels,
        double &mean,
        double &variance,
        lsst::afw::image::Image<ImageT> const &inputImage
        );

    template <typename VectorT>
    void calculateVectorStatistics(
        vw::Vector<VectorT> const &inputVector,
        double &mean,
        double &variance
        );

    template <typename PixelT, typename FunctionT>
    void addFunctionToImage(
        lsst::afw::image::Image<PixelT> &image,
        lsst::afw::math::Function2<FunctionT> const &function
        );

}}}

#endif
