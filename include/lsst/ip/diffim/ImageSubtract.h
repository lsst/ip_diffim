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

namespace lsst {
namespace ip {
namespace diffim {
    
    /**
     * @brief Class to store the summary statistics of a difference MaskedImage
     * 
     * @ingroup diffim
     */
    template <typename ImageT, typename MaskT>
    class DifferenceImageStatistics : public lsst::daf::base::Persistable,
                                      public lsst::daf::data::LsstBase {
    public:
        DifferenceImageStatistics();
        DifferenceImageStatistics(const lsst::afw::image::MaskedImage<ImageT, MaskT> differenceMaskedImage);

        virtual ~DifferenceImageStatistics() {};
        void setResidualMean(double mean) {_residualMean = mean;}
        void setResidualVariance(double variance) {_residualVariance = variance;}
        double getResidualMean() {return _residualMean;}
        double getResidualVariance() {return _residualVariance;}
    private:
        double _residualMean;
        double _residualVariance;
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
    class DifferenceImageFootprintInformation : public lsst::daf::base::Persistable,
                                                public lsst::daf::data::LsstBase {
    public:
        typedef boost::shared_ptr<DifferenceImageFootprintInformation<ImageT, MaskT> > Ptr;
        typedef std::vector<typename DifferenceImageFootprintInformation<ImageT, MaskT>::Ptr> DifiList;

        // or boost::shared_ptr<lsst::afw::image::MaskedImage<ImageT, MaskT> >
        /* SWIG DOES NOT LIKE BELOW */
        //typedef lsst::afw::image::MaskedImage<ImageT, MaskT> MaskedImage;
        //typedef boost::shared_ptr<MaskedImage> MaskedImagePtr;
        typedef boost::shared_ptr<lsst::afw::image::MaskedImage<ImageT, MaskT> > MaskedImagePtr; 

        /* No empty constructor */
        //DifferenceImageFootprintInformation();

        DifferenceImageFootprintInformation(lsst::detection::Footprint::PtrType footprintPtr,
                                            MaskedImagePtr imageToConvolvePtr,
                                            MaskedImagePtr imageToNotConvolvePtr);
        virtual ~DifferenceImageFootprintInformation() {};

        void setID(int id) {_id = id;};
        int getID() {return _id;};
        void setColcNorm(double colc) {_colcNorm = colc;};
        void setRowcNorm(double rowc) {_rowcNorm = rowc;};
        double getColcNorm() {return _colcNorm;};
        double getRowcNorm() {return _rowcNorm;};

        /* put this functionality in the constructor */
        //void setFootprintPtr(lsst::detection::Footprint::PtrType ptr) {_footprintPtr = ptr;};
        //void setImageToNotConvolvePtr(MaskedImagePtr ptr) {_imageToNotConvolvePtr = ptr;};
        //void setImageToConvolvePtr(MaskedImagePtr ptr) {_imageToConvolvePtr = ptr;};

        lsst::detection::Footprint::PtrType getFootprintPtr() {return _footprintPtr;};
        MaskedImagePtr getImageToNotConvolvePtr() {return _imageToNotConvolvePtr;};
        MaskedImagePtr getImageToConvolvePtr() {return _imageToConvolvePtr;};

        void setSingleKernelPtr(boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> ptr) {_singleKernelPtr = ptr;};
        boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> getSingleKernelPtr() {return _singleKernelPtr;};

        void setSingleKernelSum(double kernelSum) {_singleKernelSum = kernelSum;};
        double getSingleKernelSum() {return _singleKernelSum;};

        void setSingleBackground(double background) {_singleBackground = background;};
        double getSingleBackground() {return _singleBackground;};

        void setSingleStats(DifferenceImageStatistics<ImageT, MaskT> stats) {_singleKernelStats = stats;};
        DifferenceImageStatistics<ImageT, MaskT> getSingleStats() {return _singleKernelStats;};

        DifferenceImageStatistics<ImageT, MaskT> computeImageStatistics(boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> kernelPtr,
                                                                        double background);

        void setStatus(bool status) {_isGood = status;};
        bool getStatus() {return _isGood;};

    private:
        /* running ID */
        int _id;

        /* position of the Footprint in the image, -1 to 1 */
        double _colcNorm;
        double _rowcNorm;

        /* footprint assocated with the object we're building the kernels around */
        lsst::detection::Footprint::PtrType _footprintPtr;

        /* subimages associated with the Footprint */
        MaskedImagePtr _imageToConvolvePtr;    /* Typically the template image */
        MaskedImagePtr _imageToNotConvolvePtr; /* Typically the science image */

        /* results from individual kernel fit */
        boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> _singleKernelPtr;
        double _singleKernelSum;
        double _singleBackground;
        DifferenceImageStatistics<ImageT, MaskT> _singleKernelStats;

        bool _isGood;
    };
    
    template <typename ImageT, typename MaskT>
    typename DifferenceImageFootprintInformation<ImageT, MaskT>::DifiList
    getGoodFootprints(std::vector<boost::shared_ptr<DifferenceImageFootprintInformation<ImageT, MaskT> > > & difiList );
                                                                     
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
        boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> const &convolutionKernelPtr,
        double background
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
