// -*- lsst-c++ -*-
/**
 * @file ImageSubtract.h
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
#include <lsst/afw/Kernel.h>
#include <lsst/afw/KernelFunctions.h>
#include <lsst/afw/Mask.h>
#include <lsst/afw/MaskedImage.h>
#include <lsst/afw/Function.h>
#include <lsst/detection/Footprint.h>
#include <lsst/ip/diffim/PCA.h>

namespace lsst {
namespace ip {
namespace diffim {
    
    /** @struct
     *
     * @brief Structure to store the summary statistics of a difference MaskedImage
     * 
     * @ingroup diffim
     */
    struct MaskedImageDiffimStats {
        double kernelResidual;
        double kernelResidualVariance;
        double footprintResidualMean;
        double footprintResidualVariance;
    };

    /** @class
     *
     * @brief Class containing Image Subtraction information for a given
     * Footprint
     * 
     * For each set of input template,science MaskedImages, the Detection
     * pipeline or a database query will return a list of acceptible coordinates
     * around which to build a difference imaging Kernel.  Each of these
     * positions will be assigned a DiffImContainer containing pointers to the
     * Footprint generated at this position, the row,col position itself, and
     * associated subimages in the MaskedImages.
     * 
     * During each iteration of Kernel building, we add the Kernel itself to a
     * vector, the associated difference image statistics, background
     * calculation, and kernel sum.
     * 
     * The class also contains a boolean member isGood that is set to False if
     * it fails statistical tests during the processing.
     * 
     * @ingroup diffim
     */
    template <typename ImageT, typename MaskT>
    class DiffImContainer {
    public:
        /* Running ID */
        int id;
        
        typedef lsst::afw::MaskedImage<ImageT, MaskT> maskedImageType;
        typedef boost::shared_ptr<maskedImageType> maskedImagePtrType;

        /* The Footprint assocated with the object we're building the kernels around */
        lsst::detection::Footprint::PtrType diffImFootprintPtr;
        /* Position in full image */
        double colcNorm;
        double rowcNorm;

        /* Subimages associated with the Footprint */
        maskedImagePtrType imageToNotConvolve; /* Typically the template image */
        maskedImagePtrType imageToConvolve;    /* Typically the science image */

        /* Store the kernel at each stage of build */
        int nKernel;
        /* Vector of kernels : single image kernel, PCA model, spatial model */
        std::vector<boost::shared_ptr<lsst::afw::LinearCombinationKernel> > kernelList;
        std::vector<lsst::ip::diffim::MaskedImageDiffimStats> diffimStats;
        std::vector<double> backgrounds;
        std::vector<double> kernelSums;

        /* Flags */
        bool isGood;
        
        explicit DiffImContainer()
            {
                isGood = true;
                nKernel = -1;
            }
        virtual ~DiffImContainer() {};

        /* These functions hack around a SWIG issue -- the swigged code thinks the mask type of */
        /* DiffImContainerD.imageToNotConvolve's is boost::uint16_t, not lsst::afw::maskType */
        /* (even though they are identical and it knows about lsst::afw::maskType everywhere else) */
        void setImageToNotConvolve(maskedImageType const &maskedImage) {
            this->imageToNotConvolve.reset(new maskedImageType(maskedImage));
        }
        void setImageToConvolve(maskedImageType const &maskedImage) {
            this->imageToConvolve.reset(new maskedImageType(maskedImage));
        }
        void addKernel(boost::shared_ptr<lsst::afw::LinearCombinationKernel> newKernel,
                       lsst::ip::diffim::MaskedImageDiffimStats newStats,
                       double background,
                       double kernelSum) {
            this->kernelList.push_back(newKernel);
            this->diffimStats.push_back(newStats);
            this->backgrounds.push_back(background);
            this->kernelSums.push_back(kernelSum);
            this->nKernel += 1;
        }
        double getFootprintResidualMean() { return this->diffimStats[this->nKernel].footprintResidualMean; }
        double getFootprintResidualVariance() { return this->diffimStats[this->nKernel].footprintResidualVariance; }
    };

    template <typename ImageT, typename MaskT>
    void computeDiffImStats(
        lsst::ip::diffim::DiffImContainer<ImageT, MaskT> &diffImContainer,
        int const kernelID,
        lsst::pex::policy::Policy &policy
        );

    template <typename ImageT, typename MaskT>
    std::vector<lsst::detection::Footprint::PtrType> getCollectionOfFootprintsForPsfMatching(
        lsst::afw::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::pex::policy::Policy &policy
        );

    template <typename ImageT, typename MaskT>
    std::vector<double> computePsfMatchingKernelForPostageStamp(
        double &background,
        lsst::afw::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::afw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::afw::KernelList<lsst::afw::Kernel> const &kernelInBasisList, ///< Input kernel basis set
        lsst::pex::policy::Policy &policy
        );
    
    template <typename ImageT, typename MaskT>
    lsst::afw::KernelList<lsst::afw::Kernel> computePcaKernelBasis(
        std::vector<lsst::ip::diffim::DiffImContainer<ImageT, MaskT> > &diffImContainerList,
        lsst::pex::policy::Policy &policy
        );

    template <typename ImageT, typename MaskT>
    boost::shared_ptr<lsst::afw::LinearCombinationKernel> computeSpatiallyVaryingPsfMatchingKernel(
        boost::shared_ptr<lsst::afw::function::Function2<double> > &kernelFunctionPtr,
        boost::shared_ptr<lsst::afw::function::Function2<double> > &backgroundFunctionPtr,
        std::vector<lsst::ip::diffim::DiffImContainer<ImageT, MaskT> > &diffImContainerList,
        lsst::afw::KernelList<lsst::afw::Kernel> const &kernelBasisList,
        lsst::pex::policy::Policy &policy
        );

    lsst::afw::KernelList<lsst::afw::Kernel> generateDeltaFunctionKernelSet(
        unsigned int nCols,
        unsigned int nRows
        );

    lsst::afw::KernelList<lsst::afw::Kernel> generateAlardLuptonKernelSet(
        unsigned int nCols,
        unsigned int nRows,
        std::vector<double> const &sigGauss,
        std::vector<double> const &degGauss
        );

    template <typename MaskT>
    bool maskOk(
        lsst::afw::Mask<MaskT> const &inputMask,
        MaskT const badPixelMask
        );

    template <typename ImageT, typename MaskT>
    void calculateMaskedImageResiduals(
        int &nGoodPixels,
        double &meanOfResiduals,
        double &varianceOfResiduals,
        lsst::afw::MaskedImage<ImageT, MaskT> const &inputImage,
        MaskT const badPixelMask
        );

    template <typename ImageT>
    void calculateImageResiduals(
        int &nGoodPixels,
        double &meanOfResiduals,
        double &varianceOfResiduals,
        lsst::afw::Image<ImageT> const &inputImage
        );

    template <typename VectorT>
    void calculateVectorStatistics(
        vw::Vector<VectorT> const &inputVector,
        double &mean,
        double &variance
        );

    template <typename PixelT, typename FunctionT>
    void addFunction(
        lsst::afw::Image<PixelT> &image,
        lsst::afw::function::Function2<FunctionT> const &function
        );

}}}

#endif
