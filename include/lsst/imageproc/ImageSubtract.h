// -*- lsst-c++ -*-
/**
 * @file ImageSubtract.h
 *
 * @brief Implementation of Image Subtraction
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup imageproc
 */

#ifndef LSST_IMAGEPROC_IMAGESUBTRACT_H
#define LSST_IMAGEPROC_IMAGESUBTRACT_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include <lsst/mwi/policy/Policy.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/KernelFunctions.h>
#include <lsst/fw/Mask.h>
#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/Function.h>
#include <lsst/detection/Footprint.h>
#include <lsst/imageproc/PCA.h>

namespace lsst {
namespace imageproc {
    
    /** @struct
     *
     * @brief Structure to store the summary statistics of a difference MaskedImage
     * 
     * @ingroup imageproc
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
     * @ingroup imageproc
     */
    template <typename ImageT, typename MaskT, typename KernelT>
    class DiffImContainer {
    public:
        /* Running ID */
        int id;
        
        typedef lsst::fw::MaskedImage<ImageT, MaskT> maskedImageType;
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
        std::vector<boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > > kernelList;
        std::vector<lsst::imageproc::MaskedImageDiffimStats> diffimStats;
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
        /* DiffImContainerD.imageToNotConvolve's is boost::uint16_t, not lsst::fw::maskType */
        /* (even though they are identical and it knows about lsst::fw::maskType everywhere else) */
        void setImageToNotConvolve(maskedImageType const &maskedImage) {
            this->imageToNotConvolve.reset(new maskedImageType(maskedImage));
        }
        void setImageToConvolve(maskedImageType const &maskedImage) {
            this->imageToConvolve.reset(new maskedImageType(maskedImage));
        }
        void addKernel(boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > newKernel,
                       lsst::imageproc::MaskedImageDiffimStats newStats,
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

    template <typename ImageT, typename MaskT, typename KernelT>
    void computeDiffImStats(
        lsst::imageproc::DiffImContainer<ImageT, MaskT, KernelT> &diffImContainer,
        int const kernelID,
        lsst::mwi::policy::Policy &policy
        );
       

    /* DEPRECATED to python
    template <typename ImageT, typename MaskT, typename KernelT>
    boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > computePsfMatchingKernelForMaskedImage(
        boost::shared_ptr<lsst::fw::function::Function2<double> > &kernelFunctionPtr,
        boost::shared_ptr<lsst::fw::function::Function2<double> > &backgroundFunctionPtr,
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelInBasisList,
        std::vector<lsst::detection::Footprint::PtrType> const &footprintList,
        lsst::mwi::policy::Policy &policy
        );
    */

    template <typename ImageT, typename MaskT>
    std::vector<lsst::detection::Footprint::PtrType> getCollectionOfFootprintsForPsfMatching(
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::mwi::policy::Policy &policy
        );

    template <typename ImageT, typename MaskT, typename KernelT>
    std::vector<double> computePsfMatchingKernelForPostageStamp(
        double &background,
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelInBasisList,
        lsst::mwi::policy::Policy &policy
        );
    
    template <typename ImageT, typename MaskT, typename KernelT>
    std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > computePcaKernelBasis(
        std::vector<lsst::imageproc::DiffImContainer<ImageT, MaskT, KernelT> > &diffImContainerList,
        lsst::mwi::policy::Policy &policy
        );

    template <typename ImageT, typename MaskT, typename KernelT>
    boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > computeSpatiallyVaryingPsfMatchingKernel(
        boost::shared_ptr<lsst::fw::function::Function2<double> > &kernelFunctionPtr,
        boost::shared_ptr<lsst::fw::function::Function2<double> > &backgroundFunctionPtr,
        std::vector<lsst::imageproc::DiffImContainer<ImageT, MaskT, KernelT> > &diffImContainerList,
        std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelBasisList,
        lsst::mwi::policy::Policy &policy
        );

    template <typename KernelT>
    std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > generateDeltaFunctionKernelSet(
        unsigned int nCols,
        unsigned int nRows
        );

    template <typename KernelT>
    std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > generateAlardLuptonKernelSet(
        unsigned int nCols,
        unsigned int nRows,
        std::vector<double> const &sigGauss,
        std::vector<double> const &degGauss
        );

    template <typename MaskT>
    bool maskOk(
        lsst::fw::Mask<MaskT> const &inputMask,
        MaskT const badPixelMask
        );

    template <typename ImageT, typename MaskT>
    void calculateMaskedImageResiduals(
        int &nGoodPixels,
        double &meanOfResiduals,
        double &varianceOfResiduals,
        lsst::fw::MaskedImage<ImageT, MaskT> const &inputImage,
        MaskT const badPixelMask
        );

    template <typename ImageT>
    void calculateImageResiduals(
        int &nGoodPixels,
        double &meanOfResiduals,
        double &varianceOfResiduals,
        lsst::fw::Image<ImageT> const &inputImage
        );

    template <typename VectorT>
    void calculateVectorStatistics(
        vw::Vector<VectorT> const &inputVector,
        double &mean,
        double &variance
        );

    template <typename PixelT, typename FunctionT>
    void addFunction(
        lsst::fw::Image<PixelT> &image,
        lsst::fw::function::Function2<FunctionT> const &function
        );

}}

#endif
