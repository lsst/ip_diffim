// -*- lsst-c++ -*-
#ifndef LSST_Imageproc_ImageSubtract_H
#define LSST_Imageproc_ImageSubtract_H
/**
 * \file
 *
 * Implementation of Image Subtraction
 *
 * \author Andrew Becker
 *
 * \ingroup imageproc
 */
#include <vector>

#include <boost/shared_ptr.hpp>

#include <lsst/mwi/policy/Policy.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/Function.h>
#include <lsst/detection/Footprint.h>

namespace lsst {
namespace imageproc {
    
    // This should be replaced later with a more general Info class
    // Perhaps a DataProperty collection?
    template <typename KernelT>
    struct DiffImContainer {
        // Running ID
        int id;

        // The footprint assocated with the structure
        lsst::detection::Footprint::PtrType diffImFootprintPtr;
        double colcNorm; // -1 to 1
        double rowcNorm; // -1 to 1
        
        // From single kernel fit
        double background;
        double kSum;
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > diffImKernelPtr;
        
        // If using PCA, this holds the PCA Kernel
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > diffImPcaKernelPtr;
        
        /* Goodness of fit metrics */

        /* From initial kernel fit */
        double footprintResidualMean;
        double footprintResidualVariance;

        /* From PCA kernel fit */
        double kernelResidual;
        double kernelResidualVariance;

        /* From spatial kernel fit */
        double spatialKernelResidual;
        double spatialKernelResidualVariance;
        double footprintResidualMean2;
        double footprintResidualVariance2;
        
        // Flags
        bool isGood;
        
        DiffImContainer()
            {
                isGood = true;
            }
    };

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
    
    template <typename ImageT, typename MaskT, typename KernelT>
    std::vector<double> computePsfMatchingKernelForPostageStamp(
        double &background,
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelInBasisList,
        lsst::mwi::policy::Policy &policy
        );
    
    template <typename ImageT, typename MaskT>
    std::vector<lsst::detection::Footprint::PtrType> getCollectionOfFootprintsForPsfMatching(
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        lsst::mwi::policy::Policy &policy
        );

    std::vector<lsst::detection::Footprint::PtrType> getCollectionOfMaskedImagesForPsfMatching(
        );

    template <typename KernelT>
    std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > computePcaKernelBasis(
        std::vector<lsst::imageproc::DiffImContainer<KernelT> > &diffImContainerList,
        lsst::mwi::policy::Policy &policy
        );

    template <typename KernelT>
    boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > computeSpatiallyVaryingPsfMatchingKernel(
        boost::shared_ptr<lsst::fw::function::Function2<double> > &kernelFunctionPtr,
        boost::shared_ptr<lsst::fw::function::Function2<double> > &backgroundFunctionPtr,
        std::vector<lsst::imageproc::DiffImContainer<KernelT> > &diffImContainerList,
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

}} // lsst::imageproc

#ifndef SWIG // don't bother SWIG with .cc files
#include <lsst/imageproc/ImageSubtract.cc>
#endif

#endif // !defined(LSST_Imageproc_ImageSubtract_H)
