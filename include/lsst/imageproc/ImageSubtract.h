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

#include <lsst/fw/Kernel.h>
#include <lsst/fw/MaskedImage.h>

#include <lsst/mwi/policy/Policy.h>

#include <lsst/detection/Footprint.h>

namespace lsst {
namespace imageproc {

    using namespace std;
    
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
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > diffImKernelPtr;
        
        // If using PCA, this holds the PCA Kernel
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > diffImPcaKernelPtr;
        
        // Goodness of fit metrics
        double footprintResidualMean;
        double footprintResidualVariance;
        double kernelResidual;
        double kernelResidualVariance;
        double spatialKernelResidual;
        double spatialKernelResidualVariance;
        
        // Flags
        bool isGood;
        
        DiffImContainer()
            {
                isGood = true;
            }
    };

    template <typename ImageT, typename MaskT, typename KernelT, typename FuncT>
    void computePsfMatchingKernelForMaskedImage(
        lsst::fw::MaskedImage<ImageT,MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<ImageT,MaskT> const &imageToNotConvolve,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelInBasisList,
        vector<lsst::detection::Footprint::PtrType> const &footprintList,
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > &kernelPtr,
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > &kernelFunctionPtr,
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > &backgroundFunctionPtr,
        lsst::mwi::policy::Policy &policy
        );
    
    template <typename ImageT, typename MaskT, typename KernelT>
    void computePsfMatchingKernelForPostageStamp(
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelInBasisList,
        vector<double> &kernelCoeffs,
        double &background,
        lsst::mwi::policy::Policy &policy
        );
    
    template <typename ImageT, typename MaskT>
    void getCollectionOfFootprintsForPsfMatching(
        lsst::fw::MaskedImage<ImageT,MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<ImageT,MaskT> const &imageToNotConvolve,
        vector<lsst::detection::Footprint::PtrType> footprintListOut,
        lsst::mwi::policy::Policy &policy
        );

    void getCollectionOfMaskedImagesForPsfMatching(
        vector<lsst::detection::Footprint::PtrType> &footprintList
        );

    template <typename KernelT>
    void computePcaKernelBasis(
        vector<lsst::imageproc::DiffImContainer<KernelT> > &diffImContainerList,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelPcaBasisList,
        lsst::mwi::policy::Policy &policy
        );

    template <typename KernelT, typename FuncT>
    void computeSpatiallyVaryingPsfMatchingKernel(
        vector<lsst::imageproc::DiffImContainer<KernelT> > &diffImContainerList,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelOutBasisList,
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > &spatiallyVaryingKernelPtr,
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > &kernelFunctionPtr,
        lsst::mwi::policy::Policy &policy
        );







    template <typename KernelT>
    void generateDeltaFunctionKernelSet(
        unsigned int const nRows,
        unsigned int const nCols,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelBasisList
        );

    template <typename KernelT>
    void generateAlardLuptonKernelSet(
        unsigned int const nRows,
        unsigned int const nCols,
        vector<double> const sigGauss,
        vector<double> const degGauss,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelBasisList
        );

    template <typename ImageT, typename MaskT>
    bool checkMaskedImageForDiffim(
        lsst::fw::MaskedImage<ImageT,MaskT> const &inputImage
        );

    template <typename ImageT, typename MaskT>
    void calculateMaskedImageResiduals(
        lsst::fw::MaskedImage<ImageT,MaskT> const &inputImage,
        int &nGoodPixels,
        double &meanOfResiduals,
        double &varianceOfResiduals
        );

    template <typename ImageT>
    void calculateImageResiduals(
        lsst::fw::Image<ImageT> const &inputImage,
        int &nGoodPixels,
        double &meanOfResiduals,
        double &varianceOfResiduals
        );
    

}
}

#include <lsst/imageproc/ImageSubtract.cc>

#endif // !defined(LSST_Imageproc_ImageSubtract_H)
