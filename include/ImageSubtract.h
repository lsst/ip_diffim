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

#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/Source.h>

namespace lsst {
namespace imageproc {

    using namespace std;

    template <typename KernelT>
    struct DiffImContainer {
        // Running ID
        int id;

        // The source assocated with the structure
        lsst::fw::Source diffImSource;
        
        // From single kernel fit
        double background;
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > diffImKernelPtr;
        
        // If using PCA, this holds the PCA Kernel
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > diffImPCAKernelPtr;
        
        // Goodness of fit metrics
        double sourceResidualMean;
        double sourceResidualVariance;
        double kernelResidual;
        double spatialKernelResidual;
        
        // Flags
        bool isGood;
        
        DiffImContainer()
            {
                isGood = true;
            }
    };
    
    template <typename ImageT, typename MaskT, typename KernelT, typename FuncT>
    void computePSFMatchingKernelForMaskedImage(
        lsst::fw::MaskedImage<ImageT,MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<ImageT,MaskT> const &imageToNotConvolve,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelInBasisVec,
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > &kernelPtr,
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > &kernelFunctionPtr,
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > &backgroundFunctionPtr
        );
    
    template <typename ImageT, typename MaskT, typename KernelT>
    void computePSFMatchingKernelForPostageStamp(
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelInBasisVec,
        vector<double> &kernelCoeffs,
        double &background
        );

    void getCollectionOfMaskedImagesForPSFMatching(
        vector<lsst::fw::Source> &sourceCollection
        );

    template <typename KernelT>
    void computePCAKernelBasis(
        vector<lsst::imageproc::DiffImContainer<KernelT> > &diffImContainerVec,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelPCABasisVec
        );

    template <typename KernelT, typename FuncT>
    void computeSpatiallyVaryingPSFMatchingKernel(
        vector<lsst::imageproc::DiffImContainer<KernelT> > &diffImContainerVec,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelOutBasisVec,
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > &spatiallyVaryingKernelPtr,
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > &kernelFunctionPtr
        );

    template <typename KernelT>
    void generateDeltaFunctionKernelSet(
        unsigned int const nRows,
        unsigned int const nCols,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelBasisVec
        );

    template <typename KernelT>
    void generateAlardLuptonKernelSet(
        unsigned int const nRows,
        unsigned int const nCols,
        vector<double> const sigGauss,
        vector<double> const degGauss,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelBasisVec
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

#include <ImageSubtract.cc>

#endif // !defined(LSST_Imageproc_ImageSubtract_H)
