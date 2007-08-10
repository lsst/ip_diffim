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
        vector<lsst::fw::LinearCombinationKernel<KernelT> > const &kernelVec,
        vector<double> &kernelResidualsVec,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelPCABasisVec,
        vw::math::Matrix<double> &kernelCoefficients
        );

    template <typename KernelT, typename FuncT>
    void computeSpatiallyVaryingPSFMatchingKernel(
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelOutBasisVec,
        vw::math::Matrix<double> const &kernelCoefficients,
        vector<lsst::fw::Source> const &sourceCollection,
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

}
}

#include <ImageSubtract.cc>

#endif // !defined(LSST_Imageproc_ImageSubtract_H)
