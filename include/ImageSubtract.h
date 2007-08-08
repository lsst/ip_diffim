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
    
    template <typename ImageT, typename MaskT, typename KernelT>
    void computePSFMatchingKernelForMaskedImage(
        lsst::fw::MaskedImage<ImageT,MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<ImageT,MaskT> const &imageToNotConvolve,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelBasisVec);
    
    template <typename ImageT, typename MaskT, typename KernelT>
    void computePSFMatchingKernelForPostageStamp(
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelBasisVec,
        vector<double> &kernelCoeffs);

    void getCollectionOfMaskedImagesForPSFMatching(
        vector<lsst::fw::Source> &sourceCollection
        );

    template <typename KernelT>
    void computePCAKernelBasis(
        vector<lsst::fw::LinearCombinationKernel<KernelT> > const &kernelVec,
        vector<double> &kernelResidualsVec,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelPCABasisVec,
        lsst::fw::Image<KernelT> &meanImage,
        vw::math::Matrix<double> &kernelCoefficients
        );

    template <typename KernelT, typename ReturnT>
    void computeSpatiallyVaryingPSFMatchingKernel(
        lsst::fw::Image<KernelT> const &meanImage,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelPCABasisVec,
        vw::math::Matrix<double> const &kernelCoefficients,
        vector<lsst::fw::Source> const &sourceCollection,
        boost::shared_ptr<lsst::fw::function::Function2<ReturnT> > spatiallyVaryingFunctionPtr,
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > spatiallyVaryingKernelPtr
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
