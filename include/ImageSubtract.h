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

#include <lsst/fw/LsstBase.h>
#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/Source.h>

namespace lsst {
namespace imageproc {

    using namespace std;
    
    template <typename PixelT, typename MaskT, typename KernelT>
    void computePSFMatchingKernelForMaskedImage(
        lsst::fw::MaskedImage<PixelT,MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<PixelT,MaskT> const &imageToNotConvolve,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelBasisVec);
    
    template <typename PixelT, typename MaskT, typename KernelT>
    void computePSFMatchingKernelForPostageStamp(
        lsst::fw::MaskedImage<PixelT, MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<PixelT, MaskT> const &imageToNotConvolve,
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelBasisVec,
        vector<double> &kernelCoeffs);

    void getCollectionOfMaskedImagesForPSFMatching(
        vector<lsst::fw::Source> &sourceCollection
        );

    template <typename KernelT>
    void computePCAKernelBasis(
        vector<lsst::fw::LinearCombinationKernel<KernelT> > const &kernelVec,
        vector<lsst::fw::Kernel<KernelT> > &kernelPCABasisVec
        );

    void getTemplateChunkExposureFromTemplateExposure();
    void wcsMatchExposure();
    void computeSpatiallyVaryingPSFMatchingKernel();

}
}

#include <ImageSubtract.cc>

#endif // !defined(LSST_Imageproc_ImageSubtract_H)
