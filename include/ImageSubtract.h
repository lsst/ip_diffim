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
    
    template <typename PixelT, typename MaskT, typename KernelT>
    void computePSFMatchingKernelForMaskedImage(
        lsst::fw::MaskedImage<PixelT,MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<PixelT,MaskT> const &imageToNotConvolve,
        lsst::fw::LinearCombinationKernel<KernelT> &kernelBasisSet);
    
    template <typename PixelT, typename MaskT, typename KernelT>
    void computePSFMatchingKernelForPostageStamp(
        lsst::fw::MaskedImage<PixelT, MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<PixelT, MaskT> const &imageToNotConvolve,
        lsst::fw::LinearCombinationKernel<KernelT> &kernelBasisSet,
        std::vector<KernelT> &kernelCoeffs);

    void getCollectionOfMaskedImagesForPSFMatching(
        vector<lsst::fw::Source> &sourceCollection
        );

    void getTemplateChunkExposureFromTemplateExposure();
    void wcsMatchExposure();
    void computeSpatiallyVaryingPSFMatchingKernel();
    void fitKernelsUsingPrincipalComponentAnalysis();
    void fitArraysUsingPrincipalComponentAnalysis();

}
}

#include <ImageSubtract.cc>

#endif // !defined(LSST_Imageproc_ImageSubtract_H)
