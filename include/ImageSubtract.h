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

namespace lsst {
namespace fw {

    class Object : private LsstBase{
    public:
        explicit Object();
        virtual ~Object() {};
        explicit Object(
            unsigned rowc,
            unsigned colc,
            unsigned drow,
            unsigned dcol,
        );
       inline unsigned getColc() const;
       inline unsigned getRowc() const;
       inline unsigned getDcol() const;
       inline unsigned getDrow() const;
    private:
       unsigned _rowc;
       unsigned _colc;
       unsigned _drow;
       unsigned _dcol;
    }
    
    template <PixelT, MaskT, KernelT>
    static void computePSFMatchingKernelForMaskedImage(
        MaskedImage<PixelT,MaskT> const &imageToConvolve,
        MaskedImage<PixelT,MaskT> const &imageToNotConvolve,
        LinearCombinationKernel<KernelT> &kernelBasisSet);
    
    template <PixelT, MaskT, KernelT>
    static void computePSFMatchingKernelForPostageStamp(
        MaskedImage<PixelT, MaskT> const &imageToConvolve,
        MaskedImage<PixelT, MaskT> const &imageToNotConvolve,
        LinearCombinationKernel<KernelT> &kernelBasisSet,
        std::vector<KernelT> &kernelCoeffs);

    static void getCollectionOfMaskedImagesForPSFMatching(
        vector<Object> &objectCollection);

    static void getTemplateChunkExposureFromTemplateExposure(
        );
    static void wcsMatchExposure(
        );
    static void computeSpatiallyVaryingPSFMatchingKernel(
        );
    static void fitKernelsUsingPrincipalComponentAnalysis(
        );
    static void fitArraysUsingPrincipalComponentAnalysis(
        );

}
}

#include <ImageSubtract.cc>

#endif // !defined(LSST_Imageproc_ImageSubtract_H)
