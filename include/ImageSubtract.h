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
namespace imageproc {

    class Source : private lsst::fw::LsstBase
    {

    public:
        explicit Source();
        virtual ~Source() {};
        explicit Source(
            unsigned rowc,
            unsigned colc,
            unsigned drow,
            unsigned dcol
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
    };
    
    template <class PixelT, class MaskT, class KernelT>
    void computePSFMatchingKernelForMaskedImage(
        lsst::fw::MaskedImage<PixelT,MaskT> const &imageToConvolve,
        lsst::fw::MaskedImage<PixelT,MaskT> const &imageToNotConvolve,
        lsst::fw::LinearCombinationKernel<KernelT> &kernelBasisSet);
    
    template <class PixelT, class MaskT, class KernelT>
    void computePSFMatchingKernelForPostageStamp(
        lsst::fw::MaskedImage<PixelT, MaskT>  &imageToConvolve,
        lsst::fw::MaskedImage<PixelT, MaskT> const &imageToNotConvolve,
        lsst::fw::LinearCombinationKernel<KernelT> &kernelBasisSet,
        std::vector<KernelT> &kernelCoeffs);

    void getCollectionOfMaskedImagesForPSFMatching(
        vector<Source> &sourceCollection
        );

    void getTemplateChunkExposureFromTemplateExposure();
    void wcsMatchExposure();
    void computeSpatiallyVaryingPSFMatchingKernel();
    void fitKernelsUsingPrincipalComponentAnalysis();
    void fitArraysUsingPrincipalComponentAnalysis();

}
}

//#include <ImageSubtract.cc>

#endif // !defined(LSST_Imageproc_ImageSubtract_H)
