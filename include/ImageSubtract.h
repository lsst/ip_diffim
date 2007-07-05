// -*- lsst-c++ -*-
#ifndef LSST_FW_ImageSubtract_H
#define LSST_FW_ImageSubtract_H

namespace lsst {
namespace fw {

    template <PixelT, MaskT, KernelT>
    static void computePSFMatchingKernelForPostageStamp(
        MaskedImage<PixelT, MaskT> const &convolveImage,
        MaskedImage<PixelT, MaskT> const &nonconvolveImage,
        LinearCombinationKernel<KernelT> &kernelSet,
        std::vector<KernelT> &kernelCoeffs);
    
}
}

#endif // !defined(LSST_FW_ImageSubtract_H)
