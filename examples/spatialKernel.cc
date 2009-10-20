#include <cmath> 
#include <lsst/afw/math.h>
#include <lsst/afw/image.h>
#include <lsst/ip/diffim.h>

namespace afwImage = lsst::afw::image;
namespace afwMath  = lsst::afw::math;
using namespace lsst::ip::diffim;

typedef float PixelT;

int main() {
    int kSize = 5;
    unsigned int nBases = kSize * kSize;
    int spatialKernelOrder = 2;
    
    lsst::afw::math::KernelList basisList = generateDeltaFunctionBasisSet(kSize, kSize);
    /* Spatial Kernel */
    afwMath::Kernel::SpatialFunctionPtr spatialKernelFunction(
        new afwMath::PolynomialFunction2<double>(spatialKernelOrder)
        );
    std::vector<afwMath::Kernel::SpatialFunctionPtr> spatialFunctionList;
    for (unsigned int i = 0; i < nBases; i++) {
        afwMath::Kernel::SpatialFunctionPtr spatialFunction(spatialKernelFunction->copy());
        spatialFunctionList.push_back(spatialFunction);
    }
    afwMath::LinearCombinationKernel::Ptr spatialKernel(
        new afwMath::LinearCombinationKernel(basisList, spatialFunctionList)
        );
    /* Set up some fake terms */
    std::vector<std::vector<double> > kCoeffs;
    kCoeffs.reserve(nBases);
    for (unsigned int i = 0, idx = 0; i < nBases; i++) {
        kCoeffs.push_back(std::vector<double>(spatialKernelFunction->getNParameters()));
        for (unsigned int j = 0; j < spatialKernelFunction->getNParameters(); j++, idx++) {
            kCoeffs[i][j] = std::sqrt(std::sqrt(std::sqrt(idx)));
        }
    }
    spatialKernel->setSpatialParameters(kCoeffs);
    
    unsigned int loc = 50;
    afwImage::MaskedImage<PixelT>::Ptr mimg1(
        new afwImage::MaskedImage<PixelT>(100,100)
        );
    *mimg1->at(loc, loc) = afwImage::MaskedImage<PixelT>::Pixel(1, 0x0, 1);
    afwImage::MaskedImage<PixelT>::Ptr mimg2(
        new afwImage::MaskedImage<PixelT>(mimg1->getDimensions())
        );
    afwMath::convolve(*mimg2, *mimg1, *spatialKernel, false);
    mimg1->writeFits("mimg1");
    mimg2->writeFits("mimg2");
    
    afwImage::Image<double> kImage(spatialKernel->getDimensions());
    (void)spatialKernel->computeImage(kImage, 
                                      false, 
                                      afwImage::indexToPosition(loc), 
                                      afwImage::indexToPosition(loc));
    kImage.writeFits("kernel.fits");
    
}
