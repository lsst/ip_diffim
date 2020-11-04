/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#include <cmath>

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/geom.h"
#include "lsst/ip/diffim.h"

namespace geom = lsst::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
using namespace lsst::ip::diffim;

typedef float PixelT;

int main() {
    int kSize = 5;
    unsigned int nBases = kSize * kSize;
    int spatialKernelOrder = 2;

    lsst::afw::math::KernelList basisList = makeDeltaFunctionBasisList(kSize, kSize);
    /* Spatial Kernel */
    afwMath::Kernel::SpatialFunctionPtr spatialKernelFunction(
        new afwMath::PolynomialFunction2<double>(spatialKernelOrder)
        );
    std::vector<afwMath::Kernel::SpatialFunctionPtr> spatialFunctionList;
    for (unsigned int i = 0; i < nBases; i++) {
        afwMath::Kernel::SpatialFunctionPtr spatialFunction(spatialKernelFunction->clone());
        spatialFunctionList.push_back(spatialFunction);
    }
    std::shared_ptr<afwMath::LinearCombinationKernel> spatialKernel(
        new afwMath::LinearCombinationKernel(basisList, spatialFunctionList)
        );
    /* Set up some fake terms */
    std::vector<std::vector<double> > kCoeffs;
    kCoeffs.reserve(nBases);
    unsigned int idx = 0;
    for (unsigned int i = 0; i < nBases; i++) {
        kCoeffs.push_back(std::vector<double>(spatialKernelFunction->getNParameters()));
        for (unsigned int j = 0; j < spatialKernelFunction->getNParameters(); j++, idx++) {
            kCoeffs[i][j] = std::sqrt(std::sqrt(std::sqrt(idx)));
        }
    }
    spatialKernel->setSpatialParameters(kCoeffs);

    unsigned int loc = 50;
    std::shared_ptr<afwImage::MaskedImage<PixelT>> mimg1(
        new afwImage::MaskedImage<PixelT>(geom::Extent2I(100,100))
        );
    *mimg1->at(loc, loc) = afwImage::MaskedImage<PixelT>::Pixel(1, 0x0, 1);
    std::shared_ptr<afwImage::MaskedImage<PixelT>> mimg2(
        new afwImage::MaskedImage<PixelT>(mimg1->getDimensions())
        );
    afwMath::ConvolutionControl convolutionControl;
    convolutionControl.setDoNormalize(false);
    afwMath::convolve(*mimg2, *mimg1, *spatialKernel, convolutionControl);
    mimg1->writeFits("mimg1");
    mimg2->writeFits("mimg2");

    afwImage::Image<double> kImage(spatialKernel->getDimensions());
    (void)spatialKernel->computeImage(kImage,
                                      false,
                                      afwImage::indexToPosition(loc),
                                      afwImage::indexToPosition(loc));
    kImage.writeFits("kernel.fits");

}
