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
 
#include <stdexcept>

#include <boost/shared_ptr.hpp>

#include <lsst/afw/image.h>
#include <lsst/afw/math.h>
#include <lsst/ip/diffim.h>

namespace afwImage = lsst::afw::image;
namespace afwMath  = lsst::afw::math;
namespace ipDiffim = lsst::ip::diffim;

void test() {
    typedef float PixelT;

    int cellSize   = 10;
    int fullSize   = 100;
    int stampSize  = 10;
    float coord    = 7.;
    
    boost::shared_ptr<afwImage::MaskedImage<PixelT> > tmi( 
        new afwImage::MaskedImage<PixelT>(stampSize, stampSize) 
        );
    boost::shared_ptr<afwImage::MaskedImage<PixelT> > smi( 
        new afwImage::MaskedImage<PixelT>(stampSize, stampSize) 
        );

    afwImage::BBox          bbox    = afwImage::BBox(afwImage::PointI(0, 0), fullSize, fullSize);
    afwMath::SpatialCellSet cellSet = afwMath::SpatialCellSet(bbox, cellSize, cellSize);

    ipDiffim::KernelCandidate<PixelT>::Ptr cand = ipDiffim::makeKernelCandidate(coord, coord, tmi, smi);
    std::cout << cand->getStatus() << std::endl;
    cellSet.insertCandidate(cand);
}

int main() {
    try {
        ::test();
    } catch (std::exception const &e) {
        std::clog << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    std::clog << "OK" << std::endl;
    return EXIT_SUCCESS;
}
