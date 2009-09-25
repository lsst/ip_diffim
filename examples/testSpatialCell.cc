#include <stdexcept>

#include <boost/shared_ptr.hpp>

#include <lsst/afw/image.h>
#include <lsst/afw/math.h>
#include <lsst/ip/diffim.h>

namespace afwImage = lsst::afw::image;
namespace afwMath  = lsst::afw::math;
namespace ipDiffim = lsst::ip::diffim;

/*
 * Make this a subroutine so that locals go out of scope as part of test
 * of memory management
 */
void test() {
    typedef float PixelT;

    int cellSize   = 10;
    int fullSize   = 100;
    int stampSize  = 10;
    float coord    = 7.;
    
    boost::shared_ptr<afwImage::MaskedImage<PixelT> > tmi( new afwImage::MaskedImage<PixelT>(stampSize, stampSize) );
    boost::shared_ptr<afwImage::MaskedImage<PixelT> > smi( new afwImage::MaskedImage<PixelT>(stampSize, stampSize) );
    *tmi->getImage() = 0;
    *tmi->getMask() = 0x0;
    *tmi->getVariance() = 0;
    *smi->getImage() = 0;
    *smi->getMask() = 0x0;
    *smi->getVariance() = 0;


    afwImage::BBox          bbox    = afwImage::BBox(afwImage::PointI(0, 0), fullSize, fullSize);
    afwMath::SpatialCellSet cellSet = afwMath::SpatialCellSet(bbox, cellSize);

    ipDiffim::KernelCandidate<PixelT>::Ptr cand = ipDiffim::makeKernelCandidate(coord, coord, tmi, smi);
    cellSet.insertCandidate(cand);
}

int main() {
    try {
        test();
    } catch (std::exception const &e) {
        std::clog << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    std::clog << "OK" << std::endl;
    return EXIT_SUCCESS;
}
