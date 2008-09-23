#include "lsst/detection/Footprint.h"
#include "lsst/ip/diffim/ImageSubtract.h"

using namespace lsst::afw::image;
using namespace lsst::ip::diffim;
typedef float ImagePixelT;

int main() {
    typedef DifferenceImageFootprintInformation<ImagePixelT, lsst::afw::image::maskPixelType> Difi;

    boost::shared_ptr<MaskedImage<ImagePixelT, lsst::afw::image::maskPixelType> > miPtr (
        new MaskedImage<ImagePixelT, lsst::afw::image::maskPixelType>(100,100)
        );

    vw::BBox2i const region(10,10,5,5);
    lsst::detection::Footprint::PtrType fpPtr(
        new lsst::detection::Footprint(region)
        ); 
    
    Difi difi(fpPtr, miPtr, miPtr);
    Difi::Ptr difiPtr (new Difi(fpPtr, miPtr, miPtr));

    Difi::DifiList difiList;
    difiList.push_back(difiPtr);

    Difi::DifiList difiList2 = getGoodFootprints(difiList);

    return 0;
}

