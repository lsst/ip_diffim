// -*- lsst-c++ -*-
#include <iostream>
#include <string>

#include "lsst/daf/base/DataProperty.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image/MaskedImage.h"

using lsst::pex::logging::Trace;

template <typename ImagePixelT, typename MaskPixelT> 
class synthesizeCfhtPixProcFunc : public lsst::afw::image::PixelProcessingFunc<ImagePixelT, MaskPixelT> {
public:
    typedef lsst::afw::image::PixelLocator<ImagePixelT> ImageIteratorT;
    typedef lsst::afw::image::PixelLocator<MaskPixelT> MaskIteratorT;
    typedef typename lsst::afw::image::Mask<MaskPixelT>::MaskChannelT MaskChannelT;
    
    synthesizeCfhtPixProcFunc(lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT>& m) : PixelProcessingFunc<ImagePixelT, MaskPixelT>(m), initCount(0) {}
    
    void init() {
        satMask = PixelProcessingFunc<ImagePixelT, MaskPixelT>::_maskPtr->getPlaneBitMask("SAT");
        if (satMask == 0) {
            std::cout << "Warning: saturated mask plane not found" << std::endl;
        }
        badMask = PixelProcessingFunc<ImagePixelT, MaskPixelT>::_maskPtr->getPlaneBitMask("BAD");
        if (badMask == 0) {
            std::cout << "Warning: zerovalued mask plane not found" << std::endl;
        }
        lsst::daf::base::DataProperty::PtrType metaDataPtr = PixelProcessingFunc<ImagePixelT, MaskPixelT>::_imagePtr->getMetaData();
        lsst::daf::base::DataProperty::PtrType satPtr = metaDataPtr->findUnique("MAXLIN");
        satValue = boost::any_cast<const int>(satPtr->getValue());
        
        // Mask anything within 90% of saturation for now
        satFrac = 0.9;
        satValue *= satFrac;
        
        // Bad (zero-valued) pixels
        badValue = 0;

        satCount = 0;
        badCount = 0;
        initCount++;
    }
    
    void operator ()(ImageIteratorT &i, MaskIteratorT &m) {
        if (*i >= satValue) {
            *m = *m | satMask;
            satCount++;
        }
        
        if (*i <= badValue) {
            *m = *m | badMask;
            badCount++;
        }
    }
    
    int getSatCount() { return satCount; }
    int getBadCount() { return badCount; }
    
private:
    MaskChannelT satMask;
    MaskChannelT badMask;

    int initCount;

    int satCount;
    float satValue;
    float satFrac;

    int badCount;
    int badValue;
};



int main( int argc, char** argv )
{
    Trace::setDestination(std::cout);
    Trace::setVerbosity(".", 0);
    
    typedef lsst::afw::image::maskPixelType MaskPixelType;
    typedef float ImagePixelType;

    std::string inputImage = argv[1];
    std::string outputImage = argv[2];

    lsst::afw::image::MaskedImage<ImagePixelType,MaskPixelType> cfhtMaskedImage;
    cfhtMaskedImage.readFits(inputImage);
    cfhtMaskedImage.getMask()->addMaskPlane("BAD");
    cfhtMaskedImage.getMask()->addMaskPlane("SAT");
    cfhtMaskedImage.setDefaultVariance();
    synthesizeCfhtPixProcFunc<ImagePixelType, MaskPixelType> maskFunc(cfhtMaskedImage);
    maskFunc.init();
    cfhtMaskedImage.processPixels(maskFunc);
    std::cout << "Set " << maskFunc.getSatCount() << " sat mask bits in " << inputImage << std::endl;
    std::cout << "Set " << maskFunc.getBadCount() << " bad mask bits in " << inputImage << std::endl;
    cfhtMaskedImage.writeFits(outputImage);
   
    return 0;
}
