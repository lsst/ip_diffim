// -*- lsst-c++ -*-
#include "lsst/fw/MaskedImage.h"
#include "lsst/mwi/utils/Trace.h"
#include "lsst/mwi/data/DataProperty.h"

using namespace std;
using namespace lsst::fw;

template <typename ImagePixelT, typename MaskPixelT> 
class synthesizeCfhtPixProcFunc : public PixelProcessingFunc<ImagePixelT, MaskPixelT> {
public:
    typedef PixelLocator<ImagePixelT> ImageIteratorT;
    typedef PixelLocator<MaskPixelT> MaskIteratorT;
    typedef typename PixelChannelType<MaskPixelT>::type MaskChannelT;
    
    synthesizeCfhtPixProcFunc(MaskedImage<ImagePixelT, MaskPixelT>& m) : PixelProcessingFunc<ImagePixelT, MaskPixelT>(m), initCount(0) {}
    
    void init() {
        satMask = PixelProcessingFunc<ImagePixelT, MaskPixelT>::_maskPtr->getPlaneBitMask("SAT");
        if (satMask == 0) {
            cout << "Warning: saturated mask plane not found" << endl;
        }
        badMask = PixelProcessingFunc<ImagePixelT, MaskPixelT>::_maskPtr->getPlaneBitMask("BAD");
        if (badMask == 0) {
            cout << "Warning: zerovalued mask plane not found" << endl;
        }
        lsst::mwi::data::DataProperty::PtrType metaDataPtr = PixelProcessingFunc<ImagePixelT, MaskPixelT>::_imagePtr->getMetaData();
        lsst::mwi::data::DataProperty::PtrType satPtr = metaDataPtr->findUnique("MAXLIN");
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
    Trace::setDestination(cout);
    Trace::setVerbosity(".", 0);
    
    typedef lsst::fw::maskPixelType MaskPixelType;
    typedef float32 ImagePixelType;

    string inputImage = argv[1];
    string outputImage = argv[2];

    MaskedImage<ImagePixelType,MaskPixelType> cfhtMaskedImage;
    cfhtMaskedImage.readFits(inputImage);
    cfhtMaskedImage.getMask()->addMaskPlane("BAD");
    cfhtMaskedImage.getMask()->addMaskPlane("SAT");
    cfhtMaskedImage.setDefaultVariance();
    synthesizeCfhtPixProcFunc<ImagePixelType, MaskPixelType> maskFunc(cfhtMaskedImage);
    maskFunc.init();
    cfhtMaskedImage.processPixels(maskFunc);
    cout << "Set " << maskFunc.getSatCount() << " sat mask bits in " << inputImage << endl;
    cout << "Set " << maskFunc.getBadCount() << " bad mask bits in " << inputImage << endl;
    cfhtMaskedImage.writeFits(outputImage);
   
    return 0;
}
