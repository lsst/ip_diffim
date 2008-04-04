// -*- lsst-c++ -*-
#include "lsst/afw/MaskedImage.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/daf/base/DataProperty.h"

using namespace std;
using namespace lsst::afw;

template <typename ImagePixelT, typename MaskPixelT> 
class synthesizeCfhtPixProcFunc : public PixelProcessingFunc<ImagePixelT, MaskPixelT> {
public:
    typedef PixelLocator<ImagePixelT> ImageIteratorT;
    typedef PixelLocator<MaskPixelT> MaskIteratorT;
    typedef typename PixelChannelType<MaskPixelT>::type MaskChannelT;
    
    synthesizeCfhtPixProcFunc(MaskedImage<ImagePixelT, MaskPixelT>& m) : PixelProcessingFunc<ImagePixelT, MaskPixelT>(m), initCount(0) {}
    
    void init() {
        PixelProcessingFunc<ImagePixelT, MaskPixelT>::_maskPtr->getPlaneBitMask("saturated", satBit);
        PixelProcessingFunc<ImagePixelT, MaskPixelT>::_maskPtr->getPlaneBitMask("zerovalued", badBit);
        DataPropertyPtrT metaDataPtr = PixelProcessingFunc<ImagePixelT, MaskPixelT>::_imagePtr->getMetaData();
        DataPropertyPtrT satPtr = metaDataPtr->find("MAXLIN");
        satValue = boost::any_cast<const int>(satPtr->getValue());

        satFrac = 0.9;
        satValue *= satFrac;

        badValue = 0;

        satCount = 0;
        badCount = 0;
        initCount++;
    }
    
    void operator ()(ImageIteratorT &i, MaskIteratorT &m) {
        if (*i >= satValue) {
            *m = *m | satBit;
            satCount++;
        }

        if (*i <= badValue) {
            *m = *m | badBit;
            badCount++;
        }
    }
    
    int getSatCount() { return satCount; }
    int getBadCount() { return badCount; }
    
private:
    MaskChannelT satBit;
    MaskChannelT badBit;

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
    
    typedef uint8 MaskPixelType;
    typedef float32 ImagePixelType;

    string scienceInputImage = argv[1];
    string scienceOutputImage = argv[2];
    string templateInputImage = argv[3];
    string templateOutputImage = argv[4];

    MaskedImage<ImagePixelType,MaskPixelType> cfhtScienceMaskedImage;
    cfhtScienceMaskedImage.readFits(scienceInputImage);
    cfhtScienceMaskedImage.getMask()->addMaskPlane("saturated");
    cfhtScienceMaskedImage.getMask()->addMaskPlane("zerovalued");
    cfhtScienceMaskedImage.setDefaultVariance();
    synthesizeCfhtPixProcFunc<ImagePixelType, MaskPixelType> maskScienceFunc(cfhtScienceMaskedImage);
    maskScienceFunc.init();
    cfhtScienceMaskedImage.processPixels(maskScienceFunc);
    cout << "Set " << maskScienceFunc.getSatCount() << " sat mask bits in " << scienceInputImage << endl;
    cout << "Set " << maskScienceFunc.getBadCount() << " bad mask bits in " << scienceInputImage << endl;
    cfhtScienceMaskedImage.writeFits(scienceOutputImage);
    
    MaskedImage<ImagePixelType,MaskPixelType> cfhtTemplateMaskedImage;
    cfhtTemplateMaskedImage.readFits(templateInputImage);
    cfhtTemplateMaskedImage.getMask()->addMaskPlane("saturated");
    cfhtTemplateMaskedImage.getMask()->addMaskPlane("zerovalued");
    cfhtTemplateMaskedImage.setDefaultVariance();
    // Can I reuse maskScienceFunc?
    synthesizeCfhtPixProcFunc<ImagePixelType, MaskPixelType> maskTemplateFunc(cfhtTemplateMaskedImage);
    maskTemplateFunc.init();
    cfhtTemplateMaskedImage.processPixels(maskTemplateFunc);
    cout << "Set " << maskTemplateFunc.getSatCount() << " sat mask bits in " << templateInputImage << endl;
    cout << "Set " << maskTemplateFunc.getBadCount() << " bad mask bits in " << templateInputImage << endl;
    cfhtTemplateMaskedImage.writeFits(templateOutputImage);
    
    return 1;
}
