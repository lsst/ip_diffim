// -*- lsst-c++ -*-
#include "lsst/fw/MaskedImage.h"
#include "lsst/fw/Trace.h"
#include "lsst/fw/DiskImageResourceFITS.h"

#include <iostream>
using namespace lsst;

template <typename ImagePixelT, typename MaskPixelT> class synthesizeCfhtPixProcFunc : public PixelProcessingFunc<ImagePixelT, MaskPixelT> 
{
public:
    typedef PixelLocator<ImagePixelT> ImageIteratorT;
    //typedef PixelLocator<VariancePixelT> VarianceIteratorT;
    typedef PixelLocator<MaskPixelT> MaskIteratorT;
    typedef typename PixelChannelType<MaskPixelT>::type MaskChannelT;
    
    //synthesizeCfhtPixProcFunc(MaskedImage<ImagePixelT, VariancePixelT, MaskPixelT>& m) : PixelProcessingFunc<ImagePixelT, VariancePixelT, MaskPixelT>(m), initCount(0) {}
    synthesizeCfhtPixProcFunc(MaskedImage<ImagePixelT, MaskPixelT>& m) : PixelProcessingFunc<ImagePixelT, MaskPixelT>(m), initCount(0) {}
    
    void init() {
        //PixelProcessingFunc<ImagePixelT, VariancePixelT, MaskPixelT>::_maskPtr->getPlaneBitMask("manual", bitSat);
        PixelProcessingFunc<ImagePixelT, MaskPixelT>::_maskPtr->getPlaneBitMask("manual", bitSat);
        //DataProperty::DataPropertyPtrT metaDataPtr = PixelProcessingFunc<ImagePixelT, VariancePixelT, MaskPixelT>::_imagePtr->getMetaData();
        DataProperty::DataPropertyPtrT metaDataPtr = PixelProcessingFunc<ImagePixelT, MaskPixelT>::_imagePtr->getMetaData();
        DataProperty::DataPropertyPtrT satPtr = metaDataPtr->find(boost::regex("MAXLIN"));
        satValue = boost::any_cast<const float>(satPtr->getValue());
        maskCount = 0;
        initCount++;
    }
    
    //void operator ()(ImageIteratorT &i, VarianceIterator &v, MaskIteratorT &m) {
    void operator ()(ImageIteratorT &i, MaskIteratorT &m) {
        if (*i > satValue) {
            *m = *m | bitSat;
            maskCount++;
        }
    }
    
    int getCount() { return maskCount; }
    
private:
    MaskChannelT bitSat;
    int maskCount;
    int initCount;
    float satValue;
};



int main( int argc, char** argv )
{
    fw::Trace::setDestination(std::cout);
    fw::Trace::setVerbosity(".", 0);
    
    typedef uint8 MaskPixelType;
    typedef float32 ImagePixelType;
    
    MaskedImage<ImagePixelType,MaskPixelType> cfhtScienceMaskedImage;
    cfhtScienceMaskedImage.readFits(argv[1]);
    cfhtScienceMaskedImage.getMask()->addMaskPlane("manual");
    cfhtScienceMaskedImage.setDefaultVariance();
    synthesizeCfhtPixProcFunc<ImagePixelType, MaskPixelType> maskScienceFunc(cfhtScienceMaskedImage);
    maskScienceFunc.init();
    cfhtScienceMaskedImage.processPixels(maskScienceFunc);
    std::cout << "Set " << maskScienceFunc.getCount() << " mask bits" << std::endl;
    
    MaskedImage<ImagePixelType,MaskPixelType> cfhtTemplateMaskedImage;
    cfhtTemplateMaskedImage.readFits(argv[2]);
    cfhtTemplateMaskedImage.getMask()->addMaskPlane("manual");
    cfhtTemplateMaskedImage.setDefaultVariance();
    // Can I reuse maskScienceFunc?
    synthesizeCfhtPixProcFunc<ImagePixelType, MaskPixelType> maskTemplateFunc(cfhtTemplateMaskedImage);
    maskTemplateFunc.init();
    cfhtTemplateMaskedImage.processPixels(maskTemplateFunc);
    std::cout << "Set " << maskTemplateFunc.getCount() << " mask bits" << std::endl;
    
    return 1;
}
