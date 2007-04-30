#include "lsst/fw/MaskedImage.h"
#include "lsst/fw/Trace.h"

#include <iostream>
using namespace std;
using namespace lsst::fw;

template <typename ImagePixelT, typename MaskPixelT> class synthesizeCfhtPixProcFunc : public PixelProcessingFunc<ImagePixelT, MaskPixelT> 
{
public:
   typedef PixelLocator<ImagePixelT> ImageIteratorT;
   typedef PixelLocator<VariancePixelT> VarianceIteratorT;
   typedef PixelLocator<MaskPixelT> MaskIteratorT;
   typedef typename PixelChannelType<MaskPixelT>::type MaskChannelT;

   setmaskPixProcFunc(MaskedImage<ImagePixelT, VariancePixelT, MaskPixelT>& m) : PixelProcessingFunc<ImagePixelT, VariancePixelT, MaskPixelT>(m), initCount(0) {}

   void init() {
      PixelProcessingFunc<ImagePixelT, VariancePixelT, MaskPixelT>::_maskPtr->getPlaneBitMask("manual", bitSat);
      DataProperty::DataPropertyPtrT metaDataPtr = PixelProcessingFunc<ImagePixelT, VariancePixelT, MaskPixelT>::_imagePtr->getMetaData();
      DataProperty::DataPropertyPtrT satPtr = metaDataPtr->find(boost::regex("MAXLIN"));
      DataProperty::DataPropertyPtrT gainPtr = metaDataPtr->find(boost::regex("GAIN"));
      DataProperty::DataPropertyPtrT rdnoisePtr = metaDataPtr->find(boost::regex("RDNOISE"));
      maskCount = 0;
      initCount++;
   }

   void operator ()(ImageIteratorT &i, VarianceIterator &v, MaskIteratorT &m) {
      if (*i > satValue) {
	 *m = *m | bitSat;
	 maskCount++;
      }
      // synthesize variance :
      // noise_e**2   = pixel_e + rdnoise**2
      // pixel_e      = pixel_adu * gain
      // noise_adu**2 = (noise_e / gain)**2
      //              = pixel_e / gain**2 + (rdnoise/gain)**2
      //              = pixel_adu / gain  + (rdnoise/gain)**2
      *v = *i * igainValue + rdnoiseOverGainSquared;
   }

   int getCount() { return maskCount; }

private:
   MaskChannelT bitSat;
   int maskCount;
   int initCount;
   float satValue;
   float igainValue;   // inverse gain
   float rdnoiseValue; 
   float rdnoiseOverGainSquared;
};



int main( int argc, char** argv )
{
   Trace::setDestination(std::cout);
   Trace::setVerbosity(".", 0);

   typedef uint8 MaskPixelType;
   typedef float32 ImagePixelType;

   MaskedImage<ImagePixelType,MaskPixelType> cfhtScienceMaskedImage;
   cfhtScienceMaskedImage.readFits(argv[1]);
   cfhtScienceMaskedImage.getMask()->addMaskPlane("manual");
   synthesizeCfhtPixProcFunc<ImagePixelType, MaskPixelType> maskScienceFunc(cfhtScienceMaskedImage);
   maskScienceFunc.init();
   cfhtScienceMaskedImage.processPixels(maskScienceFunc);
   std::cout << "Set " << maskScienceFunc.getCount() << " mask bits" << std::endl;

   MaskedImage<ImagePixelType,MaskPixelType> cfhtTemplateMaskedImage;
   cfhtTemplateMaskedImage.readFits(argv[2]);
   cfhtTemplateMaskedImage.getMask()->addMaskPlane("manual");
   // Can I reuse maskScienceFunc?
   synthesizeCfhtPixProcFunc<ImagePixelType, MaskPixelType> maskTemplateFunc(cfhtTemplateMaskedImage);
   maskTemplateFunc.init();
   cfhtTemplateMaskedImage.processPixels(maskTemplateFunc);
   std::cout << "Set " << maskTemplateFunc.getCount() << " mask bits" << std::endl;

   return 1;
}
