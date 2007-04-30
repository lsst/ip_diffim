#include "lsst/fw/MaskedImage.h"
#include "lsst/fw/Trace.h"

#include <iostream>
using namespace std;
using namespace lsst::fw;

template <typename ImagePixelT, typename MaskPixelT> class setmaskPixProcFunc : public PixelProcessingFunc<ImagePixelT, MaskPixelT> 
{
public:
   typedef PixelLocator<ImagePixelT> ImageIteratorT;
   typedef PixelLocator<MaskPixelT> MaskIteratorT;
   typedef typename PixelChannelType<MaskPixelT>::type MaskChannelT;

   setmaskPixProcFunc(MaskedImage<ImagePixelT, MaskPixelT>& m) : PixelProcessingFunc<ImagePixelT, MaskPixelT>(m), initCount(0) {}

   void init() {
      PixelProcessingFunc<ImagePixelT, MaskPixelT>::_maskPtr->getPlaneBitMask("manual", bitSat);
      maskCount = 0;
      initCount++;
   }

   void operator ()(ImageIteratorT &i,MaskIteratorT &m) {
      if (*i > goodHighPixel) {
	 *m = *m | bitSat;
	 maskCount++;
      }
   }

   int getCount() { return maskCount; }

private:
   MaskChannelT bitSat;
   int maskCount;
   int initCount;
   int goodHighPixel = 20000;
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
   setMaskPixProcFunc<ImagePixelType, MaskPixelType> maskScienceFunc(cfhtScienceMaskedImage);
   maskScienceFunc.init();
   cfhtScienceMaskedImage.processPixels(maskScienceFunc);
   std::cout << "Set " << maskScienceFunc.getCount() << " mask bits" << std::endl;

   MaskedImage<ImagePixelType,MaskPixelType> cfhtTemplateMaskedImage;
   cfhtTemplateMaskedImage.readFits(argv[2]);
   cfhtTemplateMaskedImage.getMask()->addMaskPlane("manual");
   // Can I reuse maskScienceFunc?
   setMaskPixProcFunc<ImagePixelType, MaskPixelType> maskTemplateFunc(cfhtTemplateMaskedImage);
   maskTemplateFunc.init();
   cfhtTemplateMaskedImage.processPixels(maskTemplateFunc);
   std::cout << "Set " << maskTemplateFunc.getCount() << " mask bits" << std::endl;

   return 1;
}
