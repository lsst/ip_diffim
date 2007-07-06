#include "lsst/fw/MaskedImage.h"
#include "lsst/fw/Kernel.h"
#include "ImageSubtract.h"

using namespace std;
using namespace lsst::fw;

int main( int argc, char** argv )
{
    Trace::setDestination(cout);
    Trace::setVerbosity(".", 0);
    
    typedef uint8 MaskT;
    typedef float32 PixelT;
    typedef float32 KernelT;

    // Read input images
    string scienceInputImage = argv[1];
    MaskedImage<PixelT,MaskT> scienceMaskedImage;
    scienceMaskedImage.readFits(templateInputImage);
    
    string templateInputImage = argv[2];
    MaskedImage<PixelT,MaskT> templateMaskedImage;
    templateMaskedImage.readFits(templateInputImage);

    // set up basis of delta functions for kernel
    vector<Kernel<KernelT> > kernelVec;
    int kernelRows = 9;
    int kernelCols = 9;
    int colCtr = (kernelCols - 1) / 2;
    int rowCtr = (kernelRows - 1) / 2;
    for (unsigned row = 0; row < kernelRows; ++row) {
        int y = static_cast<int>(row) - rowCtr;
        
        for (unsigned col = 0; col < kernelCols; ++col) {
            int x = static_cast<int>(col) - colCtr;
            
            Kernel<kernelPixelType>::Function2PtrType kfuncPtr(
                new IntegerDeltaFunction2<kernelPixelType>(x, y)
                );
            Kernel<kernelPixelType> kernelPtr(
                new AnalyticKernel<kernelPixelType>(kfuncPtr, kernelCols, kernelRows)
                );
            kernelVec.push_back(kernelPtr);
        }
    }
    vector<PixelT> kernelParams(kernelRows * kernelCols);
    LinearCombinationKernel<KernelT> deltaFunctionKernelSet(kernelVec, kernelParams);

    // Currently does nothing
    getTemplateChunkExposureFromTemplateExposure();

    // This has some functionality!  Lets at least get it to compile.
    PSFMatchMaskedImages(scienceMaskedImage, templateMaskedImage, deltaFunctionKernelSet);

    // Currently does nothing
    subtractMatchedImage();

}
