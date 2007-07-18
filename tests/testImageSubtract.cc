#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/FunctionLibrary.h>
#include <ImageSubtract.h>
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace lsst::fw;

int main( int argc, char** argv )
{
    Trace::setDestination(cout);
    Trace::setVerbosity(".", 5);
    
    typedef uint8 MaskT;
    typedef float ImageT; // have to make sure this jibes with the input data!
    typedef double KernelT;

    // Read input images
    string scienceInputImage = argv[1];
    MaskedImage<ImageT,MaskT> scienceMaskedImage;
    try {
        scienceMaskedImage.readFits(scienceInputImage);
    } catch (lsst::fw::Exception &e) {
        cerr << "Failed to open " << scienceInputImage << ": " << e.what() << endl;
        return 1;
    }
    
    string templateInputImage = argv[2];
    MaskedImage<ImageT,MaskT> templateMaskedImage;
    try {
        templateMaskedImage.readFits(templateInputImage);
    } catch (lsst::fw::Exception &e) {
        cerr << "Failed to open " << templateInputImage << ": " << e.what() << endl;
        return 1;
    }
    
    // set up basis of delta functions for kernel
    // One way
    //vector<Kernel<KernelT> > kernelVec;
    // Another way with pointers
    vector<boost::shared_ptr<Kernel<KernelT> > > kernelBasisVec;
    unsigned kernelRows = 9;
    unsigned kernelCols = 9;
    int colCtr = (kernelCols - 1) / 2;
    int rowCtr = (kernelRows - 1) / 2;
    for (unsigned row = 0; row < kernelRows; ++row) {
        int y = static_cast<int>(row) - rowCtr;
        
        for (unsigned col = 0; col < kernelCols; ++col) {
            int x = static_cast<int>(col) - colCtr;
            
            Kernel<KernelT>::KernelFunctionPtrType kfuncPtr(
                new IntegerDeltaFunction2<KernelT>(x, y)
                );
 
            // One way
            //AnalyticKernel<KernelT> kernel(kfuncPtr, kernelCols, kernelRows);
            // Another way with pointers
            boost::shared_ptr<Kernel<KernelT> > kernelPtr(
                new AnalyticKernel<KernelT>(kfuncPtr, kernelCols, kernelRows)
            );
 
            kernelBasisVec.push_back(kernelPtr);
        }
    }

    // This has some functionality!  Lets at least get it to compile.
    lsst::imageproc::computePSFMatchingKernelForMaskedImage
        (scienceMaskedImage, templateMaskedImage, kernelBasisVec);

    // Currently does nothing
    //subtractMatchedImage();
    cout << endl;

}
