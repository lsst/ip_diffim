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
    unsigned kernelRows = 3;
    unsigned kernelCols = 3;
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

    if (0) {
        
        // And, lets convolve the science image with a gaussian to see what we get back out.
        double sigmaX = 2.0;
        double sigmaY = 2.5;
        lsst::fw::Kernel<ImageT>::KernelFunctionPtrType gaussFuncPtr(
            new lsst::fw::GaussianFunction2<ImageT>(sigmaX, sigmaY));
        lsst::fw::AnalyticKernel<ImageT> gaussKernel(gaussFuncPtr, kernelCols, kernelRows);
        cout << boost::format("Convolving science image with a gaussian Kernel with sigmaX=%.1f, sigmaY=%.1f\n\n") % sigmaX % sigmaY;
        
        const float threshold = 0.0;
        const int edgeMaskBit = -1;
        lsst::fw::MaskedImage<ImageT, MaskT> convolvedScienceMaskedImage = 
            lsst::fw::kernel::convolve(scienceMaskedImage, gaussKernel, threshold, vw::NoEdgeExtension(), edgeMaskBit);
        cout << "in " << scienceMaskedImage.getRows() << " " << scienceMaskedImage.getCols() << endl;
        cout << "out " << convolvedScienceMaskedImage.getRows() << " " << convolvedScienceMaskedImage.getCols() << endl;
        convolvedScienceMaskedImage.writeFits( (boost::format("convolvedS")).str() );
        
        ////////////////
        
        lsst::fw::Kernel<ImageT>::KernelFunctionPtrType deltaFuncPtr(new IntegerDeltaFunction2<ImageT>(colCtr, rowCtr));
        lsst::fw::AnalyticKernel<ImageT> deltaKernel(deltaFuncPtr, kernelCols, kernelRows);
        cout << boost::format("Convolving template image with a delta function\n\n") << endl;
        lsst::fw::MaskedImage<ImageT, MaskT> convolvedTemplateMaskedImage = 
            lsst::fw::kernel::convolve(templateMaskedImage, deltaKernel, threshold, vw::NoEdgeExtension(), edgeMaskBit);
        cout << "in " << templateMaskedImage.getRows() << " " << templateMaskedImage.getCols() << endl;
        cout << "out " << convolvedTemplateMaskedImage.getRows() << " " << convolvedTemplateMaskedImage.getCols() << endl;
        convolvedTemplateMaskedImage.writeFits( (boost::format("convolvedT")).str() );
        
        ////////////////
        
        
        
        lsst::imageproc::computePSFMatchingKernelForMaskedImage
            (convolvedScienceMaskedImage, convolvedTemplateMaskedImage, kernelBasisVec);
    }
    else {
        lsst::imageproc::computePSFMatchingKernelForMaskedImage
            (scienceMaskedImage, templateMaskedImage, kernelBasisVec);
    }

    // Currently does nothing
    //subtractMatchedImage();
    cout << endl;

}
