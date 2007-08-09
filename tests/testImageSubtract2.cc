#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/FunctionLibrary.h>
#include <lsst/mwi/utils/Trace.h>
#include <lsst/mwi/exceptions/Exception.h>
#include <ImageSubtract.h>
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace lsst::fw;

int main( int argc, char** argv )
{
    lsst::mwi::utils::Trace::setDestination(cout);
    lsst::mwi::utils::Trace::setVerbosity(".", 5);
    
    typedef uint8 MaskT;
    typedef float ImageT; // have to make sure this jibes with the input data!
    typedef double KernelT;

    // Read input images
    if (argc < 2) {
        cout << "This program takes a single input image on the command line" << endl;
        cout << "  and uses it as the template image." << endl;
        cout << "  The science image is derived from the template convolved with a non-spatially varying Gaussian." << endl;
        cout << "  Your output kernel should the input Gaussian." << endl;
        cout << "  Basis function set is delta functions." << endl;
        cout << "  There is no spatial variation." << endl;
        exit(1);
    }
    string inputImage = argv[1];
    MaskedImage<ImageT,MaskT> scienceMaskedImage;
    try {
        scienceMaskedImage.readFits(inputImage);
    } catch (lsst::mwi::exceptions::Exception &e) {
        cerr << "Failed to open science image " << inputImage << ": " << e.what() << endl;
        return 1;
    }
    
    MaskedImage<ImageT,MaskT> templateMaskedImage;
    try {
        templateMaskedImage.readFits(inputImage);
    } catch (lsst::mwi::exceptions::Exception &e) {
        cerr << "Failed to open template image " << inputImage << ": " << e.what() << endl;
        return 1;
    }

    // Hardcoded
    unsigned int kernelRows = 7;
    unsigned int kernelCols = 7;
    
    // The kernel to convolve the template image with to yield the science image
    double sigmaX = 2.0;
    double sigmaY = 2.5;

    lsst::fw::Kernel<ImageT>::KernelFunctionPtrType gaussFuncPtr(
        new lsst::fw::function::GaussianFunction2<ImageT>(sigmaX, sigmaY));
    lsst::fw::AnalyticKernel<ImageT> gaussKernel(gaussFuncPtr, kernelCols, kernelRows);

    // Convolved science image
    lsst::mwi::utils::Trace("testImageSubtract2", 2, "Convolving input image for testing");
    const float threshold = 0.0;
    const int edgeMaskBit = 1;
    lsst::fw::MaskedImage<ImageT, MaskT> convolvedScienceMaskedImage =
        lsst::fw::kernel::convolve(scienceMaskedImage, gaussKernel, threshold, edgeMaskBit);
       
    convolvedScienceMaskedImage.writeFits( (boost::format("%s_test2") % inputImage).str() );

    // set up basis of delta functions for kernel
    vector<boost::shared_ptr<Kernel<KernelT> > > kernelBasisVec;
    lsst::imageproc::generateDeltaFunctionKernelSet(kernelRows, kernelCols, kernelBasisVec);

    lsst::imageproc::computePSFMatchingKernelForMaskedImage
        (convolvedScienceMaskedImage, templateMaskedImage, kernelBasisVec);

}
