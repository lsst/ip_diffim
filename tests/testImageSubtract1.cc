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
        cout << "  and uses it as both the template and the science image." << endl;
        cout << "  Your output kernel should be a delta function." << endl;
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
    
    // set up basis of delta functions for kernel
    vector<boost::shared_ptr<Kernel<KernelT> > > kernelBasisVec;
    unsigned int kernelRows = 11;
    unsigned int kernelCols = 11;
    lsst::imageproc::generateDeltaFunctionKernelSet(kernelRows, kernelCols, kernelBasisVec);

    lsst::imageproc::computePSFMatchingKernelForMaskedImage
        (scienceMaskedImage, templateMaskedImage, kernelBasisVec);

}
