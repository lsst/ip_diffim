#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/FunctionLibrary.h>
#include <lsst/fw/PixelAccessors.h>
#include <lsst/mwi/utils/Trace.h>
#include <lsst/mwi/exceptions/Exception.h>
#include <lsst/mwi/data/Citizen.h>
#include <lsst/imageproc/ImageSubtract.h>
#include <lsst/detection/Footprint.h>
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace lsst::fw;

typedef uint8 MaskT;
typedef float ImageT;
typedef double KernelT;
typedef double FuncT;

int main( int argc, char** argv )
{
    {
        lsst::mwi::utils::Trace::setDestination(cout);
        lsst::mwi::utils::Trace::setVerbosity(".", 4);
        
        // Read input images
        if (argc < 2) {
            cout << "This program takes 2 input images on the command line" << endl;
            cout << "  and finds its own detection." << endl;
            cout << "  Basis function set is delta functions." << endl;
            cout << "  There is spatial variation." << endl;
            exit(1);
        }
        string templateImage = argv[1];
        MaskedImage<ImageT,MaskT> templateMaskedImage;
        try {
            templateMaskedImage.readFits(templateImage);
        } catch (lsst::mwi::exceptions::Exception &e) {
            cerr << "Failed to open template image " << templateImage << ": " << e.what() << endl;
            return 1;
        }

        string scienceImage = argv[2];
        MaskedImage<ImageT,MaskT> scienceMaskedImage;
        try {
            scienceMaskedImage.readFits(scienceImage);
        } catch (lsst::mwi::exceptions::Exception &e) {
            cerr << "Failed to open science image " << scienceImage << ": " << e.what() << endl;
            return 1;
        }
       
        // Generate basis of delta functions for kernel
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > kernelBasisVec;
        unsigned int kernelRows = 7;
        unsigned int kernelCols = 7;
        lsst::imageproc::generateDeltaFunctionKernelSet(kernelRows, kernelCols, kernelBasisVec);
        
        // Output kernel
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > kernelPtr(
            new lsst::fw::LinearCombinationKernel<KernelT>
            );
        
        // Function for spatially varying kernel.  Make null here for this test.
        unsigned int kernelSpatialOrder = 0;
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > kernelFunctionPtr(
            new lsst::fw::function::PolynomialFunction2<FuncT>(kernelSpatialOrder)
            );
        
        // Function for spatially varying background.  
        unsigned int backgroundSpatialOrder = 0;
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > backgroundFunctionPtr(
            new lsst::fw::function::PolynomialFunction2<FuncT>(backgroundSpatialOrder)
            );

        // Now you let the code find the peaks!
        lsst::imageproc::computePSFMatchingKernelForMaskedImage
            (templateMaskedImage, scienceMaskedImage, kernelBasisVec, 
             kernelPtr, kernelFunctionPtr, backgroundFunctionPtr);

    }
    
    if (Citizen::census(0) == 0) {
        cerr << "No leaks detected" << endl;
    } else {
        cerr << "Leaked memory blocks:" << endl;
        Citizen::census(cerr);
    } 
}
