#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/FunctionLibrary.h>
#include <lsst/mwi/utils/Trace.h>
#include <lsst/mwi/exceptions/Exception.h>
#include <lsst/mwi/data/Citizen.h>
#include <lsst/imageproc/ImageSubtract.h>
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace lsst::fw;

int main( int argc, char** argv )
{
    {
        lsst::mwi::utils::Trace::setDestination(cout);
        lsst::mwi::utils::Trace::setVerbosity(".", 4);
        
        typedef uint8 MaskT;
        typedef float ImageT; // have to make sure this jibes with the input data!
        typedef double KernelT;
        typedef double FuncT;

        const KernelT threshold = 0.0;
        const int edgeMaskBit = 1;
        
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
        
        // Use hard-coded positions for now
        vector<lsst::detection::Footprint::PtrType> footprintVector;
        lsst::imageproc::getCollectionOfMaskedImagesForPSFMatching(footprintVector);
        
        
        lsst::imageproc::computePSFMatchingKernelForMaskedImage
            (templateMaskedImage, scienceMaskedImage, kernelBasisVec, footprintVector,
             kernelPtr, kernelFunctionPtr, backgroundFunctionPtr);
        
        lsst::fw::MaskedImage<ImageT, MaskT> convolvedTemplateMaskedImage =
            lsst::fw::kernel::convolve(templateMaskedImage, *kernelPtr, threshold, edgeMaskBit);
        
        scienceMaskedImage -= convolvedTemplateMaskedImage;
        scienceMaskedImage.writeFits( (boost::format("%s_diff1") % inputImage).str() );

        // TEST : the output kernel is a delta function.  The kernel coefficients of all bases other than the first (mean) are 0.
    }
    
    if (Citizen::census(0) == 0) {
        cerr << "No leaks detected" << endl;
    } else {
        cerr << "Leaked memory blocks:" << endl;
        Citizen::census(cerr);
    } 
    
}
