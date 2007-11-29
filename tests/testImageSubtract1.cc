#include <fstream>

#include <boost/shared_ptr.hpp>

#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/FunctionLibrary.h>

#include <lsst/mwi/data/Citizen.h>
#include <lsst/mwi/exceptions/Exception.h>
#include <lsst/mwi/utils/Trace.h>
#include <lsst/mwi/policy/Policy.h>
#include <lsst/mwi/policy/paf/PAFParser.h>

#include <lsst/imageproc/ImageSubtract.h>

#define DEBUG_IO 1

using namespace std;

int main( int argc, char** argv )
{
    {
        lsst::mwi::utils::Trace::setVerbosity("lsst.imageproc", 4);
        
        typedef lsst::fw::maskPixelType MaskT;
        typedef float ImageT; // have to make sure this jibes with the input data!
        typedef double KernelT;
        typedef double FuncT;

        // Read in Policy
        ifstream is("tests/ImageSubtract_policy.paf");
        lsst::mwi::policy::Policy p;
        lsst::mwi::policy::paf::PAFParser pp(p);
        pp.parse(is);
        is.close();

        // Parse policy
        unsigned int kernelCols = p.getInt("kernelCols");
        unsigned int kernelRows = p.getInt("kernelRows");
        unsigned int kernelSpatialOrder = p.getInt("kernelSpatialOrder");
        unsigned int backgroundSpatialOrder = p.getInt("backgroundSpatialOrder");
        
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
        lsst::fw::MaskedImage<ImageT, MaskT> scienceMaskedImage;
        scienceMaskedImage.readFits(inputImage);

        lsst::fw::MaskedImage<ImageT, MaskT> templateMaskedImage;
        templateMaskedImage.readFits(inputImage);
        
        // Generate basis of delta functions for kernel
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > kernelBasisVec =
            lsst::imageproc::generateDeltaFunctionKernelSet<KernelT>(kernelCols, kernelRows);
        
        // Function for spatially varying kernel.  Make null here for this test.
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > kernelFunctionPtr(
            new lsst::fw::function::PolynomialFunction2<FuncT>(kernelSpatialOrder)
            );
        
        // Function for spatially varying background.  
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > backgroundFunctionPtr(
            new lsst::fw::function::PolynomialFunction2<FuncT>(backgroundSpatialOrder)
            );
        
        // Use hard-coded positions for now
        vector<lsst::detection::Footprint::PtrType> footprintList =
            lsst::imageproc::getCollectionOfMaskedImagesForPsfMatching();
        
        
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > kernelPtr =
            lsst::imageproc::computePsfMatchingKernelForMaskedImage(
                kernelFunctionPtr, backgroundFunctionPtr,
                templateMaskedImage, scienceMaskedImage, kernelBasisVec, footprintList, p);
    }
    
    if (Citizen::census(0) == 0) {
        cerr << "No leaks detected" << endl;
    } else {
        cerr << "Leaked memory blocks:" << endl;
        Citizen::census(cerr);
    } 
    
}
