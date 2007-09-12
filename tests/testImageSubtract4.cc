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
using namespace lsst::fw;

int main( int argc, char** argv )
{
    {
        lsst::mwi::utils::Trace::setDestination(cout);
        lsst::mwi::utils::Trace::setVerbosity(".", 5);
        
        typedef uint8 MaskT;
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
        KernelT convolveThreshold = static_cast<KernelT>(p.getDouble("convolveThreshold"));
        int edgeMaskBit = p.getInt("edgeMaskBit");
        unsigned int kernelRows = p.getInt("kernelRows");
        unsigned int kernelCols = p.getInt("kernelCols");
        unsigned int kernelSpatialOrder = p.getInt("kernelSpatialOrder");
        unsigned int backgroundSpatialOrder = p.getInt("backgroundSpatialOrder");
        
        // Read input images
        if (argc < 2) {
            cout << "This program takes a single input image on the command line" << endl;
            cout << "  and uses it as the template image." << endl;
            cout << "  The science image is derived from the template convolved with a spatially varying Gaussian." << endl;
            cout << "  Your output kernel should the input Gaussian." << endl;
            cout << "  Basis function set is delta functions." << endl;
            cout << "  There is no spatial variation." << endl;
            exit(1);
        }
        string inputImage = argv[1];
        MaskedImage<ImageT,MaskT> scienceMaskedImage;
        try {
            scienceMaskedImage.readFits(inputImage);
        } catch (lsst::mwi::exceptions::ExceptionStack &e) {
            cerr << "Failed to open science image " << inputImage << ": " << e.what() << endl;
            return 1;
        }
        
        MaskedImage<ImageT,MaskT> templateMaskedImage;
        try {
            templateMaskedImage.readFits(inputImage);
        } catch (lsst::mwi::exceptions::ExceptionStack &e) {
            cerr << "Failed to open template image " << inputImage << ": " << e.what() << endl;
            return 1;
        }
        
        // The kernel to convolve the template image with to yield the science image
        double sigmaX = 1.0;
        double sigmaY = 2.0;
        unsigned int polyOrder = 2;
        double minSigma = 0.1;
        double maxSigma = 3.0;
        
        lsst::fw::Kernel<KernelT>::KernelFunctionPtrType gaussFuncPtr(
            new lsst::fw::function::GaussianFunction2<FuncT>(sigmaX, sigmaY));
        lsst::fw::Kernel<KernelT>::SpatialFunctionPtrType polyFuncPtr(
            new lsst::fw::function::PolynomialFunction2<double>(polyOrder));
        lsst::fw::AnalyticKernel<KernelT> gaussSpVarKernel(
            gaussFuncPtr, kernelCols, kernelRows, polyFuncPtr);
        
        // Get copy of spatial parameters (all zeros), set and feed back to the kernel
        vector<vector<double> > polyParams = gaussSpVarKernel.getSpatialParameters();
        // Set spatial parameters for kernel parameter 0
        polyParams[0][0] = minSigma;
        polyParams[0][1] = (maxSigma - minSigma) / static_cast<double>(scienceMaskedImage.getCols());
        polyParams[0][2] = 0.0;
        // Set spatial function parameters for kernel parameter 1
        polyParams[1][0] = minSigma;
        polyParams[1][1] = 0.0;
        polyParams[1][2] = (maxSigma - minSigma) / static_cast<double>(scienceMaskedImage.getRows());
        gaussSpVarKernel.setSpatialParameters(polyParams);
        
        // Convolved science image
        lsst::mwi::utils::Trace("testImageSubtract4", 2, "Convolving input image for testing");
        lsst::fw::MaskedImage<ImageT, MaskT> convolvedScienceMaskedImage =
            lsst::fw::kernel::convolve(scienceMaskedImage, gaussSpVarKernel, convolveThreshold, edgeMaskBit);
        
        convolvedScienceMaskedImage.writeFits( (boost::format("%s_test4") % inputImage).str() );
        
        // Generate basis of delta functions for kernel
        vector<boost::shared_ptr<Kernel<KernelT> > > kernelBasisVec;
        lsst::imageproc::generateDeltaFunctionKernelSet(kernelRows, kernelCols, kernelBasisVec);
        
        
        // Output kernel
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > kernelPtr(
            new lsst::fw::LinearCombinationKernel<KernelT>
            );
        
        // Function for spatially varying kernel. 
        // Give it less power than the Fcn we convolved the image with, see what happens
        // OVERRIDE POLICY HERE
        kernelSpatialOrder = polyOrder;
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > kernelFunctionPtr(
            new lsst::fw::function::PolynomialFunction2<FuncT>(kernelSpatialOrder)
            );
        
        // Function for spatially varying background.  
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > backgroundFunctionPtr(
            new lsst::fw::function::PolynomialFunction2<FuncT>(backgroundSpatialOrder)
            );

        // Use hard-coded positions for now
        vector<lsst::detection::Footprint::PtrType> footprintList;
        lsst::imageproc::getCollectionOfMaskedImagesForPsfMatching(footprintList);

        lsst::imageproc::computePsfMatchingKernelForMaskedImage
            (templateMaskedImage, convolvedScienceMaskedImage, kernelBasisVec, footprintList,
             kernelPtr, kernelFunctionPtr, backgroundFunctionPtr, p);
        
        lsst::fw::MaskedImage<ImageT, MaskT> convolvedTemplateMaskedImage =
            lsst::fw::kernel::convolve(templateMaskedImage, *kernelPtr, convolveThreshold, edgeMaskBit);
        
        convolvedScienceMaskedImage -= convolvedTemplateMaskedImage;
        convolvedScienceMaskedImage.writeFits( (boost::format("%s_diff4") % inputImage).str() );

        // TEST : The output kernel is has same spatial variation as we input above
        //      : The background function has no spatial variation and value equal to 0
        
    }
    
    if (Citizen::census(0) == 0) {
        cerr << "No leaks detected" << endl;
    } else {
        cerr << "Leaked memory blocks:" << endl;
        Citizen::census(cerr);
    } 
    
}
