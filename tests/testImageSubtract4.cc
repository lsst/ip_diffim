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
        lsst::mwi::utils::Trace::setVerbosity(".", 5);
        
        typedef uint8 MaskT;
        typedef float ImageT; // have to make sure this jibes with the input data!
        typedef double KernelT;
        typedef double FuncT;
        
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
        const KernelT threshold = 0.0;
        const int edgeMaskBit = 1;
        lsst::fw::MaskedImage<ImageT, MaskT> convolvedScienceMaskedImage =
            lsst::fw::kernel::convolve(scienceMaskedImage, gaussSpVarKernel, threshold, edgeMaskBit);
        
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
        unsigned int kernelSpatialOrder = polyOrder-1;
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > kernelFunctionPtr(
            new lsst::fw::function::PolynomialFunction2<FuncT>(kernelSpatialOrder)
            );
        
        // Function for spatially varying background.  
        unsigned int backgroundSpatialOrder = 0;
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > backgroundFunctionPtr(
            new lsst::fw::function::PolynomialFunction2<FuncT>(backgroundSpatialOrder)
            );
        
        lsst::imageproc::computePSFMatchingKernelForMaskedImage
            (templateMaskedImage, convolvedScienceMaskedImage, kernelBasisVec,
             kernelPtr, kernelFunctionPtr, backgroundFunctionPtr);
        
        lsst::fw::MaskedImage<ImageT, MaskT> convolvedTemplateMaskedImage =
            lsst::fw::kernel::convolve(templateMaskedImage, *kernelPtr, threshold, edgeMaskBit);
        
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
