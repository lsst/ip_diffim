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
        double sigmaX = 1.0;
        double sigmaY = 2.0;
        
        lsst::fw::Kernel<KernelT>::KernelFunctionPtrType gaussFuncPtr(
            new lsst::fw::function::GaussianFunction2<FuncT>(sigmaX, sigmaY));
        lsst::fw::AnalyticKernel<KernelT> gaussKernel(gaussFuncPtr, kernelCols, kernelRows);
        
        // Convolved science image
        lsst::mwi::utils::Trace("testImageSubtract3", 2, "Convolving input image for testing");
        const KernelT threshold = 0.0;
        const int edgeMaskBit = 1;
        lsst::fw::MaskedImage<ImageT, MaskT> convolvedScienceMaskedImage =
            lsst::fw::kernel::convolve(scienceMaskedImage, gaussKernel, threshold, edgeMaskBit);
        
        convolvedScienceMaskedImage.writeFits( (boost::format("%s_test3") % inputImage).str() );
        
        // Generate basis of delta functions for kernel
        vector<boost::shared_ptr<Kernel<KernelT> > > kernelBasisVec;
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
        
        // Hack some positions in for /lsst/becker/lsst_devel/DC2/fw/tests/data/871034p_1_MI_img.fits
        vector<lsst::fw::Source> sourceCollection;
        double radius = 20;
        lsst::fw::Source src1(0, 78.654, 3573.945, radius, radius);
        lsst::fw::Source src2(1, 341.149, 2753.536, radius, radius);
        lsst::fw::Source src3(2, 353.237, 2755.959, radius, radius);
        lsst::fw::Source src4(3, 367.756, 3827.671, radius, radius);
        lsst::fw::Source src5(4, 381.062, 3212.948, radius, radius);
        lsst::fw::Source src6(5, 404.433, 573.462, radius, radius);
        lsst::fw::Source src7(6, 420.967, 3306.310, radius, radius);
        lsst::fw::Source src8(7, 518.953, 2035.000, radius, radius);
        lsst::fw::Source src9(8, 546.657, 285.079, radius, radius); // This one is a CR!
        //lsst::fw::Source src10(10, 779.549, 1672.990, radius, radius);
        //lsst::fw::Source src11(11, 1010.618, 2375.691, radius, radius);
        //lsst::fw::Source src12(12, 1219.023, 1683.485, radius, radius);
        //lsst::fw::Source src13(13, 1457.960, 2282.024, radius, radius);
        //lsst::fw::Source src14(14, 1588.583, 3536.200, radius, radius);
        //lsst::fw::Source src15(15, 1604.816, 4070.769, radius, radius);
        //lsst::fw::Source src16(16, 1609.222, 4071.211, radius, radius);
        //lsst::fw::Source src17(17, 1686.953, 1880.928, radius, radius);
        //lsst::fw::Source src18(18, 1698.308, 3860.842, radius, radius);
        //lsst::fw::Source src19(19, 1709.934, 3861.217, radius, radius);
        //lsst::fw::Source src20(20, 1737.637, 3139.729, radius, radius);
        //lsst::fw::Source src21(21, 1794.293, 2023.711, radius, radius);
        //lsst::fw::Source src22(22, 1799.249, 733.596, radius, radius);
        //lsst::fw::Source src23(23, 1959.672, 4232.035, radius, radius);
        sourceCollection.push_back(src1);
        sourceCollection.push_back(src2);
        sourceCollection.push_back(src3);
        sourceCollection.push_back(src4);
        sourceCollection.push_back(src5);
        sourceCollection.push_back(src6);
        sourceCollection.push_back(src7);
        sourceCollection.push_back(src8);
        sourceCollection.push_back(src9);
        //sourceCollection.push_back(src10);
        //sourceCollection.push_back(src11);
        //sourceCollection.push_back(src12);
        //sourceCollection.push_back(src13);
        //sourceCollection.push_back(src14);
        //sourceCollection.push_back(src15);
        //sourceCollection.push_back(src16);
        //sourceCollection.push_back(src17);
        //sourceCollection.push_back(src18);
        //sourceCollection.push_back(src19);
        //sourceCollection.push_back(src20);
        //sourceCollection.push_back(src21);
        //sourceCollection.push_back(src22);
        //sourceCollection.push_back(src23);
        
        lsst::imageproc::computePSFMatchingKernelForMaskedImage
            (templateMaskedImage, convolvedScienceMaskedImage, kernelBasisVec, 
             sourceCollection, kernelPtr, kernelFunctionPtr, backgroundFunctionPtr);
        
        lsst::fw::MaskedImage<ImageT, MaskT> convolvedTemplateMaskedImage =
            lsst::fw::kernel::convolve(templateMaskedImage, *kernelPtr, threshold, edgeMaskBit);
        
        convolvedScienceMaskedImage -= convolvedTemplateMaskedImage;
        convolvedScienceMaskedImage.writeFits( (boost::format("%s_diff3") % inputImage).str() );
        
        // TEST : The output kernel is a gaussian with sigmaX = 1 and sigmaY = 2.
        //      : The kernel coefficients of all bases other than the first (mean) are 0.
        //      : The background function has no spatial variation and value equal to 0
        
    }
    
    if (Citizen::census(0) == 0) {
        cerr << "No leaks detected" << endl;
    } else {
        cerr << "Leaked memory blocks:" << endl;
        Citizen::census(cerr);
    } 
    
}
