#include <fstream>

#include <boost/shared_ptr.hpp>

#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/math/FunctionLibrary.h>

#include <lsst/daf/base/Citizen.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/loggin/Trace.h>
#include <lsst/pex/policy/Policy.h>
#include <lsst/pex/policy/paf/PAFParser.h>

#include <lsst/ip/diffim/ImageSubtract.h>

#define DEBUG_IO 1

using namespace std;
using namespace lsst::afw;

int main( int argc, char** argv )
{
    {
        lsst::pex::loggin::Trace::setVerbosity("lsst.ip.diffim", 4);
        
        typedef lsst::afw::image::maskPixelType MaskT;
        typedef float ImageT; // have to make sure this jibes with the input data!
        typedef double KernelT;
        typedef double FuncT;

        // Read in Policy
        ifstream is("tests/ImageSubtract_policy.paf");
        lsst::pex::policy::Policy p;
        lsst::pex::policy::paf::PAFParser pp(p);
        pp.parse(is);
        is.close();

        // Parse policy
        KernelT convolveThreshold = static_cast<KernelT>(p.getDouble("convolveThreshold"));
//        int badMaskBit = p.getInt("badMaskBit");
        int edgeMaskBit = p.getInt("edgeMaskBit");
        unsigned int kernelCols = p.getInt("kernelCols");
        unsigned int kernelRows = p.getInt("kernelRows");
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
        scienceMaskedImage.readFits(inputImage);
        
        MaskedImage<ImageT,MaskT> templateMaskedImage;
        templateMaskedImage.readFits(inputImage);
        
        // The kernel to convolve the template image with to yield the science image
        double sigmaX = 1.0;
        double sigmaY = 2.0;
        unsigned int polyOrder = 2;
        double minSigma = 0.1;
        double maxSigma = 3.0;
        
        lsst::afw::math::Kernel<KernelT>::KernelFunctionPtrType gaussFuncPtr(
            new lsst::afw::math::GaussianFunction2<FuncT>(sigmaX, sigmaY));
        lsst::afw::math::Kernel<KernelT>::SpatialFunctionPtrType polyFuncPtr(
            new lsst::afw::math::PolynomialFunction2<double>(polyOrder));
        lsst::afw::math::AnalyticKernel<KernelT> gaussSpVarKernel(
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
        lsst::pex::loggin::Trace("testImageSubtract4", 2, "Convolving input image for testing");
        lsst::afw::image::MaskedImage<ImageT, MaskT> convolvedScienceMaskedImage =
            lsst::afw::math::convolve(scienceMaskedImage, gaussSpVarKernel, convolveThreshold, edgeMaskBit, false);
        
        convolvedScienceMaskedImage.writeFits( (boost::format("%s_test4") % inputImage).str() );
        
        // Generate basis of delta functions for kernel
        vector<boost::shared_ptr<Kernel<KernelT> > > kernelBasisVec =
            lsst::ip::diffim::generateDeltaFunctionKernelSet<KernelT>(kernelCols, kernelRows);
                
        // Function for spatially varying kernel. 
        // Give it less power than the Fcn we convolved the image with, see what happens
        // OVERRIDE POLICY HERE
        kernelSpatialOrder = polyOrder;
        boost::shared_ptr<lsst::afw::math::Function2<FuncT> > kernelFunctionPtr(
            new lsst::afw::math::PolynomialFunction2<FuncT>(kernelSpatialOrder)
            );
        
        // Function for spatially varying background.  
        boost::shared_ptr<lsst::afw::math::Function2<FuncT> > backgroundFunctionPtr(
            new lsst::afw::math::PolynomialFunction2<FuncT>(backgroundSpatialOrder)
            );

        // Use hard-coded positions for now
        vector<lsst::detection::Footprint::PtrType> footprintList =
            lsst::ip::diffim::getCollectionOfMaskedImagesForPsfMatching();

        boost::shared_ptr<lsst::afw::math::LinearCombinationKernel<KernelT> > kernelPtr =
            lsst::ip::diffim::computePsfMatchingKernelForMaskedImage(
                kernelFunctionPtr, backgroundFunctionPtr, templateMaskedImage, convolvedScienceMaskedImage,
                kernelBasisVec, footprintList, p);
        
        lsst::afw::image::MaskedImage<ImageT, MaskT> convolvedTemplateMaskedImage(templateMaskedImage.getCols(), 
                                                                          templateMaskedImage.getRows());
        lsst::afw::math::convolveLinear(convolvedTemplateMaskedImage, templateMaskedImage, *kernelPtr, edgeMaskBit);
        
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
