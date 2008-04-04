#include <fstream>

#include <boost/shared_ptr.hpp>

#include <lsst/afw/MaskedImage.h>
#include <lsst/afw/Kernel.h>
#include <lsst/afw/FunctionLibrary.h>

#include <lsst/daf/base/Citizen.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/policy/Policy.h>
#include <lsst/pex/policy/paf/PAFParser.h>

#include <lsst/ip/diffim/ImageSubtract.h>

#define DEBUG_IO 1

using namespace std;
using namespace lsst::afw;

int main( int argc, char** argv )
{
    {
        lsst::pex::logging::Trace::setVerbosity("lsst.ip.diffim", 4);
        
        typedef lsst::afw::maskPixelType MaskT;
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
            cout << "  The science image is derived from the template convolved with a non-spatially varying Gaussian." << endl;
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
        
        lsst::afw::Kernel<KernelT>::KernelFunctionPtrType gaussFuncPtr(
            new lsst::afw::function::GaussianFunction2<FuncT>(sigmaX, sigmaY));
        lsst::afw::AnalyticKernel<KernelT> gaussKernel(gaussFuncPtr, kernelCols, kernelRows);

        // write out kernel
        double imSum;
        lsst::afw::Image<KernelT> kImage = gaussKernel.computeNewImage(imSum, 0.0, 0.0, false);
        kImage.writeFits( (boost::format("k3Fits.fits")).str() );
        
        // Convolved science image
        lsst::pex::logging::Trace("testImageSubtract3", 2, "Convolving input image for testing");
        lsst::afw::MaskedImage<ImageT, MaskT> convolvedScienceMaskedImage =
            lsst::afw::kernel::convolve(scienceMaskedImage, gaussKernel, convolveThreshold, edgeMaskBit, false);

        // Generate basis of delta functions for kernel
        vector<boost::shared_ptr<Kernel<KernelT> > > kernelBasisVec =
            lsst::ip::diffim::generateDeltaFunctionKernelSet<KernelT>(kernelCols, kernelRows);
        
        // Function for spatially varying kernel.  Make null here for this test.
        boost::shared_ptr<lsst::afw::function::Function2<FuncT> > kernelFunctionPtr(
            new lsst::afw::function::PolynomialFunction2<FuncT>(kernelSpatialOrder)
            );
        
        // Function for spatially varying background.  
        boost::shared_ptr<lsst::afw::function::Function2<FuncT> > backgroundFunctionPtr(
            new lsst::afw::function::PolynomialFunction2<FuncT>(backgroundSpatialOrder)
            );
        
        // Use hard-coded positions for now
        vector<lsst::detection::Footprint::PtrType> footprintList =
            lsst::ip::diffim::getCollectionOfMaskedImagesForPsfMatching();

        boost::shared_ptr<lsst::afw::LinearCombinationKernel<KernelT> > kernelPtr =        
            lsst::ip::diffim::computePsfMatchingKernelForMaskedImage(
            kernelFunctionPtr, backgroundFunctionPtr, templateMaskedImage, convolvedScienceMaskedImage,
            kernelBasisVec, footprintList, p);
        
        lsst::afw::MaskedImage<ImageT, MaskT> convolvedTemplateMaskedImage(templateMaskedImage.getCols(), 
                                                                          templateMaskedImage.getRows());
        lsst::afw::kernel::convolveLinear(convolvedTemplateMaskedImage, templateMaskedImage, *kernelPtr, edgeMaskBit);

        // Add background
        lsst::afw::MaskedPixelAccessor<ImageT, MaskT> accessorCol(convolvedTemplateMaskedImage);
        for (unsigned int col = 0; col < convolvedTemplateMaskedImage.getCols(); ++col) {
            lsst::afw::MaskedPixelAccessor<ImageT, MaskT> accessorRow = accessorCol;
            for (unsigned int row = 0; row < convolvedTemplateMaskedImage.getRows(); ++row) {
                *accessorRow.image += (*backgroundFunctionPtr)(col, row);
                accessorRow.nextRow();
            }
            accessorCol.nextCol();
        }

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
