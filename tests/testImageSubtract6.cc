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
        
        lsst::fw::Kernel<KernelT>::KernelFunctionPtrType gaussFuncPtr(
            new lsst::fw::function::GaussianFunction2<FuncT>(sigmaX, sigmaY));
        lsst::fw::AnalyticKernel<KernelT> gaussKernel(gaussFuncPtr, kernelCols, kernelRows);
        
        // Convolved science image
        lsst::mwi::utils::Trace("testImageSubtract6", 2, "Convolving input image for testing");
        lsst::fw::MaskedImage<ImageT, MaskT> convolvedScienceMaskedImage =
            lsst::fw::kernel::convolve(scienceMaskedImage, gaussKernel, convolveThreshold, edgeMaskBit, false);
        
        scienceMaskedImage *= 1.7;
        scienceMaskedImage -= 75;
        
        // Generate basis of delta functions for kernel
        vector<boost::shared_ptr<Kernel<KernelT> > > kernelBasisVec =
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
                kernelFunctionPtr, backgroundFunctionPtr, templateMaskedImage, convolvedScienceMaskedImage,
                kernelBasisVec, footprintList, p);
        
        lsst::fw::MaskedImage<ImageT, MaskT> convolvedTemplateMaskedImage(templateMaskedImage.getCols(), 
                                                                          templateMaskedImage.getRows());
        lsst::fw::kernel::convolveLinear(convolvedTemplateMaskedImage, templateMaskedImage, *kernelPtr, edgeMaskBit);

        // Add background
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> accessorCol(convolvedTemplateMaskedImage);
        for (unsigned int col = 0; col < convolvedTemplateMaskedImage.getCols(); ++col) {
            lsst::fw::MaskedPixelAccessor<ImageT, MaskT> accessorRow = accessorCol;
            for (unsigned int row = 0; row < convolvedTemplateMaskedImage.getRows(); ++row) {
                *accessorRow.image += (*backgroundFunctionPtr)(col, row);
                accessorRow.nextRow();
            }
            accessorCol.nextCol();
        }
        
        convolvedScienceMaskedImage -= convolvedTemplateMaskedImage;
        convolvedScienceMaskedImage.writeFits( (boost::format("%s_diff6") % inputImage).str() );
        
    }
    
    if (Citizen::census(0) == 0) {
        cerr << "No leaks detected" << endl;
    } else {
        cerr << "Leaked memory blocks:" << endl;
        Citizen::census(cerr);
    } 
    
}
