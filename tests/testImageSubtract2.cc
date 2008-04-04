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
//        KernelT convolveThreshold = static_cast<KernelT>(p.getDouble("convolveThreshold"));
        int badMaskBit = p.getInt("badMaskBit");
        int edgeMaskBit = p.getInt("edgeMaskBit");
        unsigned int kernelCols = p.getInt("kernelCols");
        unsigned int kernelRows = p.getInt("kernelRows");
        unsigned int kernelSpatialOrder = p.getInt("kernelSpatialOrder");
        unsigned int backgroundSpatialOrder = p.getInt("backgroundSpatialOrder");

        MaskT edgePixelMask = (edgeMaskBit < 0) ? 0 : (1 << edgeMaskBit);
        MaskT badPixelMask = (badMaskBit < 0) ? 0 : (1 << badMaskBit);
        badPixelMask |= edgePixelMask;

        // Read input images
        if (argc < 3) {
            cout << "This program takes a single input image on the command line (and DC offset)" << endl;
            cout << "  and uses it as both the template and the science image." << endl;
            cout << "  Your output kernel should be a delta function." << endl;
            cout << "  Basis function set is delta functions." << endl;
            cout << "  There is no spatial variation." << endl;
            cout << "  I add DC offset counts to the background to test for this" << endl;
            exit(1);
        }
        string inputImage = argv[1];
        MaskedImage<ImageT,MaskT> scienceMaskedImage;
        scienceMaskedImage.readFits(inputImage);

        float offset = atof(argv[2]);
        
        MaskedImage<ImageT,MaskT> templateMaskedImage;
        templateMaskedImage.readFits(inputImage);
        
        // Generate basis of delta functions for kernel
        vector<boost::shared_ptr<lsst::afw::Kernel<KernelT> > > kernelBasisVec =
            lsst::ip::diffim::generateDeltaFunctionKernelSet<KernelT>(kernelCols, kernelRows);
        
        // Function for spatially varying kernel.  Make null here for this test.
        boost::shared_ptr<lsst::afw::function::Function2<FuncT> > kernelFunctionPtr(
            new lsst::afw::function::PolynomialFunction2<FuncT>(kernelSpatialOrder)
            );
        
        // Function for spatially varying background.  
        boost::shared_ptr<lsst::afw::function::Function2<FuncT> > backgroundFunctionPtr(
            new lsst::afw::function::PolynomialFunction2<FuncT>(backgroundSpatialOrder)
            );

        scienceMaskedImage -= offset;

        // Use hard-coded positions for now
        vector<lsst::detection::Footprint::PtrType> footprintList =
            lsst::ip::diffim::getCollectionOfMaskedImagesForPsfMatching();
        
        boost::shared_ptr<lsst::afw::LinearCombinationKernel<KernelT> > kernelPtr =
            lsst::ip::diffim::computePsfMatchingKernelForMaskedImage(
                kernelFunctionPtr, backgroundFunctionPtr, templateMaskedImage, scienceMaskedImage,
                kernelBasisVec, footprintList, p);
        
        lsst::afw::MaskedImage<ImageT, MaskT> convolvedTemplateMaskedImage(templateMaskedImage.getCols(), 
                                                                          templateMaskedImage.getRows());
        lsst::afw::kernel::convolveLinear(convolvedTemplateMaskedImage, templateMaskedImage, *kernelPtr, edgeMaskBit);

        // Add background
        lsst::afw::MaskedPixelAccessor<ImageT, MaskT> rowAcc(convolvedTemplateMaskedImage);
        for (unsigned int row = 0; row < convolvedTemplateMaskedImage.getRows(); ++row, rowAcc.nextRow()) {
            lsst::afw::MaskedPixelAccessor<ImageT, MaskT> colAcc = rowAcc;
            for (unsigned int col = 0; col < convolvedTemplateMaskedImage.getCols(); ++col, colAcc.nextCol()) {
                *colAcc.image += (*backgroundFunctionPtr)(col, row);
            }
        }

        // Subtract off template
        scienceMaskedImage -= convolvedTemplateMaskedImage;
        scienceMaskedImage.writeFits( (boost::format("%s_diff2") % inputImage).str() );
    }
    
    if (Citizen::census(0) == 0) {
        cerr << "No leaks detected" << endl;
    } else {
        cerr << "Leaked memory blocks:" << endl;
        Citizen::census(cerr);
    } 
    
}
