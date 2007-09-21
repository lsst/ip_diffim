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
        lsst::mwi::utils::Trace::setVerbosity(".", 4);
        
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
            cout << "  and uses it as both the template and the science image." << endl;
            cout << "  Your output kernel should be a delta function." << endl;
            cout << "  Basis function set is delta functions." << endl;
            cout << "  There is no spatial variation." << endl;
            cout << "  I add 100 counts to the background to test for this" << endl;
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
        
        // Generate basis of delta functions for kernel
        vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > kernelBasisVec;
        lsst::imageproc::generateDeltaFunctionKernelSet(kernelRows, kernelCols, kernelBasisVec);
        
        // Output kernel
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > kernelPtr(
            new lsst::fw::LinearCombinationKernel<KernelT>
            );
        
        // Function for spatially varying kernel.  Make null here for this test.
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > kernelFunctionPtr(
            new lsst::fw::function::PolynomialFunction2<FuncT>(kernelSpatialOrder)
            );
        
        // Function for spatially varying background.  
        boost::shared_ptr<lsst::fw::function::Function2<FuncT> > backgroundFunctionPtr(
            new lsst::fw::function::PolynomialFunction2<FuncT>(backgroundSpatialOrder)
            );

        scienceMaskedImage -= 100;

        // Use hard-coded positions for now
        vector<lsst::detection::Footprint::PtrType> footprintList;
        lsst::imageproc::getCollectionOfMaskedImagesForPsfMatching(footprintList);
        

        lsst::imageproc::computePsfMatchingKernelForMaskedImage
            (templateMaskedImage, scienceMaskedImage, kernelBasisVec, footprintList,
             kernelPtr, kernelFunctionPtr, backgroundFunctionPtr, p);
        
        lsst::fw::MaskedImage<ImageT, MaskT> convolvedTemplateMaskedImage(templateMaskedImage.getCols(), 
                                                                          templateMaskedImage.getRows());
        lsst::fw::kernel::convolveLinear(convolvedTemplateMaskedImage, templateMaskedImage, *kernelPtr, edgeMaskBit);
        convolvedTemplateMaskedImage.writeFits( (boost::format("%s_diff2XXX") % inputImage).str() );

        lsst::fw::MaskedImage<ImageT, MaskT> 
            convolvedTemplateMaskedImage2 = lsst::fw::kernel::convolve(templateMaskedImage, *kernelPtr, convolveThreshold, edgeMaskBit);
        convolvedTemplateMaskedImage2.writeFits( (boost::format("%s_diff2YYY") % inputImage).str() );

        // Subtract off template
        scienceMaskedImage -= convolvedTemplateMaskedImage;

        scienceMaskedImage.writeFits( (boost::format("%s_diff2ZZZ") % inputImage).str() );

        // Subtract off background
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> accessorCol(scienceMaskedImage);
        for (unsigned int col = 0; col < scienceMaskedImage.getCols(); ++col) {
            lsst::fw::MaskedPixelAccessor<ImageT, MaskT> accessorRow = accessorCol;
            for (unsigned int row = 0; row < scienceMaskedImage.getRows(); ++row) {
                *accessorRow.image += (*backgroundFunctionPtr)(col, row);
                accessorRow.nextRow();
            }
            accessorCol.nextCol();
        }

        // NOTE - Might be slicker to make an image from the background model.
        //        Then add in science image, and subtract off convolved image

        scienceMaskedImage.writeFits( (boost::format("%s_diff2") % inputImage).str() );
    }
    
    if (Citizen::census(0) == 0) {
        cerr << "No leaks detected" << endl;
    } else {
        cerr << "Leaked memory blocks:" << endl;
        Citizen::census(cerr);
    } 
    
}
