#include <fstream>

#include <boost/shared_ptr.hpp>
#include <boost/timer.hpp>

#include <lsst/afw/MaskedImage.h>
#include <lsst/afw/Kernel.h>
#include <lsst/afw/FunctionLibrary.h>

#include <lsst/daf/base/Citizen.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/policy/Policy.h>
#include <lsst/pex/policy/paf/PAFParser.h>

#include <lsst/ip/diffim/ImageSubtract.h>

using namespace std;
using namespace lsst::afw;

typedef lsst::afw::maskPixelType MaskT;
typedef float ImageT;
typedef double KernelT;
typedef double FuncT;

int main( int argc, char** argv )
{
    {
        lsst::pex::logging::Trace::setVerbosity("lsst.ip.diffim", 5);

        boost::timer t;
        double time;
        t.restart();

        // Read in Policy
        ifstream is("pipeline/ImageSubtractStageDictionary.paf");
        lsst::pex::policy::Policy policy;
        lsst::pex::policy::paf::PAFParser pp(policy);
        pp.parse(is);
        is.close();

        // Parse policy
        unsigned int kernelCols = policy.getInt("kernelCols");
        unsigned int kernelRows = policy.getInt("kernelRows");
        unsigned int kernelSpatialOrder = policy.getInt("kernelSpatialOrder");
        unsigned int backgroundSpatialOrder = policy.getInt("backgroundSpatialOrder");

       
        // Read input images
        if (argc < 2) {
            cout << "This program takes 2 input images on the command line" << endl;
            cout << "  and finds its own detection." << endl;
            cout << "  Basis function set is delta functions." << endl;
            cout << "  There is spatial variation." << endl;
            cout << " E.g. " << endl;
            cout << " ./tests/testImageSubtract5 $FWDATA_DIR/SM/sme9_10_tmpl $FWDATA_DIR/SM/sme9_10_im " << endl;
            exit(1);
        }
        string templateImage = argv[1];
        MaskedImage<ImageT,MaskT> templateMaskedImage;
        try {
            templateMaskedImage.readFits(templateImage);
        } catch (lsst::pex::exceptions::ExceptionStack &e) {
            cerr << "Failed to open template image " << templateImage << ": " << e.what() << endl;
            return 1;
        }

        string scienceImage = argv[2];
        MaskedImage<ImageT,MaskT> scienceMaskedImage;
        try {
            scienceMaskedImage.readFits(scienceImage);
        } catch (lsst::pex::exceptions::ExceptionStack &e) {
            cerr << "Failed to open science image " << scienceImage << ": " << e.what() << endl;
            return 1;
        }

        /* grab mask bits from the image to convolve, since that is what we'll be operating on */
        int badMaskBit = templateMaskedImage.getMask()->getMaskPlane("BAD");
        int edgeMaskBit = templateMaskedImage.getMask()->getMaskPlane("EDGE");
        MaskT badPixelMask = (badMaskBit < 0) ? 0 : (1 << badMaskBit);
        
        /* Debugging */
//        lsst::afw::MaskedPixelAccessor<ImageT, MaskT> rowAcc2(templateMaskedImage);
//        for (unsigned int row = 0; row < templateMaskedImage.getRows(); ++row, rowAcc2.nextRow()) {
//            lsst::afw::MaskedPixelAccessor<ImageT, MaskT> colAcc2 = rowAcc2;
//            for (unsigned int col = 0; col < templateMaskedImage.getCols(); ++col, colAcc2.nextCol()) {
//                cout << row << " " << col << " " << *rowAcc2.image << " " << *rowAcc2.variance << " " << *rowAcc2.mask << endl;
//            }
//        }
       
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

        // Now you let the code find the peaks!
        vector<lsst::detection::Footprint::PtrType> footprintList =
            lsst::ip::diffim::getCollectionOfFootprintsForPsfMatching(
                templateMaskedImage, scienceMaskedImage, policy);

        // Do it
        boost::shared_ptr<lsst::afw::LinearCombinationKernel<KernelT> > kernelPtr =
            lsst::ip::diffim::computePsfMatchingKernelForMaskedImage(
                kernelFunctionPtr, backgroundFunctionPtr, templateMaskedImage, scienceMaskedImage,
                kernelBasisVec, footprintList, policy);

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
        scienceMaskedImage.writeFits( (boost::format("%s_diff") % scienceImage).str() );

        // Find quality metrics.
        double meanOfResiduals = 0.0;
        double varianceOfResiduals = 0.0;
        int nGoodPixels = 0;
        lsst::ip::diffim::calculateMaskedImageResiduals(
            nGoodPixels, meanOfResiduals, varianceOfResiduals,
            scienceMaskedImage, badPixelMask);

        lsst::pex::logging::Trace("lsst.ip.diffim.computePsfMatchingKernelForMaskedImage", 4, 
                                (boost::format("Mean and variance of residuals in difference image : %.3f %.3f (%d pixels)") 
                                 % meanOfResiduals % varianceOfResiduals % nGoodPixels));

        time = t.elapsed();
        cout << "Difference imaging took " << time << "s" << endl;

    }
    
    if (Citizen::census(0) == 0) {
        cerr << "No leaks detected" << endl;
    } else {
        cerr << "Leaked memory blocks:" << endl;
        Citizen::census(cerr);
    } 
}
