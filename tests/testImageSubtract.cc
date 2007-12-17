#include <fstream>

#include <boost/shared_ptr.hpp>
#include <boost/timer.hpp>

#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/FunctionLibrary.h>

#include <lsst/mwi/data/Citizen.h>
#include <lsst/mwi/exceptions/Exception.h>
#include <lsst/mwi/utils/Trace.h>
#include <lsst/mwi/policy/Policy.h>
#include <lsst/mwi/policy/paf/PAFParser.h>

#include <lsst/imageproc/ImageSubtract.h>

using namespace std;
using namespace lsst::fw;

typedef lsst::fw::maskPixelType MaskT;
typedef float ImageT;
typedef double KernelT;
typedef double FuncT;

int main( int argc, char** argv )
{
    {
        lsst::mwi::utils::Trace::setVerbosity("lsst.imageproc", 5);

        boost::timer t;
        double time;
        t.restart();

        // Read in Policy
        ifstream is("pipeline/ImageSubtractStageDictionary.paf");
        lsst::mwi::policy::Policy policy;
        lsst::mwi::policy::paf::PAFParser pp(policy);
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
        } catch (lsst::mwi::exceptions::ExceptionStack &e) {
            cerr << "Failed to open template image " << templateImage << ": " << e.what() << endl;
            return 1;
        }

        string scienceImage = argv[2];
        MaskedImage<ImageT,MaskT> scienceMaskedImage;
        try {
            scienceMaskedImage.readFits(scienceImage);
        } catch (lsst::mwi::exceptions::ExceptionStack &e) {
            cerr << "Failed to open science image " << scienceImage << ": " << e.what() << endl;
            return 1;
        }

        /* grab mask bits from the image to convolve, since that is what we'll be operating on */
        int badMaskBit = templateMaskedImage.getMask()->getMaskPlane("BAD");
        int edgeMaskBit = templateMaskedImage.getMask()->getMaskPlane("EDGE");
        MaskT badPixelMask = (badMaskBit < 0) ? 0 : (1 << badMaskBit);
        
        /* Debugging */
//        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> rowAcc2(templateMaskedImage);
//        for (unsigned int row = 0; row < templateMaskedImage.getRows(); ++row, rowAcc2.nextRow()) {
//            lsst::fw::MaskedPixelAccessor<ImageT, MaskT> colAcc2 = rowAcc2;
//            for (unsigned int col = 0; col < templateMaskedImage.getCols(); ++col, colAcc2.nextCol()) {
//                cout << row << " " << col << " " << *rowAcc2.image << " " << *rowAcc2.variance << " " << *rowAcc2.mask << endl;
//            }
//        }
       
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

        // Now you let the code find the peaks!
        vector<lsst::detection::Footprint::PtrType> footprintList =
            lsst::imageproc::getCollectionOfFootprintsForPsfMatching(
                templateMaskedImage, scienceMaskedImage, policy);

        // Do it
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > kernelPtr =
            lsst::imageproc::computePsfMatchingKernelForMaskedImage(
                kernelFunctionPtr, backgroundFunctionPtr, templateMaskedImage, scienceMaskedImage,
                kernelBasisVec, footprintList, policy);

        lsst::fw::MaskedImage<ImageT, MaskT> convolvedTemplateMaskedImage(templateMaskedImage.getCols(), 
                                                                          templateMaskedImage.getRows());
        lsst::fw::kernel::convolveLinear(convolvedTemplateMaskedImage, templateMaskedImage, *kernelPtr, edgeMaskBit);
        
        // Add background
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> rowAcc(convolvedTemplateMaskedImage);
        for (unsigned int row = 0; row < convolvedTemplateMaskedImage.getRows(); ++row, rowAcc.nextRow()) {
            lsst::fw::MaskedPixelAccessor<ImageT, MaskT> colAcc = rowAcc;
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
        lsst::imageproc::calculateMaskedImageResiduals(
            nGoodPixels, meanOfResiduals, varianceOfResiduals,
            scienceMaskedImage, badPixelMask);

        lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForMaskedImage", 4, 
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
