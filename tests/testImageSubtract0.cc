#include <lsst/afw/MaskedImage.h>
#include <lsst/afw/Kernel.h>
#include <lsst/afw/KernelFunctions.h>
#include <lsst/afw/FunctionLibrary.h>
#include <lsst/daf/base/Citizen.h>
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace lsst::afw;

// This is for debugging MaskedImage; has no testing other than it has to succeed

int main( int argc, char** argv )
{
    // Cool, for Citizen!
    { 
        lsst::pex::logging::Trace::setVerbosity("lsst.ip.diffim", 4);
        
        typedef lsst::afw::maskPixelType MaskT;
        typedef double ImageT; 
        typedef double KernelT;
        const KernelT CONVOLVE_THRESHOLD = 0;
        const int EDGE_MASK_BIT = 1;
        
        string inputImage = argv[1];
        MaskedImage<ImageT,MaskT> scienceMaskedImage;
        scienceMaskedImage.readFits(inputImage);
        
        Kernel<KernelT>::KernelFunctionPtrType kfuncPtr(
            new lsst::afw::function::IntegerDeltaFunction2<KernelT>(0, 0)
            );
        
        boost::shared_ptr<Kernel<KernelT> > kernelPtr(
            new AnalyticKernel<KernelT>(kfuncPtr, 3, 3)
            );
        
        
        BBox2i stamp(100, 100, 10, 10);
        MaskedImage<ImageT,MaskT>::MaskedImagePtrT stampPtr = scienceMaskedImage.getSubImage(stamp);
        
        // Works
        //MaskedImage<ImageT, MaskT>
        //convolvedStamp = (*stampPtr);
        
        // Now Works
        MaskedImage<ImageT, MaskT>
            convolvedStamp = lsst::afw::kernel::convolve(*stampPtr, *kernelPtr, CONVOLVE_THRESHOLD, EDGE_MASK_BIT);        
        
        // Write out the stamp
        stampPtr->writeFits( "test0_stamp" );
        
        // Write out the kernel
        double imSum;
        Image<KernelT> kImage = kernelPtr->computeNewImage(imSum, 0.0, 0.0, false);
        kImage.writeFits( "test0_kernel.fits" );
        
        // Write out the convolved image
        convolvedStamp.writeFits( "test0_convolved" );
        
        // Subtract the stamp from the convolved image
        convolvedStamp -= (*stampPtr);
        
        // Write out the convolved image
        convolvedStamp.writeFits( "test0_subtracted" );
    }
    
    if (Citizen::census(0) == 0) {
        cerr << "No leaks detected" << endl;
    } else {
        cerr << "Leaked memory blocks:" << endl;
        Citizen::census(cerr);
    } 

    return 1;
}
