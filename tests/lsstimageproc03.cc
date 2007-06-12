// -*- lsst-c++ -*-
#include "lsst/fw/Trace.h"
#include "lsst/fw/Kernel.h"

using namespace std;
using namespace lsst::fw;

// This piece of code includes the use case names as subroutines.
void getTemplateChunkExposureFromTemplateExposure() {
   wcsMatchExposure();
}
void wcsMatchExposure() {
}

void PSFMatchMaskedImage() {
   getCollectionOfMaskedImagesForPSFMatching();
   computePSFMatchingKernelForEachMaskedImage();
   computeSpatiallyVaryingPSFMatchingKernel();
}
void getCollectionOfMaskedImagesForPSFMatching() {
}
void computePSFMatchingKernelForEachMaskedImage() {
}
void computeSpatiallyVaryingPSFMatchingKernel() {
   fitKernelsUsingPrincipalComponentAnalysis();
}
void fitKernelsUsingPrincipalComponentAnalysis() {
   fitArraysUsingPrincipalComponentAnalysis();
}
void fitArraysUsingPrincipalComponentAnalysis() {
}


int main( int argc, char** argv )
{
   // This is effectively Subtract Template From Science Chunk Exposure

    Trace::setDestination(cout);
    Trace::setVerbosity(".", 0);
    
    typedef uint8 MaskPixelType;
    typedef float32 ImagePixelType;

    getTemplateChunkExposureFromTemplateExposure();
    PSFMatchMaskedImage();
    subtractMatchedImage();

}
