#include "lsst/fw/MaskedImage.h"
#include "lsst/fw/Kernel.h"

// Since we don't have real objects
struct Object {
    double rowc, colc;
    double drow, dcol;
};

void getTemplateChunkExposureFromTemplateExposure();
void wcsMatchExposure();
void PSFMatchMaskedImages(MaskedImage convolveMaskedImage, MaskedImage noncololveMaskedImage, LinearCombinationKernel deltaFunctionKernelSet);
void getCollectionOfMaskedImagesForPSFMatching(vector<Object> objectCollection);

void computePSFMatchingKernelForEachMaskedImage(MaskedImagePtrT convolvePtr,
                                                MaskedImagePtrT nonconvolvePtr,
                                                LinearCombinationKernel kernelSet,
                                                boost::shared_ptr kernelVec);
void computeSpatiallyVaryingPSFMatchingKernel();
void fitKernelsUsingPrincipalComponentAnalysis();
void fitArraysUsingPrincipalComponentAnalysis();

