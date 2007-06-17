// -*- lsst-c++ -*-
#include "lsst/fw/MaskedImage.h"
#include "lsst/fw/Trace.h"
#include "lsst/fw/Kernel.h"
#include "lsstimageproc03.h"
#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 

using namespace std;
using namespace lsst::fw;

// This piece of code includes the use case names as subroutines.
void getTemplateChunkExposureFromTemplateExposure() {
    wcsMatchExposure();
}
void wcsMatchExposure() {
}

void PSFMatchMaskedImages(MaskedImage convolveMaskedImage,
                          MaskedImage noncololveMaskedImage,
                          LinearCombinationKernel deltaFunctionKernelSet) {

    vector<Object> objectCollection;
    getCollectionOfMaskedImagesForPSFMatching(objectCollection);

    // Reusable view around each object
    MaskedImage<ImagePixelType,MaskPixelType>::MaskedImagePtrT convolvePtr;
    MaskedImage<ImagePixelType,MaskPixelType>::MaskedImagePtrT nonconvolvePtr;

    // Output kernels
    vector<boost::shared_ptr<Kernel<kernelPixelType> > > kernelVec;

    // Iterate over object; use iterator instead?
    for (unsigned nobj = 0; nobj < objectCollection.size(); nobj++) {
        Object diffImObject = objectCollection[nobj];
        
        // grab view around each object
        // do i really want a new stamp or just a view?
        BBox2i stamp(diffImObject.rowc - diffImObject.drow, diffImObject.rowc + diffImObject.drow, 
                     diffImObject.colc - diffImObject.dcol, diffImObject.colc + diffImObject.dcol);
        convolvePtr    = convolveMaskedImage.getSubImage(stamp);
        nonconvolvePtr = nonconvolveMaskedImage.getSubImage(stamp);

        boost::shared_ptr<lsst::fw::Kernel<pixelType> > kernelPtr(
            new Kernel<kernelPixelType>()
            );

        computePSFMatchingKernelForEachMaskedImage(convolvePtr, nonconvolvePtr, deltaFunctionKernelSet, kernelPtr);
        kernelVec.push_back(kernelPtr);
    }

    // Does nothing currently
    computeSpatiallyVaryingPSFMatchingKernel(kernelVec);
}
void getCollectionOfMaskedImagesForPSFMatching(vector<Object> objectCollection) {
    Object obj1;
    obj1.rowc = 100;
    obj1.colc = 100;
    obj1.drow = 10;
    obj1.dcol = 10;
    objectCollection.push_back(obj1);
}

void computePSFMatchingKernelForEachMaskedImage(MaskedImagePtrT convolvePtr,
                                                MaskedImagePtrT nonconvolvePtr,
                                                LinearCombinationKernel kernelSet,
                                                boost::shared_ptr kernelVec) {

    unsigned nKernelParameters=0, nBackgroundParameters=0, nParameters=0;
    for (unsigned i = 0; i < kernelSet.size(); i++)
        nKernelParameters += kernelSet[i].getNKernelParameters();
    
    // XXX; need to set up a Function for the background
    // Or, we just assume that across a single kernel, its 0th order.  This quite makes sense.
    nBackgroundParameters = 1;
    nParameters = nKernelParameters + nBackgroundParameters;
    
    vw::Vector<double> B(nParameters);
    vw::Matrix<double> M(nParameters, nParameters);

    // Calculate convolution of Reference image with Kernel
    // We can make this faster for delta function kernels
    // We might also want to pre-compute these things for each kernel, once we get to more complicated (Alard-Lupton) kernels.
    vector<boost::shared_ptr<MaskedImage<ImagePixelType, MaskPixelType> > > convolvedVec;
    for (unsigned ki1 = 0; ki1 < kernelSet.size(); ki1++) {
        for (unsigned ki2=0; ki2 < kernelSet[ki1].getNKernelParameters(); ki2++) {
            // The above loop is a placeholder for when we have multiple params to fit for per kernel
            // We need a kernel method that returns pixelized representations of all the meshes we need to convolve
            //   the image with at this stage.  For delta functions, or any single parameters functions, its just 
            //   itself so we are OK.
        }
        boost::shared_ptr<MaskedImage<ImagePixelType, MaskPixelType> >
            cImagePtr = Kernel::convolve(convolvePtr, kernelSet[ki1], 0, vw::NoEdgeExtension(), -1);
        convolvedVec.push_back(cImagePtr);
    }

    // integral over image's dx and dy
    for (unsigned row = 0; row < convolvePtr->getImage().nrow; row++) {
        for (unsigned col = 0; col < convolvePtr->getImage().ncol; col++) {
            
            // What is the correct way to access a pixel again?
            double ncCamera = nonconvolvePtr[row,col].camera();
            double ncVariance = nonconvolvePtr[row,col].variance();

            // Its an approximation to use this since you don't really know the kernel yet
            // You really want the post-convolved variance
            double cVariance = convolvePtr[row,col].variance();

            // Quicker computation?
            double iVariance = 1.0 / (ncVariance + cVariance);

            // kernel index i
            unsigned kidxi = 0;
            for (unsigned ki1 = 0; ki1 < kernelSet.size(); ki1++) {
                for (unsigned ki2 = 0; ki2 < kernelSet[ki1].getNKernelParameters(); ki2++, kidxi++) {
                    
                    double cdCamerai = convolvedVec[kidxi][row,col].camera();
                    B[kidxi] += ncCamera * cdCamerai * iVariance;
                    
                    // kernel index j 
                    unsigned kidxj = 0;
                    for (unsigned kj1 = ki1; kj1 < kernelSet.size(); kj1++) {
                        for (unsigned kj2 = 0; kj2 < kernelSet[kj2].getNKernelParameters(); kj2++, kidxj++) {
                            double cdCameraj = convolvedVec[kidxj][row,col].camera();
                            M[kidxi, kidxj] += cdCamerai * cdCameraj * iVariance;
                        } // kj2
                    } // kj1, kidxj

                } // ki2
            } // ki1, kidxi

            // Here we have a single background term
            B[kidxi+1] += ncCamera * iVariance;
            M[kidxi+1,kixdj+1] += 1.0 * iVariance;

        } // col
    } // row

    // Fill in rest of M
    for (unsigned kidxi=0; kidxi < nKernelParameters; kidxi++) 
        for (unsigned kidxj=kidxi+1; kidxj < nKernelParameters; kidxj++) 
            M[kidxj, kidxi] = M[kidxi, kidxj];
    
    // Invert M
    vw::Matrix<double> Minv = inverse(M);
    // Solve for x in Mx = B
    vw::Vector<double> Soln = B * Minv;

    // rearrange vector into kernel Image
    vw::MatrixProxy<double> Kimage(Ksize, Ksize, Soln);
    kernelVec = Kimage;  // somehow

    // clean up
    delete B;
    delete M;
    delete convolvedVec;
    delete Minv;
    delete Soln;
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
    
    typedef uint8 maskPixelType;
    typedef float32 imagePixelType;
    typedef float32 kernelPixelType;

    // Read input images
    string scienceInputImage = argv[1];
    MaskedImage<ImagePixelType,MaskPixelType> scienceMaskedImage;
    scienceMaskedImage.readFits(templateInputImage);

    string templateInputImage = argv[2];
    MaskedImage<ImagePixelType,MaskPixelType> templateMaskedImage;
    templateMaskedImage.readFits(templateInputImage);

    // set up basis of delta functions for kernel
    vector<boost::shared_ptr<Kernel<kernelPixelType> > > kernelVec;
    int kernelRows = 9;
    int kernelCols = 9;
    int colCtr = (kernelCols - 1) / 2;
    int rowCtr = (kernelRows - 1) / 2;
    for (unsigned row = 0; row < kernelRows; ++row) {
        int y = static_cast<int>(row) - rowCtr;
        
        for (unsigned col = 0; col < kernelCols; ++col) {
            int x = static_cast<int>(col) - colCtr;

            Kernel<kernelPixelType>::Function2PtrType kfuncPtr(
                new IntegerDeltaFunction2<kernelPixelType>(x, y)
            );
            boost::shared_ptr<Kernel<kernelPixelType> > kernelPtr(
                new AnalyticKernel<kernelPixelType>(kfuncPtr, kernelCols, kernelRows)
            );
            kernelVec.push_back(kernelPtr);

        }
    }
    vector<imagePixelType> kernelParams(kernelRows * kernelCols);
    LinearCombinationKernel<KernelPixelType> deltaFunctionKernelSet(kernelVec, kernelParams);
    
    // Currently does nothing
    getTemplateChunkExposureFromTemplateExposure();

    // This has some functionality!  Lets at least get it to compile.
    PSFMatchMaskedImages(scienceMaskedImage, templateMaskedImage, deltaFunctionKernelSet);

    // Currently does nothing
    subtractMatchedImage();
    
}
