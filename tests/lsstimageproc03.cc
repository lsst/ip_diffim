// -*- lsst-c++ -*-
#include "lsst/fw/Trace.h"
#include "lsst/fw/Kernel.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

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

    // Reusable view around each object
    MaskedImage<ImagePixelType,MaskPixelType>::MaskedImagePtrT convolvePtr;
    MaskedImage<ImagePixelType,MaskPixelType>::MaskedImagePtrT nonconvolvePtr;

    // Output kernels
    vector<boost::shared_ptr<Kernel<kernelPixelType> > > kernelVec;

    // Iterate over object; use iterator instead?
    for (unsigned nobj = 0; nobj < objectCollection.size(); nobj++) {
        Object<XXX> diffImObject = objectCollection[nobj];
        
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


    computeSpatiallyVaryingPSFMatchingKernel();
}
void getCollectionOfMaskedImagesForPSFMatching() {
}

void computePSFMatchingKernelForEachMaskedImage(MaskedImagePtrT convolvePtr,
                                                MaskedImagePtrT nonconvolvePtr,
                                                LinearCombinationKernel deltaFunctionKernelSet,
                                                boost::shared_ptr kernelVec) {
    boost::numeric::ublas::vector<double> B(XXX);
    boost::numeric::ublas::matrix<double> M(XXX, XXX);

    // integral over image's dx and dy
    for (unsigned row = 0; row < convolvePtr->getImage().nrow; row++) {
        for (unsigned col = 0; col < convolvePtr->getImage().ncol; col++) {

            double ncCamera   = nonconvolvePtr[row][col].camera();
            double ncVariance = nonconvolvePtr[row][col].variance();

            // kernel index i
            for (unsigned kerni = 0; kerni < deltaFunctionKernelSet.size(); kerni++) {
                kerneli = deltaFunctionKernelSet[kerni];

                cData   = convolvePtr[row][col];
                // This is effectively what we want to happen with the delta function kernel.
                // cData = convolvePtr[row + kerneli->drow][col + kerneli->dcol];
                // The more general case is
                // convolve
                boost::shared_ptr<MaskedImage<ImagePixelType, MaskPixelType> >
                    cImagePtri = Kernel::convolve(convolvePtr[row][col], kerneli, threshold, vw::NoEdgeExtension(), -1);
                
                B[kerneli] += ncData.camera() * cImagePtri.camera() / (ncData.variance() + cImagePtri.variance());
                
                // background term (is this cData right?)
                M[kerneli][nkernel] += cImagePtri.camera() / (cData.variance() + cImagePtri.variance());
                
                // kernel index j 
                for (unsigned kernj = 0; kernj < deltaFunctionSet.size(); kernj++) {
                    kernelj = deltaFunctionSet[kernj];
                    
                    // convolve
                    boost::shared_ptr<MaskedImage<ImagePixelType, MaskPixelType> >
                        cImagePtrj = Kernel::convolve(convolvePtr[row][col], kernelj, threshold, vw::NoEdgeExtension(), -1);
                    
                    M[kerni][kernj] += cImagePtri.camera() * cImagePtrj.camera() / (cImagePtri.variance() * cImagePtrj.variance());
                }
            }
            
            // background
            M[nkernel][nkernel] += 1. / (ncData.variance() + cData.variance());
            B[knernel]          += ncData.camera() / (ncData.variance() + cData.variance());
        } // col
    } // row
    
    vector LUv;
    matrix LUm = LUD(LUv, M);
    
    vector solution;
    LUSolve(LUm, B, LUv);


    // create a working copy of the input
    matrix<T> I(M);
 
    // create a permutation matrix for the LU-factorization
    boost::numeric::ublas::permutation_matrix<std::size_t> pm(I.size1());
    
    // perform LU-factorization
    int res = lu_factorize(I,pm);
    if( res != 0 ) return false;
    
    // create identity matrix of "inverse"
    inverse.assign(boost::numeric::ublas::identity_matrix<pixelType>(I.size1()));
    
    // backsubstitute to get the inverse
    lu_substitute(I, pm, inverse);



    // rearrange vector into kernel Image
    kernelVec = LUv;
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
    
    
    getTemplateChunkExposureFromTemplateExposure();
    PSFMatchMaskedImage();
    subtractMatchedImage();
    
}
