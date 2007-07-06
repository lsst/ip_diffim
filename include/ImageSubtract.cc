// -*- lsst-c++ -*-
/**
 * \file
 *
 * Implementation of image subtraction
 *
 * Implementation of image subtraction
 *
 * \author Andrew Becker
 *
 * \ingroup imageproc
 */

#include "lsst/fw/MaskedImage.h"
#include "lsst/fw/Trace.h"
#include "lsst/fw/Kernel.h"
#include "lsstimageproc03.h"
#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 
#include "ImageSubtract.h"
using namespace std;
using namespace lsst::fw;

// Inline Member Functions
inline unsigned
lsst::fw::Obejct::getColc() const {
    return _colc;
}
inline unsigned
lsst::fw::Obejct::getRowc() const {
    return _rowc;
}
inline unsigned
lsst::fw::Obejct::getDcol() const {
    return _dcol;
}
inline unsigned
lsst::fw::Obejct::getDrow() const {
    return _drow;
}

/**
 * Computes spatially varying PSF matching kernel for image subtraction
 *
 * Note: Longer description here
 *
 * \return This describes the return value if there is one
 * \throw Any exceptions thrown must be described here
 * \throw Here too
 * \ingroup imageproc
 */
template <PixelT, MaskT, KernelT>
static void computePSFMatchingKernelForMaskedImage(
    MaskedImage<PixelT,MaskT> const &imageToConvolve, ///< Template image; convolved
    MaskedImage<PixelT,MaskT> const &imageToNotConvolve, ///< Science image; not convolved
    LinearCombinationKernel<KernelT> &kernelBasisSet, ///< Input set of basis kernels; modified on return
    ) {

    vector<Object> objectCollection;
    getCollectionOfMaskedImagesForPSFMatching(&objectCollection);

    // Reusable view around each object
    MaskedImage<PixelT,MaskT>::MaskedImagePtrT imageToConvolvePtr;
    MaskedImage<PixelT,MaskT>::MaskedImagePtrT imageToNotConvolvePtr;

    // Output kernels
    vector<boost::shared_ptr<Kernel<kernelPixelType> > > kernelVec;

    // Iterate over object; use iterator instead?
    for (unsigned nobj = 0; nobj < objectCollection.size(); nobj++) {
        Object diffImObject = objectCollection[nobj];
        
        // grab view around each object
        // do i really want a new stamp or just a view?
        BBox2i stamp(diffImObject.rowc - diffImObject.drow, diffImObject.rowc + diffImObject.drow, 
                     diffImObject.colc - diffImObject.dcol, diffImObject.colc + diffImObject.dcol);
        imageToConvolvePtr    = imageToConvolve.getSubImage(stamp);
        imageToNotConvolvePtr = imageToNotConvolve.getSubImage(stamp);

        std::vector<KernelT> kernelCoeffs(kernelBasisSet.getNKernelParameters());

        computePSFMatchingKernelForPostageStamp(imageToConvolvePtr, imageToNotConvolvePtr, kernelBasisSet, kernelCoeffs);
        kernelCoeffsVec.push_back(kernelCoeffs);
    }

    // Does nothing currently
    computeSpatiallyVaryingPSFMatchingKernel(kernelBasisSet, kernelCoeffsVec);
}

/**
 * Single line description with no period
 *
 * Note: Long description
 *
 * \return This describes the return value if there is one
 * \throw Any exceptions thrown must be described here
 * \throw Here too
 * \ingroup I guess imageproc for me
 */
template <PixelT, MaskT, KernelT>
static void computePSFMatchingKernelForPostageStamp(
    MaskedImage<PixelT, MaskT> const &imageToConvolve, ///< Goes with the code
    MaskedImage<PixelT, MaskT> const &imageToNotConvolve, ///< This is for doxygen
    LinearCombinationKernel<KernelT> &kernelBasisSet, ///< This is for doxygen
    std::vector<KernelT> &kernelCoeffs) { ///< This is for doxygen

    unsigned nKernelParameters=0, nBackgroundParameters=0, nParameters=0;
    
    // We assume that each kernel in the Set has 1 parameter you fit for
    nKernelParameters = kernelBasisSet.getNKernelParameters();
    // Or, we just assume that across a single kernel, background 0th order.  This quite makes sense.
    nBackgroundParameters = 1;
    // Total number of parameters
    nParameters = nKernelParameters + nBackgroundParameters;
    
    vw::Vector<double> B(nParameters);
    vw::Matrix<double> M(nParameters, nParameters);

    // Calculate convolution of Reference image with Kernel
    // We can make this faster for delta function kernels
    vector<MaskedImage<PixelT, MaskT> > convolvedImageVec(nKernelParameters);
    for (unsigned ki = 0; ki < nKernelParameters; ki++) {
        MaskedImage<PixelT, MaskT>
            convolvedImage = Kernel::convolve(imageToConvolve, kernelBasisSet[ki], 0, vw::NoEdgeExtension());
        convolvedImageVec.push_back(convolvedImage);
    }

    // An accessor for each convolution plane
    vector<MaskedPixelAccessor<PixelT, MaskT> > convolvedAccessorRowVec(nKernelParameters);
    for (unsigned ki = 0; ki < nKernelParameters; ki++) {
        MaskedPixelAccessor<PixelT, MaskT> convolvedAccessorRow(convolvedImageVec[i]);
        convolvedAccessorRowVec.push_back(convolvedAccessorRow);
    }

    // An accessor for each input image
    MaskedPixelAccessor<PixelT, MaskT> imageToConvolve(imageToConvolve);
    MaskedPixelAccessor<PixelT, MaskT> imageToNotConvolveRow(imageToNotConvolve);

    // integral over image's dx and dy
    for (unsigned row = 0; row < imageToConvolve.getRows(); row++) {

        // An accessor for each convolution plane
        vector<MaskedPixelAccessor<PixelT, MaskT> > convolvedAccessorColVec(nKernelParameters);
        for (unsigned ki = 0; ki < nKernelParameters; ki++) {
            MaskedPixelAccessor<PixelT, MaskT> convolvedAccessorCol = convolvedAccessorRow;
            convolvedAccessorColVec.push_back(convolvedAccessorCol);
        }

        // An accessor for each input image
        MaskedPixelAccessor<PixelT, MaskT> imageToConvolveCol = imageToConvolveRow;
        MaskedPixelAccessor<PixelT, MaskT> imageToNotConvolveCol = imageToNotConvolveRow;

        for (unsigned col = 0; col < imageToConvolve.getCols(); col++) {
            
            PixelT ncCamera = *nonconvolveCol.image;
            PixelT ncVariance = *nonconvolveCol.variance;

            // Its an approximation to use this since you don't really know the kernel yet
            // You really want the post-convolved variance but this is close enough
            PixelT cVariance = *convolveCol.variance;
            
            // Quicker computation
            double iVariance = 1.0 / (ncVariance + cVariance);

            // kernel index i
            for (unsigned kidxi = 0; kidxi < nKernelParameters; kidxi++) {
                PixelT cdCamerai = *convolvedAccessorColVec[kidxi].image;
                B[kidxi] += ncCamera * cdCamerai * iVariance;
                
                // kernel index j 
                for (unsigned kidxj = kidxi; kidxj < nKernelParameters; kidxj++) {
                    PixelT cdCameraj = *convolvedAccessorColVec[kidxj].image;
                    M[kidxi, kidxj] += cdCamerai * cdCameraj * iVariance;
                } // kidxj
            } // kidxi

            // Here we have a single background term
            B[kidxi+1] += ncCamera * iVariance;
            M[kidxi+1,kixdj+1] += 1.0 * iVariance;

            // Step each accessor in column
            imageToConvolveCol.nextCol();
            imageToNotConvolveCol.nextCol();
            for (unsigned ki = 0; ki < nKernelParameters; ki++) {
                convolvedAccessorColVec[i].nextCol();
            }
        } // col
        
        // Step each accessor in row
        imageToConvolveRow.nextRow();
        imageToNotConvolveRow.nextRow();
        for (unsigned ki = 0; ki < nKernelParameters; ki++) {
            convolvedAccessorRowVec[i].nextRow();
        }

        // clean up
        delete convolvedAccessorColVec;
    } // row

    // Fill in rest of M
    for (unsigned kidxi=0; kidxi < nKernelParameters; kidxi++) 
        for (unsigned kidxj=kidxi+1; kidxj < nKernelParameters; kidxj++) 
            M[kidxj, kidxi] = M[kidxi, kidxj];
    
    // Invert M
    vw::Matrix<double> Minv = vw::math::inverse(M);
    // Solve for x in Mx = B
    vw::Vector<double> Soln = B * Minv;

    // Worry about translating here...
    kernelCoeffs = Soln;  // somehow

    // clean up
    delete B;
    delete M;
    delete convolvedImageVec;
    delete convolvedAccessorRowVec;
    delete Minv;
    delete Soln;
}

static void getCollectionOfMaskedImagesForPSFMatching(
    vector<Object> &objectCollection ///< Vector of objects to use for diffim kernel
    ) {
    Object obj1 = new Object(100, 100, 10, 10);
    objectCollection.push_back(obj1);
}

// TODO BELOW

static void getTemplateChunkExposureFromTemplateExposure() {
    wcsMatchExposure();
}
static void wcsMatchExposure() {
}




static void computeSpatiallyVaryingPSFMatchingKernel() {
    fitKernelsUsingPrincipalComponentAnalysis();
}
static void fitKernelsUsingPrincipalComponentAnalysis() {
    fitArraysUsingPrincipalComponentAnalysis();
}
static void fitArraysUsingPrincipalComponentAnalysis() {
}


