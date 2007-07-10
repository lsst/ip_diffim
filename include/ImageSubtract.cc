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
#include "lsst/fw/Kernel.h"
#include "lsst/fw/KernelFunctions.h"
#include "lsst/fw/PixelAccessors.h"
#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 

using namespace std;

// Inline Member Functions
inline unsigned
lsst::imageproc::Source::getColc() const {
    return _colc;
}
inline unsigned
lsst::imageproc::Source::getRowc() const {
    return _rowc;
}
inline unsigned
lsst::imageproc::Source::getDcol() const {
    return _dcol;
}
inline unsigned
lsst::imageproc::Source::getDrow() const {
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
template <class PixelT, class MaskT, class KernelT>
void lsst::imageproc::computePSFMatchingKernelForMaskedImage(
    lsst::fw::MaskedImage<PixelT,MaskT> const &imageToConvolve, ///< Template image; convolved
    lsst::fw::MaskedImage<PixelT,MaskT> const &imageToNotConvolve, ///< Science image; not convolved
    lsst::fw::LinearCombinationKernel<KernelT> &kernelBasisSet ///< Input set of basis kernels; modified on return
    ) {

    vector<lsst::imageproc::Source> sourceCollection;
    getCollectionOfMaskedImagesForPSFMatching(sourceCollection);

    // Reusable view around each source
    // typedef typename MaskedImage<PixelT,MaskT>::MaskedImagePtrT maskedImagePtrTType;

    typename lsst::fw::MaskedImage<PixelT,MaskT>::MaskedImagePtrT imageToConvolvePtr;
    typename lsst::fw::MaskedImage<PixelT,MaskT>::MaskedImagePtrT imageToNotConvolvePtr;

    // Output kernels
    // This vector of vectors thing is a bit tricky...
    vector<vector<KernelT> > kernelCoeffsVec( sourceCollection.size(), vector<KernelT>(kernelBasisSet.getNKernelParameters()) );

    // Iterate over source; use iterator instead?
    for (unsigned nobj = 0; nobj < sourceCollection.size(); nobj++) {
        Source diffImSource = sourceCollection[nobj];
        
        // grab view around each source
        // do i really want a new stamp or just a view?
        BBox2i stamp(diffImSource.getRowc() - diffImSource.getDrow(), diffImSource.getRowc() + diffImSource.getDrow(), 
                     diffImSource.getColc() - diffImSource.getDcol(), diffImSource.getColc() + diffImSource.getDcol());
        imageToConvolvePtr    = imageToConvolve.getSubImage(stamp);
        imageToNotConvolvePtr = imageToNotConvolve.getSubImage(stamp);

        // you need to initialize the size of this one
        //std::vector<KernelT> kernelCoeffs(kernelBasisSet.getNKernelParameters());
        computePSFMatchingKernelForPostageStamp(*imageToConvolvePtr, *imageToNotConvolvePtr, kernelBasisSet, kernelCoeffsVec[nobj]);
        //kernelCoeffsVec[nobj] = kernelCoeffs;
    }

    // Does nothing currently
    //computeSpatiallyVaryingPSFMatchingKernel(kernelBasisSet, kernelCoeffsVec);
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
template <class PixelT, class MaskT, class KernelT>
void lsst::imageproc::computePSFMatchingKernelForPostageStamp(
    lsst::fw::MaskedImage<PixelT, MaskT> const &imageToConvolve, ///< Goes with the code
    lsst::fw::MaskedImage<PixelT, MaskT> const &imageToNotConvolve, ///< This is for doxygen
    lsst::fw::LinearCombinationKernel<KernelT> &kernelBasisSet, ///< This is for doxygen
    std::vector<KernelT> &kernelCoeffs ///< This is for doxygen
    ) { 

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
    std::vector<boost::shared_ptr<lsst::fw::Kernel<PixelT> > > kernelList = kernelBasisSet.getKernelList();
    // convolve creates a MaskedImage, push it onto the back of the Vector
    vector<lsst::fw::MaskedImage<PixelT, MaskT> > convolvedImageVec;

    // Non-iterator
    //    for (unsigned ki = 0; ki < nKernelParameters; ki++) {
    //        lsst::fw::MaskedImage<PixelT, MaskT>
    //            convolvedImage = lsst::fw::kernel::convolve(imageToConvolve, *(kernelList[ki]), 0.0, vw::NoEdgeExtension());
    //        convolvedImageVec.push_back(convolvedImage);
    //    }

    // iterator
    typename std::vector<boost::shared_ptr<lsst::fw::Kernel<PixelT> > >::iterator kiter = kernelList.begin();
    for (; kiter != kernelList.end(); ++kiter) {
        lsst::fw::MaskedImage<PixelT, MaskT>
            convolvedImage = lsst::fw::kernel::convolve(imageToConvolve, *kiter, 0.0, vw::NoEdgeExtension());
        convolvedImageVec.push_back(convolvedImage);
    } 

    // An accessor for each convolution plane
    // NOTE : MaskedPixelAccessor has no empty constructor, therefore we need to push_back()
    vector<lsst::fw::MaskedPixelAccessor<PixelT, MaskT> > convolvedAccessorRowVec;
    for (unsigned ki = 0; ki < nKernelParameters; ki++) {
        lsst::fw::MaskedPixelAccessor<PixelT, MaskT> convolvedAccessorRow(convolvedImageVec[ki]);
        convolvedAccessorRowVec.push_back(convolvedAccessorRow);
    }

    // An accessor for each input image
    lsst::fw::MaskedPixelAccessor<PixelT, MaskT> imageToConvolveRow(imageToConvolve);
    lsst::fw::MaskedPixelAccessor<PixelT, MaskT> imageToNotConvolveRow(imageToNotConvolve);

    // integral over image's dx and dy
    for (unsigned row = 0; row < imageToConvolve.getRows(); row++) {

        // An accessor for each convolution plane
        vector<lsst::fw::MaskedPixelAccessor<PixelT, MaskT> > convolvedAccessorColVec;
        for (unsigned ki = 0; ki < nKernelParameters; ki++) {
            lsst::fw::MaskedPixelAccessor<PixelT, MaskT> convolvedAccessorCol = convolvedAccessorRowVec[ki];
            convolvedAccessorColVec.push_back(convolvedAccessorCol);
        }

        // An accessor for each input image
        lsst::fw::MaskedPixelAccessor<PixelT, MaskT> imageToConvolveCol = imageToConvolveRow;
        lsst::fw::MaskedPixelAccessor<PixelT, MaskT> imageToNotConvolveCol = imageToNotConvolveRow;

        for (unsigned col = 0; col < imageToConvolve.getCols(); col++) {
            
            PixelT ncCamera = *imageToNotConvolveCol.image;
            PixelT ncVariance = *imageToNotConvolveCol.variance;

            // Its an approximation to use this since you don't really know the kernel yet
            // You really want the post-convolved variance but this is close enough
            PixelT cVariance = *imageToConvolveCol.variance;
            
            // Quicker computation
            double iVariance = 1.0 / (ncVariance + cVariance);

            // kernel index i
            for (unsigned kidxi = 0; kidxi < nKernelParameters; kidxi++) {
                PixelT cdCamerai = *convolvedAccessorColVec[kidxi].image;
                B[kidxi] += ncCamera * cdCamerai * iVariance;
                
                // kernel index j 
                for (unsigned kidxj = kidxi; kidxj < nKernelParameters; kidxj++) {
                    PixelT cdCameraj = *convolvedAccessorColVec[kidxj].image;
                    M[kidxi][kidxj] += cdCamerai * cdCameraj * iVariance;
                } // kidxj
            } // kidxi

            // Here we have a single background term
            B[nParameters] += ncCamera * iVariance;
            M[nParameters][nParameters] += 1.0 * iVariance;

            // Step each accessor in column
            imageToConvolveCol.nextCol();
            imageToNotConvolveCol.nextCol();
            for (unsigned ki = 0; ki < nKernelParameters; ki++) {
                convolvedAccessorColVec[ki].nextCol();
            }
        } // col
        
        // Step each accessor in row
        imageToConvolveRow.nextRow();
        imageToNotConvolveRow.nextRow();
        for (unsigned ki = 0; ki < nKernelParameters; ki++) {
            convolvedAccessorRowVec[ki].nextRow();
        }

        // clean up
        convolvedAccessorColVec.~vector<lsst::fw::MaskedPixelAccessor<PixelT, MaskT> >();
    } // row

    // Fill in rest of M
    for (unsigned kidxi=0; kidxi < nKernelParameters; kidxi++) 
        for (unsigned kidxj=kidxi+1; kidxj < nKernelParameters; kidxj++) 
            M[kidxj][kidxi] = M[kidxi][kidxj];
    
    // Invert M
    vw::Matrix<double> Minv = vw::math::inverse(M);
    // Solve for x in Mx = B
    //vw::math::Vector<double> Soln = B * Minv;  // uh, this does not work for B
    vw::math::Vector<double> Soln = Minv * B; // will this at least compile, '*' is a method for Minv

    // Worry about translating here...
    //kernelCoeffs = Soln;  // somehow
    for (unsigned ki = 0; ki < nKernelParameters; ki++) {
        kernelCoeffs[ki] = Soln[ki];
    }

    // clean up
    //delete B;
    //delete M;
    convolvedImageVec.~vector<lsst::fw::MaskedImage<PixelT, MaskT> >();
    convolvedAccessorRowVec.~vector<lsst::fw::MaskedPixelAccessor<PixelT, MaskT> >();
    //delete Minv;
    //delete Soln;
}

void lsst::imageproc::getCollectionOfMaskedImagesForPSFMatching(
    vector<lsst::imageproc::Source> &sourceCollection ///< Vector of sources to use for diffim kernel
    ) {
    lsst::imageproc::Source src1(100, 100, 10, 10);
    sourceCollection.push_back(src1);
}

// TODO BELOW

void lsst::imageproc::getTemplateChunkExposureFromTemplateExposure() {
    wcsMatchExposure();
}
void lsst::imageproc::wcsMatchExposure() {
}
void lsst::imageproc::computeSpatiallyVaryingPSFMatchingKernel() {
    fitKernelsUsingPrincipalComponentAnalysis();
}
void lsst::imageproc::fitKernelsUsingPrincipalComponentAnalysis() {
    fitArraysUsingPrincipalComponentAnalysis();
}
void lsst::imageproc::fitArraysUsingPrincipalComponentAnalysis() {
}


