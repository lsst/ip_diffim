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

#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/KernelFunctions.h>
#include <lsst/fw/PixelAccessors.h>
#include <lsst/fw/Source.h>
#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 
#include <boost/shared_ptr.hpp>

using namespace std;

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
template <typename PixelT, typename MaskT, typename KernelT>
void lsst::imageproc::computePSFMatchingKernelForMaskedImage(
    lsst::fw::MaskedImage<PixelT,MaskT> const &imageToConvolve, ///< Template image; convolved
    lsst::fw::MaskedImage<PixelT,MaskT> const &imageToNotConvolve, ///< Science image; not convolved
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelBasisVec ///< Input set of basis kernels
    ) {

    vector<lsst::fw::Source> sourceCollection;
    getCollectionOfMaskedImagesForPSFMatching(sourceCollection);

    // Reusable view around each source
    typename lsst::fw::MaskedImage<PixelT,MaskT>::MaskedImagePtrT imageToConvolvePtr;
    typename lsst::fw::MaskedImage<PixelT,MaskT>::MaskedImagePtrT imageToNotConvolvePtr;

    // Collection of output kernels
    vector<lsst::fw::LinearCombinationKernel<KernelT> > kernelVec(sourceCollection.size());

    // Iterate over source
    typename vector<lsst::fw::Source>::iterator siter = sourceCollection.begin();
    typename vector<lsst::fw::LinearCombinationKernel<KernelT> >::iterator kiter = kernelVec.begin();

    for (; siter != sourceCollection.end(); ++siter, ++kiter) {
        lsst::fw::Source diffImSource = *siter;
        
        // grab view around each source
        // do i really want a new stamp or just a view?
        BBox2i stamp(boost::any_cast<const int>(diffImSource.getRowc() - diffImSource.getDrow()), 
                     boost::any_cast<const int>(diffImSource.getRowc() + diffImSource.getDrow()),
                     boost::any_cast<const int>(diffImSource.getColc() - diffImSource.getDcol()), 
                     boost::any_cast<const int>(diffImSource.getColc() + diffImSource.getDcol()));
        imageToConvolvePtr    = imageToConvolve.getSubImage(stamp);
        imageToNotConvolvePtr = imageToNotConvolve.getSubImage(stamp);

        vector<double> kernelCoeffs(kernelBasisVec.size());

        // Find best single kernel for this stamp
        lsst::imageproc::computePSFMatchingKernelForPostageStamp
            (*imageToConvolvePtr, *imageToNotConvolvePtr, kernelBasisVec, kernelCoeffs);
        
        // Create a linear combination kernel from this and append to kernelVec
        lsst::fw::LinearCombinationKernel<KernelT> sourceKernel(kernelBasisVec, kernelCoeffs);
        *kiter = sourceKernel;
    }
    // Hold output PCA kernel
    vector<lsst::fw::Kernel<KernelT> > kernelPCABasisVec;
    lsst::imageproc::computePCAKernelBasis(kernelVec, kernelPCABasisVec);
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
template <typename PixelT, typename MaskT, typename KernelT>
void lsst::imageproc::computePSFMatchingKernelForPostageStamp(
    lsst::fw::MaskedImage<PixelT, MaskT> const &imageToConvolve, ///< Goes with the code
    lsst::fw::MaskedImage<PixelT, MaskT> const &imageToNotConvolve, ///< This is for doxygen
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelBasisVec, ///< Input kernel basis set
    vector<double> &kernelCoeffs ///< Output kernel basis coefficients
    ) { 

    int nKernelParameters=0, nBackgroundParameters=0, nParameters=0;
    const float threshold = 0.0;
    
    // We assume that each kernel in the Set has 1 parameter you fit for
    nKernelParameters = kernelBasisVec.size();
    // Or, we just assume that across a single kernel, background 0th order.  This quite makes sense.
    nBackgroundParameters = 1;
    // Total number of parameters
    nParameters = nKernelParameters + nBackgroundParameters;
    
    vw::Vector<double> B(nParameters);
    vw::Matrix<double> M(nParameters, nParameters);

    // convolve creates a MaskedImage, push it onto the back of the Vector
    // need to use shared pointers because MaskedImage copy does not work
    vector<boost::shared_ptr<lsst::fw::MaskedImage<PixelT, MaskT> > > convolvedImageVec;

    typename vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > >::const_iterator kiter = kernelBasisVec.begin();
    typename vector<boost::shared_ptr<lsst::fw::MaskedImage<PixelT, MaskT> > >::iterator citer = convolvedImageVec.begin();

    for (; kiter != kernelBasisVec.end(); ++kiter, ++citer) {

        boost::shared_ptr<lsst::fw::MaskedImage<PixelT, MaskT> > imagePtr(
            new lsst::fw::MaskedImage<PixelT, MaskT>(lsst::fw::kernel::convolve(imageToConvolve, **kiter, threshold, vw::NoEdgeExtension(), -1))
            );

        *citer = imagePtr;
    } 

    // An accessor for each convolution plane
    // NOTE : MaskedPixelAccessor has no empty constructor, therefore we need to push_back()
    vector<lsst::fw::MaskedPixelAccessor<PixelT, MaskT> > convolvedAccessorRowVec;
    for (citer = convolvedImageVec.begin(); citer != convolvedImageVec.end(); ++citer) {
        lsst::fw::MaskedPixelAccessor<PixelT, MaskT> convolvedAccessorRow(**citer);
        convolvedAccessorRowVec.push_back(convolvedAccessorRow);
    }

    // An accessor for each input image
    lsst::fw::MaskedPixelAccessor<PixelT, MaskT> imageToConvolveRow(imageToConvolve);
    lsst::fw::MaskedPixelAccessor<PixelT, MaskT> imageToNotConvolveRow(imageToNotConvolve);

    // integral over image's dx and dy
    for (int row = 0; row < imageToConvolve.getRows(); row++) {

        // An accessor for each convolution plane
        vector<lsst::fw::MaskedPixelAccessor<PixelT, MaskT> > convolvedAccessorColVec;
        for (int ki = 0; ki < nKernelParameters; ki++) {
            lsst::fw::MaskedPixelAccessor<PixelT, MaskT> convolvedAccessorCol = convolvedAccessorRowVec[ki];
            convolvedAccessorColVec.push_back(convolvedAccessorCol);
        }

        // An accessor for each input image
        lsst::fw::MaskedPixelAccessor<PixelT, MaskT> imageToConvolveCol = imageToConvolveRow;
        lsst::fw::MaskedPixelAccessor<PixelT, MaskT> imageToNotConvolveCol = imageToNotConvolveRow;

        for (int col = 0; col < imageToConvolve.getCols(); col++) {
            
            PixelT ncCamera = *imageToNotConvolveCol.image;
            PixelT ncVariance = *imageToNotConvolveCol.variance;

            // Its an approximation to use this since you don't really know the kernel yet
            // You really want the post-convolved variance but this is close enough
            PixelT cVariance = *imageToConvolveCol.variance;
            
            // Quicker computation
            double iVariance = 1.0 / (ncVariance + cVariance);

            // kernel index i
            for (int kidxi = 0; kidxi < nKernelParameters; kidxi++) {
                PixelT cdCamerai = *convolvedAccessorColVec[kidxi].image;
                B[kidxi] += ncCamera * cdCamerai * iVariance;
                
                // kernel index j 
                for (int kidxj = kidxi; kidxj < nKernelParameters; kidxj++) {
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
            for (int ki = 0; ki < nKernelParameters; ki++) {
                convolvedAccessorColVec[ki].nextCol();
            }
        } // col
        
        // Step each accessor in row
        imageToConvolveRow.nextRow();
        imageToNotConvolveRow.nextRow();
        for (int ki = 0; ki < nKernelParameters; ki++) {
            convolvedAccessorRowVec[ki].nextRow();
        }

        // clean up
        convolvedAccessorColVec.~vector<lsst::fw::MaskedPixelAccessor<PixelT, MaskT> >();
    } // row

    // Fill in rest of M
    for (int kidxi=0; kidxi < nKernelParameters; kidxi++) 
        for (int kidxj=kidxi+1; kidxj < nKernelParameters; kidxj++) 
            M[kidxj][kidxi] = M[kidxi][kidxj];
    
    // Invert M
    vw::Matrix<double> Minv = vw::math::inverse(M);
    // Solve for x in Mx = B
    //vw::math::Vector<double> Soln = B * Minv;  // uh, this does not work for B
    vw::math::Vector<double> Soln = Minv * B; // will this at least compile, '*' is a method for Minv

    // Worry about translating here...
    //kernelCoeffs = Soln;  // somehow
    for (int ki = 0; ki < nKernelParameters; ki++) {
        kernelCoeffs[ki] = Soln[ki];
    }

    // clean up
    //delete B;
    //delete M;
    convolvedImageVec.~vector<boost::shared_ptr<lsst::fw::MaskedImage<PixelT, MaskT> > >();
    convolvedAccessorRowVec.~vector<lsst::fw::MaskedPixelAccessor<PixelT, MaskT> >();
    //delete Minv;
    //delete Soln;
}

void lsst::imageproc::getCollectionOfMaskedImagesForPSFMatching(
    vector<lsst::fw::Source> &sourceCollection ///< Vector of sources to use for diffim kernel
    ) {

    lsst::fw::Source src1(100, 100, 10, 10); // Hack some positions here
    sourceCollection.push_back(src1);
}

template <typename KernelT>
void lsst::imageproc::computePCAKernelBasis(
    vector<lsst::fw::LinearCombinationKernel<KernelT> > const &kernelVec, ///< Original input kernel basis set
    vector<lsst::fw::Kernel<KernelT> > &kernelPCABasisVec ///< Output principal components as kernel images
    ) {
    
    int xEval=0;
    int yEval=0;
    bool doNormalize = false;
    int nKernel = kernelVec.size();
    int nPixels = kernelVec[0].getCols() * kernelVec[0].getRows();
    vw::Matrix<double> A(nKernel, nPixels);

    typedef double ImageT;
    typedef typename vw::ImageView<ImageT>::pixel_accessor imageAccessorType;

    typename vector<lsst::fw::LinearCombinationKernel<KernelT> >::const_iterator kiter = kernelVec.begin();
    // do i want kiter++ or ++kiter?
    for (int ki = 0; kiter != kernelVec.end(); ki++, kiter++) {
        lsst::fw::Image<ImageT> kImage = kiter->getImage(xEval, yEval, doNormalize);
        imageAccessorType imageAccessor(kImage.origin());
        
        int nRows = kImage.getRows();
        int nCols = kImage.getCols();
        for (int row = 0; row < nRows; row++) {
            for (int col = 0; col < nCols; col++) {
                A[ki][col + row * nRows] = *imageAccessor;
                imageAccessor.next_col();
            }
            imageAccessor.next_row();
        }
    }
    
    // Use LAPACK SVD
    // http://www.netlib.org/lapack/lug/node32.html
    // Suggests to me that we use DGESVD or DGESDD (preferred)
}

// TODO BELOW

void lsst::imageproc::getTemplateChunkExposureFromTemplateExposure() {
}
void lsst::imageproc::wcsMatchExposure() {
}
void lsst::imageproc::computeSpatiallyVaryingPSFMatchingKernel() {
}


