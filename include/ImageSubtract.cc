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
#include <lsst/fw/MinimizerFunctionBase.h>
#include <lsst/mwi/utils/Trace.h>
#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 
#include <vw/Math/LinearAlgebra.h> 
#include <vw/Math/Functions.h> 
#include <boost/shared_ptr.hpp>
#include <PCA.h>
#include <Minuit/MnMigrad.h>
#include <Minuit/MnMinos.h>
#include <Minuit/MnPrint.h>
#include <Minuit/FunctionMinimum.h>
#include <string> // for upar.add
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
template <typename ImageT, typename MaskT, typename KernelT>
void lsst::imageproc::computePSFMatchingKernelForMaskedImage(
    lsst::fw::MaskedImage<ImageT,MaskT> const &imageToConvolve, ///< Template image; convolved
    lsst::fw::MaskedImage<ImageT,MaskT> const &imageToNotConvolve, ///< Science image; not convolved
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelBasisVec, ///< Input set of basis kernels
    boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > kernelPtr, ///< The output convolution kernel
    boost::shared_ptr<lsst::fw::function::Function2<KernelT> > backgroundFunctionPtr ///< Function for spatial variation of background

    vector<lsst::fw::Source> sourceCollection;
    getCollectionOfMaskedImagesForPSFMatching(sourceCollection);

    lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForMaskedImage", 2, "Entering subroutine computePSFMatchingKernelForMaskedImage");

    // Reusable view around each source
    typename lsst::fw::MaskedImage<ImageT,MaskT>::MaskedImagePtrT imageToConvolvePtr;
    typename lsst::fw::MaskedImage<ImageT,MaskT>::MaskedImagePtrT imageToNotConvolvePtr;

    // Collection of output kernels
    vector<lsst::fw::LinearCombinationKernel<KernelT> > kernelVec(sourceCollection.size());

    // Iterate over source
    typename vector<lsst::fw::Source>::iterator siter = sourceCollection.begin();
    typename vector<lsst::fw::LinearCombinationKernel<KernelT> >::iterator kiter = kernelVec.begin();

    // Vector for backgrounds
    vector<double> backgrounds;

    for (; siter != sourceCollection.end(); ++siter, ++kiter) {
        lsst::fw::Source diffImSource = *siter;
        
        // grab view around each source
        // do i really want a new stamp or just a view?
        // Bbox2i has x y cols rows
        // NOTE : we need to make sure we get the centering right with these +1, etc...
        BBox2i stamp(floor(diffImSource.getColc() - diffImSource.getDcol()), 
                     floor(diffImSource.getRowc() - diffImSource.getDrow()), 
                     ceil(2 * diffImSource.getDcol() + 1),
                     ceil(2 * diffImSource.getDrow() + 1)
                     );

        imageToConvolvePtr    = imageToConvolve.getSubImage(stamp);
        imageToNotConvolvePtr = imageToNotConvolve.getSubImage(stamp);

        vector<double> kernelCoeffs(kernelBasisVec.size());

        imageToConvolvePtr->writeFits( (boost::format("iFits_%d") % diffImSource.getSourceId()).str() );

        // Find best single kernel for this stamp
        double background;
        lsst::imageproc::computePSFMatchingKernelForPostageStamp
            (*imageToConvolvePtr, *imageToNotConvolvePtr, kernelBasisVec, kernelCoeffs, background);
        backgrounds.push_back(background);

        // Create a linear combination kernel from this and append to kernelVec
        lsst::fw::LinearCombinationKernel<KernelT> sourceKernel(kernelBasisVec, kernelCoeffs);
        *kiter = sourceKernel;

        // NOTE
        // QA - calculate the residual of the subtracted image here

        // DEBUGGING
        double imSum;
        lsst::fw::Image<KernelT> kImage = sourceKernel.computeNewImage(imSum);
        kImage.writeFits( (boost::format("kFits_%d.fits") % diffImSource.getSourceId()).str() );
        

    }

    vector<double> kernelResidualsVec;
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > kernelBasisVec;
    vw::math::Matrix<double> kernelCoefficients;

    if (type == DeltaFunction) {
        lsst::imageproc::computePCAKernelBasis(kernelVec, kernelResidualsVec, kernelBasisVec, kernelCoefficients);
    }

    // From policy
    const unsigned int spatialOrder = 2;

    boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > spatiallyVaryingKernelPtr(
        new lsst::fw::LinearCombinationKernel<KernelT>
        );
    // NOTE - do we want KernelT or another template FuncT?
    boost::shared_ptr<lsst::fw::function::Function2<KernelT> > spatiallyVaryingFunctionPtr(
        new lsst::fw::function::PolynomialFunction2<KernelT>(spatialOrder)
        );

    if (kernelPtr->getNSpatialParameters() > 1) {
        computeSpatiallyVaryingPSFMatchingKernel(kernelBasisVec, 
                                                 kernelCoefficients, 
                                                 sourceCollection,
                                                 spatiallyVaryingFunctionPtr, 
                                                 spatiallyVaryingKernelPtr);
    }

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
template <typename ImageT, typename MaskT, typename KernelT>
void lsst::imageproc::computePSFMatchingKernelForPostageStamp(
    lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve, ///< Goes with the code
    lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve, ///< This is for doxygen
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelBasisVec, ///< Input kernel basis set
    vector<double> &kernelCoeffs, ///< Output kernel basis coefficients
    double &background ///< Difference in the backgrounds
    ) { 

    int nKernelParameters=0, nBackgroundParameters=0, nParameters=0;
    const KernelT threshold = 0.0;

    lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 2, "Entering subroutine computePSFMatchingKernelForPostageStamp");
    
    // We assume that each kernel in the Set has 1 parameter you fit for
    nKernelParameters = kernelBasisVec.size();
    // Or, we just assume that across a single kernel, background 0th order.  This quite makes sense.
    nBackgroundParameters = 1;
    // Total number of parameters
    nParameters = nKernelParameters + nBackgroundParameters;
    
    vw::math::Vector<double> B(nParameters);
    vw::math::Matrix<double> M(nParameters, nParameters);
    for (unsigned int i = nParameters; i--;) {
        B(i) = 0;
        for (unsigned int j = nParameters; j--;) {
            M(i,j) = 0;
        }
    }

    // convolve creates a MaskedImage, push it onto the back of the Vector
    // need to use shared pointers because MaskedImage copy does not work
    vector<boost::shared_ptr<lsst::fw::MaskedImage<ImageT, MaskT> > > convolvedImageVec(nKernelParameters);
    // and an iterator over this
    typename vector<boost::shared_ptr<lsst::fw::MaskedImage<ImageT, MaskT> > >::iterator citer = convolvedImageVec.begin();
    
    // Iterator for input kernel basis
    typename vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > >::const_iterator kiter = kernelBasisVec.begin();
    int kId = 0; // debugging
    // Create C_ij in the formalism of Alard & Lupton
    for (; kiter != kernelBasisVec.end(); ++kiter, ++citer, ++kId) {

        lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 3, "Convolving an Object with Basis");
        
        // NOTE : we could also *precompute* the entire template image convolved with these functions
        //        and save them somewhere to avoid this step each time.  however, our paradigm is to
        //        compute whatever is needed on the fly.  hence this step here.
        boost::shared_ptr<lsst::fw::MaskedImage<ImageT, MaskT> > imagePtr(
            new lsst::fw::MaskedImage<ImageT, MaskT>(lsst::fw::kernel::convolve(imageToConvolve, **kiter, threshold, 1))
            );

        lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 3, "Convolved an Object with Basis");

        *citer = imagePtr;
        
        imagePtr->writeFits( (boost::format("cFits_%d") % kId).str() );
    } 

    // NOTE ABOUT CONVOLUTION :
    // getCtrCol:getCtrRow pixels are masked on the left:bottom side
    // getCols()-getCtrCol():getRows()-getCtrRow() masked on right/top side
    // 
    // The convolved image and the input image are by default the same size, so
    // we offset our initial pixel references by the same amount
    //
    // NOTE : the last *good* pixel is endCol - 1
    //        in your for loops, therefore say "col < endCol" not "col <= endCol"
    kiter = kernelBasisVec.begin();
    citer = convolvedImageVec.begin();
    unsigned int startCol = (*kiter)->getCtrCol();
    unsigned int startRow = (*kiter)->getCtrRow();
    unsigned int endCol   = (*citer)->getCols() - ((*kiter)->getCols() - (*kiter)->getCtrCol());
    unsigned int endRow   = (*citer)->getRows() - ((*kiter)->getRows() - (*kiter)->getCtrRow());
    // NOTE - we need to enforce that the input images are large enough
    // How about some multiple of the PSF FWHM?  Or second moments?

    // An accessor for each convolution plane
    // NOTE : MaskedPixelAccessor has no empty constructor, therefore we need to push_back()
    vector<lsst::fw::MaskedPixelAccessor<ImageT, MaskT> > convolvedAccessorRowVec;
    for (citer = convolvedImageVec.begin(); citer != convolvedImageVec.end(); ++citer) {
        convolvedAccessorRowVec.push_back(lsst::fw::MaskedPixelAccessor<ImageT, MaskT>(**citer));
    }

    // An accessor for each input image; address rows and cols separately
    lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imageToConvolveRow(imageToConvolve);
    lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imageToNotConvolveRow(imageToNotConvolve);

    // Take into account buffer for kernel images
    imageToConvolveRow.advance(startCol, startRow);
    imageToNotConvolveRow.advance(startCol, startRow);
    for (int ki = 0; ki < nKernelParameters; ki++) {
        convolvedAccessorRowVec[ki].advance(startCol, startRow);
    }

    // integral over image's dx and dy
    for (unsigned int row = startRow; row < endRow; row++) {

        // An accessor for each convolution plane
        vector<lsst::fw::MaskedPixelAccessor<ImageT, MaskT> > convolvedAccessorColVec = convolvedAccessorRowVec;

        // An accessor for each input image; places the col accessor at the correct row
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imageToConvolveCol = imageToConvolveRow;
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imageToNotConvolveCol = imageToNotConvolveRow;

        for (unsigned int col = startCol; col < endCol; col++) {

            ImageT ncCamera = *imageToNotConvolveCol.image;
            ImageT ncVariance = *imageToNotConvolveCol.variance;
            MaskT  ncMask = *imageToNotConvolveCol.mask;

            ImageT cVariance = *imageToConvolveCol.variance;

            // Do we skip these?
            if (ncMask & 0x1) {
                continue;
            }

            lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 4, 
                                    boost::format("Accessing image row %d col %d : %f %f") % row % col % ncCamera % ncVariance);

            // kernel index i
            typename vector<lsst::fw::MaskedPixelAccessor<ImageT, MaskT> >::iterator kiteri = convolvedAccessorColVec.begin();
            for (int kidxi = 0; kiteri != convolvedAccessorColVec.end(); kiteri++, kidxi++) {
                ImageT cdCamerai = *kiteri->image;
                ImageT cdVariancei = *kiteri->variance;
                MaskT cdMaski = *kiteri->mask;
                // Do we skip these?
                if (cdMaski & 0x1) {
                    continue;
                }
                lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 5, 
                                        boost::format("Accessing convolved image %d : %f %f") % kidxi % cdCamerai % cdVariancei);

                B[kidxi] += ncCamera * cdCamerai / (ncVariance + cdVariancei);
                
                // kernel index j 
                typename vector<lsst::fw::MaskedPixelAccessor<ImageT, MaskT> >::iterator kiterj = kiteri;
                for (int kidxj = kidxi; kiterj != convolvedAccessorColVec.end(); kiterj++, kidxj++) {
                    ImageT cdCameraj = *kiterj->image;
                    ImageT cdVariancej = *kiterj->variance;
                    MaskT cdMaskj = *kiterj->mask;
                    // Do we skip these?
                    if (cdMaskj & 0x1) {
                        continue;
                    }
                    lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 6, 
                                            boost::format("Accessing convolved image %d : %f %f") % kidxj % cdCameraj % cdVariancej);
                    M[kidxi][kidxj] += cdCamerai * cdCameraj / (cdVariancei + cdVariancej);
                } 
            } 

            // Here we have a single background term
            B[nParameters-1] += ncCamera / (ncVariance + cVariance); 
            M[nParameters-1][nParameters-1] += 1.0 / (ncVariance + cVariance);
            lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 4, 
                                    boost::format("Background terms : %f %f") %B[nParameters-1] % M[nParameters-1][nParameters-1]);

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
        
    } // row

    // Fill in rest of M
    for (int kidxi=0; kidxi < nKernelParameters; kidxi++) 
        for (int kidxj=kidxi+1; kidxj < nKernelParameters; kidxj++) 
            M[kidxj][kidxi] = M[kidxi][kidxj];

    cout << "B : " << B << endl;
    cout << "M : " << M << endl;

    // Invert M
    vw::math::Matrix<double> Minv = vw::math::inverse(M);

    // Solve for x in Mx = B
    vw::math::Vector<double> Soln = Minv * B;

    // Worry about translating here...
    for (int ki = 0; ki < nKernelParameters; ki++) {
        kernelCoeffs[ki] = Soln[ki];
    }
    background = Soln[nParameters-1];
}

void lsst::imageproc::getCollectionOfMaskedImagesForPSFMatching(
    vector<lsst::fw::Source> &sourceCollection ///< Vector of sources to use for diffim kernel
    ) {

    // Hack some positions in for /lsst/becker/lsst_devel/DC2/fw/tests/data/871034p_1_MI_img.fits
    lsst::fw::Source src1(1, 1010.345, 2375.548, 10., 10.); 
    lsst::fw::Source src2(2, 404.248, 573.398, 10., 10.);
    lsst::fw::Source src3(3, 1686.743, 1880.935, 10., 10.);

    sourceCollection.push_back(src1);
    sourceCollection.push_back(src2);
    sourceCollection.push_back(src3);
}

template <typename KernelT>
void lsst::imageproc::computePCAKernelBasis(
    vector<lsst::fw::LinearCombinationKernel<KernelT> > const &kernelVec, ///< Original input kernel basis set
    vector<double> &kernelResidualsVec, ///< Output residuals of reconstructed kernel
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelPCABasisVec, ///< Output principal components as kernel images
    vw::math::Matrix<double> &kernelCoefficients ///< Output coefficients for each basis function for each kernel
    ) {
    
    //typedef double CalcT;
    //typedef float ImageT;  // Watch out with multiple ImageT definitions!

    const int nKernel = kernelVec.size();
    const int nKCols = kernelVec[0].getCols();
    const int nKRows = kernelVec[0].getRows();
    const int nPixels = nKCols * nKRows;

    lsst::mwi::utils::Trace("lsst.imageproc.computePCAKernelBasis", 2, "Entering subroutine computePCAKernelBasis");

    // Matrix to invert.  Number of rows = number of pixels; number of columns = number of kernels
    // All calculations here are in double
    vw::Matrix<double> M(nPixels, nKernel); 

    typedef typename vw::ImageView<KernelT>::pixel_accessor imageAccessorType;

    // iterator over kernels
    typename vector<lsst::fw::LinearCombinationKernel<KernelT> >::const_iterator kiter = kernelVec.begin();

    // fill up matrix for PCA
    // does it matter if i use kiter++ or ++kiter?
    double imSum;
    for (int ki = 0; kiter != kernelVec.end(); ki++, kiter++) {
        lsst::fw::Image<KernelT> kImage = kiter->computeNewImage(imSum);

        //assert(nKRows == kImage.getRows());
        //assert(nKCols == kImage.getCols());
        int mIdx = 0;

        imageAccessorType imageAccessorCol(kImage.origin());
        for (int col = 0; col < nKCols; col++) {

            imageAccessorType imageAccessorRow(imageAccessorCol);
            for (int row = 0; row < nKRows; row++, mIdx++) {

                // NOTE : 
                //   arguments to matrix-related functions are given in row,col
                //   arguments to image-related functions are given in col,row
                // HOWEVER :
                //   it doesn't matter which order we put the kernel elements into the PCA
                //   since they are not correlated 
                //   as long as we do the same thing when extracting the components
                // UNLESS :
                //   we want to put in some weighting/regularlization into the PCA
                //   not sure if that is even possible...

                M(mIdx, ki) = *imageAccessorRow;
                imageAccessorRow.next_row();
            }
            imageAccessorCol.next_col();
        }
    }

    vw::math::Matrix<double> eVec(nPixels, nKernel);
    vw::math::Vector<double> eVal(nKernel);
    vw::math::Vector<double> mMean(nPixels);

    // M is mean-subtracted if subtractMean = true
    lsst::mwi::utils::Trace("lsst.imageproc.computePCAKernelBasis", 4, "Computing pricipal components");
    lsst::imageproc::computePCA(M, mMean, eVal, eVec, true);
    lsst::mwi::utils::Trace("lsst.imageproc.computePCAKernelBasis", 4, "Computed pricipal components");

    // We now have the basis functions determined
    // Next determine the coefficients that go in front of all of the individual kernels
    // The size of the coefficients will be rows=M.cols() cols=M.rows()
    // We save the first entry of kernelCoefficients for the Mean Image
    for (unsigned int i = 0; i < M.cols(); i++) {
        kernelCoefficients(0,i) = 1;
    }
    vw::math::Matrix<double> subMatrix = vw::math::submatrix(kernelCoefficients, 0, 1, M.rows(), M.cols()-1);
    lsst::imageproc::decomposeMatrixUsingBasis(M, eVec, subMatrix);

    // Write results into kernelCoefficients
    for (unsigned int i = 1; i < M.rows(); i++) {
        vw::math::select_row(kernelCoefficients, i) = vw::math::select_row(subMatrix, i);
    }

    // We next do quality control here; we reconstruct the input kernels with the truncated basis function set
    // From the policy
    const unsigned int nCoeff = 3; 
    vw::math::Matrix<double> approxM(nPixels, nKernel); 
    lsst::imageproc::approximateMatrixUsingBasis(eVec, subMatrix, nCoeff, approxM);

    for (unsigned int i = 0; i < M.cols(); i++) {
        double residual = 0;
        for (unsigned int j = 0; j < M.rows(); j++) {
            residual += M(i,j) - approxM(i,j);
        }
        kernelResidualsVec[i] = residual;
        cout << " Residual " << i << " " << residual << endl;
    }

    // Turn the Mean Image into a Kernel
    lsst::fw::Image<KernelT> meanImage(nKCols, nKRows);
    imageAccessorType imageAccessorCol(meanImage.origin());
    int mIdx = 0;
    for (int col = 0; col < nKCols; col++) {
        imageAccessorType imageAccessorRow(imageAccessorCol);
        for (int row = 0; row < nKRows; row++, mIdx++) {
            *imageAccessorRow = mMean(mIdx);
            imageAccessorRow.next_row();
        }
        imageAccessorCol.next_col();
    }
    // The mean image is the first member of kernelPCABasisVec
    kernelPCABasisVec.push_back(boost::shared_ptr<lsst::fw::Kernel<KernelT> > (new lsst::fw::FixedKernel<KernelT>(meanImage)));

    // Debugging
    meanImage.writeFits( (boost::format("mFits.fits")).str() );

    // Turn each eVec into an Image and then into a Kernel
    for (unsigned int ki = 0; ki < eVec.cols(); ki++) {
        lsst::fw::Image<KernelT> basisImage(nKCols, nKRows);

        // Not sure how to bulk load information into Image
        int kIdx = 0;

        imageAccessorType imageAccessorCol(basisImage.origin());
        for (int col = 0; col < nKCols; col++) {
            
            imageAccessorType imageAccessorRow(imageAccessorCol);
            for (int row = 0; row < nKRows; row++, kIdx++) {

                *imageAccessorRow = eVec(kIdx, ki);
                imageAccessorRow.next_row();
            }
            imageAccessorCol.next_col();
        }
        // Add to kernel basis
        kernelPCABasisVec.push_back(boost::shared_ptr<lsst::fw::Kernel<KernelT> > (new lsst::fw::FixedKernel<KernelT>(basisImage)));

        // Debugging info
        basisImage.writeFits( (boost::format("eFits_%d.fits") % ki).str() );
    }
}

template <typename KernelT, typename ReturnT>
void lsst::imageproc::computeSpatiallyVaryingPSFMatchingKernel(
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const kernelBasisVec, ///< Input basis kernel set
    vw::math::Matrix<double> const kernelCoefficients, ///< Basis coefficients for all kernels 
    vector<double> const backgrounds, ///< Background measurements for all kernels
    vector<lsst::fw::Source> const sourceCollection, ///< Needed right now for the centers of the kernels in the image; should be able to avoid this 
    boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > &spatiallyVaryingKernelPtr, ///< Output kernel
    boost::shared_ptr<lsst::fw::function::Function2<KernelT> > &backgroundFunctionPtr ///< Output background model
    )
 {
     // NOTE - For a delta function basis set, the mean image is the first entry in kernelPCABasisVec.  
     // This should *not* vary spatially.  Should we code this up or let the fitter determine that it does not vary..?
     // I think I would like some sort of typeid() discrimination here.  For later...
     typename vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > >::const_iterator kiter = kernelBasisVec.begin();

     // Hold the fit parameters
     vector<double> nParameters(spatiallyVaryingFunctionPtr->getNParameters());
     vector<vector<double> > fitParameters(kernelBasisVec.size(), nParameters);

     unsigned int ncoeff = 0;
     for (; kiter != kernelBasisVec.end(); ncoeff++, kiter++) {
         
         vector<double> measurements;
         vector<double> variances;
         vector<double> position1;
         vector<double> position2;

         typename vector<lsst::fw::Source>::const_iterator siter = sourceCollection.begin();

         unsigned int nsource = 0;
         for (; siter != sourceCollection.end(); nsource++, siter++) {
             // NOTE - need to make this work!
             //position1.push_back(siter.colc_norm()); // between -1 and 1
             //position2.push_back(siter.rowc_norm()); // between -1 and 1
             measurements.push_back(kernelCoefficients(nsource, ncoeff));
             //variances.push_back(what);
         }

         double def = 1.0;
         lsst::fw::function::MinimizerFunctionBase2<KernelT> 
             myFcn(measurements, variances, position1, position2, def, spatiallyVaryingFunctionPtr);

         // Initialize paramters
         // Name; value; uncertainty
         MnUserParameters upar;
         for (unsigned int i = 0; i < spatiallyVaryingFunctionPtr->getNParameters(); i++) {
             upar.add((boost::format("p%d") % i).str().c_str(), 0.0, 0.1);
         }
         
         MnMigrad migrad(myFcn, upar);
         FunctionMinimum min = migrad();
         //MnMinos minos(myFcn, min); // only if you want serious uncertainties

         for (unsigned int i = 0; i < spatiallyVaryingFunctionPtr->getNParameters(); i++) {
             fitParameters[ncoeff][i] = min.userState().value((boost::format("p%d") % i).str().c_str());
         }
     }
     // Set up the spatially varying kernel and we are done!
     lsst::fw::LinearCombinationKernel<KernelT> spatiallyVaryingKernel(kernelBasisVec, spatiallyVaryingFunctionPtr, fitParameters); 
}

template <typename KernelT>
void lsst::imageproc::generateDeltaFunctionKernelSet(
    unsigned int const nRows, ///< Number of rows in the kernel basis
    unsigned int const nCols, ///< Number of colunms in the kernel basis
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelBasisVec ///< Output kernel basis function, length nRows * nCols
    )
{
    int colCtr = (nCols - 1) / 2;
    int rowCtr = (nRows - 1) / 2;
    for (unsigned row = 0; row < nRows; ++row) {
        int y = static_cast<int>(row) - rowCtr;
        
        for (unsigned col = 0; col < nCols; ++col) {
            int x = static_cast<int>(col) - colCtr;
            
            typename lsst::fw::Kernel<KernelT>::KernelFunctionPtrType kfuncPtr(
                new lsst::fw::function::IntegerDeltaFunction2<KernelT>(x, y)
                );
            
            boost::shared_ptr<lsst::fw::Kernel<KernelT> > kernelPtr(
                new lsst::fw::AnalyticKernel<KernelT>(kfuncPtr, nCols, nRows)
                );
            
            kernelBasisVec.push_back(kernelPtr);
        }
    }
}

template <typename KernelT>
void lsst::imageproc::generateAlardLuptonKernelSet(
    int const nRows, ///< Number of rows in the kernel basis
    int const nCols, ///< Number of columns in the kernel basis
    vector<double> const sigGauss, ///< Width of gaussians in basis; size = number of Gaussians
    vector<double> const degGauss, ///< Degree of spatial variation within each Gaussian; size = sigGauss.size()
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelBasisVec ///< Output kernel basis function, length sum_n 0.5*(deg_n+1)*(deg_n+2)
    )
{
    // TO DO 
}
