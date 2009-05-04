// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Implementation of image subtraction functions declared in ImageSubtract.h
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */
#include <iostream>
#include <limits>
#include <boost/timer.hpp> 

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>

// NOTE -  trace statements >= 6 can ENTIRELY kill the run time
// #define LSST_MAX_TRACE 5

#if !defined(USE_VW)
#   define USE_VW 0
#endif
#if USE_VW
#include <vw/Math/Functions.h> 
#include <vw/Math/Vector.h> 
#include <vw/Math/Matrix.h> 
#include <vw/Math/LinearAlgebra.h> 
#endif

#include <lsst/ip/diffim/ImageSubtract.h>
#include <lsst/afw/image.h>
#include <lsst/afw/math.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/logging/Log.h>
#include <lsst/afw/detection/Footprint.h>
#include <lsst/afw/math/ConvolveImage.h>

#define DEBUG_MATRIX 0

namespace exceptions = lsst::pex::exceptions; 
namespace logging    = lsst::pex::logging; 
namespace image      = lsst::afw::image;
namespace math       = lsst::afw::math;
namespace diffim     = lsst::ip::diffim;

//
// Constructors
//
template <typename ImageT, typename VarT>
diffim::PsfMatchingFunctor<ImageT, VarT>::PsfMatchingFunctor(
    math::KernelList<lsst::afw::math::Kernel> const& basisList
    ) :
    _basisList(basisList),
    _background(0.),
    _backgroundError(0.),
    _kernel(boost::shared_ptr<lsst::afw::math::Kernel>()),
    _kernelError(boost::shared_ptr<lsst::afw::math::Kernel>())
{;}

//
// Public Member Functions
//

template <typename ImageT, typename VarT>
void diffim::PsfMatchingFunctor<ImageT, VarT>::reset() {
    /* HEY , FOR SOME REASON THE KERNEL RESET DOES NOT WORK AND SEG FAULTS */
    //this->_background      = 0.;
    //this->_backgroundError = 0.;
    //this->_kernel.reset();
    //this->_kernelError.reset();
}

/** Create PSF matching kernel
 */
template <typename ImageT, typename VarT>
void diffim::PsfMatchingFunctor<ImageT, VarT>::apply(
    image::Image<ImageT> const& imageToConvolve,        //!< Image to apply kernel to
    image::Image<ImageT> const& imageToNotConvolve,     //!< Image whose PSF you want to match to
    image::Image<VarT>   const& varianceEstimate,       //!< Estimate of the variance per pixel
    lsst::pex::policy::Policy  const& policy            //!< Policy file
    ) {
    
    // Make sure you do not overwrite anyone else's kernels
    this->reset();

    int const edgeMaskBit = image::Mask<unsigned short>::getMaskPlane("EDGE");
    
    int const nKernelParameters     = this->_basisList.size();
    int const nBackgroundParameters = 1;
    int const nParameters           = nKernelParameters + nBackgroundParameters;
    
    boost::timer t;
    t.restart();
    
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(nParameters, nParameters);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(nParameters);
    
    std::vector<boost::shared_ptr<image::Image<ImageT> > > convolvedImageList(nKernelParameters);
    typename std::vector<boost::shared_ptr<image::Image<ImageT> > >::iterator citer = convolvedImageList.begin();
    std::vector<boost::shared_ptr<math::Kernel> >::const_iterator kiter = this->_basisList.begin();
    
    // Create C_ij in the formalism of Alard & Lupton */
    for (; kiter != this->_basisList.end(); ++kiter, ++citer) {
        /*
         * NOTE : we could also *precompute* the entire template image convolved with these functions
         *        and save them somewhere to avoid this step each time.  however, our paradigm is to
         *        compute whatever is needed on the fly.  hence this step here.
         */
        *citer = typename image::Image<ImageT>::Ptr(new image::Image<ImageT>(imageToConvolve.getDimensions()));
        math::convolve(**citer, imageToConvolve, **kiter, false, edgeMaskBit);
    } 

    double time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                       "Total compute time to do basis convolutions : %.2f s", time);
    t.restart();
     
    kiter = this->_basisList.begin();
    citer = convolvedImageList.begin();

    // Ignore buffers around edge of convolved images :
    //
    // If the kernel has width 5, it has center pixel 2.  The first good pixel
    // is the (5-2)=3rd pixel, which is array index 2, and ends up being the
    // index of the central pixel.
    //
    // You also have a buffer of unusable pixels on the other side, numbered
    // width-center-1.  The last good usable pixel is N-width+center+1.

    // Example : the kernel is width = 5, center = 2
    //
    //     ---|---|-c-|---|---|
    //          
    //           the image is width = N
    //           convolve this with the kernel, and you get
    //
    //    |-x-|-x-|-g-|---|---| ... |---|---|-g-|-x-|-x-|
    //
    //           g = first/last good pixel
    //           x = bad
    // 
    //           the first good pixel is the array index that has the value "center", 2
    //           the last good pixel has array index N-(5-2)+1
    //           eg. if N = 100, you want to use up to index 97
    //               100-3+1 = 98, and the loops use i < 98, meaning the last
    //               index you address is 97.
   
    int const startCol = (*kiter)->getCtrX();
    int const startRow = (*kiter)->getCtrY();
    int const endCol   = (*citer)->getWidth()  - ((*kiter)->getWidth()  - (*kiter)->getCtrX()) + 1;
    int const endRow   = (*citer)->getHeight() - ((*kiter)->getHeight() - (*kiter)->getCtrY()) + 1;

    std::vector<typename image::Image<ImageT>::xy_locator> convolvedLocatorList;
    for (citer = convolvedImageList.begin(); citer != convolvedImageList.end(); ++citer) {
        convolvedLocatorList.push_back( (**citer).xy_at(startCol,startRow) );
    }
    typename image::Image<ImageT>::xy_locator imageToConvolveLocator = imageToConvolve.xy_at(startCol, startRow);
    typename image::Image<ImageT>::xy_locator imageToNotConvolveLocator = imageToNotConvolve.xy_at(startCol, startRow);
    xyi_locator varianceLocator           = varianceEstimate.xy_at(startCol, startRow);

    // Unit test ImageSubtract_1.py should show
    // Image range : 9 9 -> 31 31 : 2804.000000 2798.191162
    logging::TTrace<8>("lsst.ip.diffim.PsfMatchingFunctor.apply",
                       "Image range : %d %d -> %d %d : %f %f",
                       startCol, startRow, endCol, endRow, 
                       0 + *imageToConvolveLocator, 0 + *imageToNotConvolveLocator);

    std::pair<int, int> rowStep = std::make_pair(static_cast<int>(-(endCol-startCol)), 1);
    for (int row = startRow; row < endRow; ++row) {
        for (int col = startCol; col < endCol; ++col) {
            ImageT const ncImage          = *imageToNotConvolveLocator;
            double const iVariance        = 1.0 / *varianceLocator;
            
            // kernel index i
            typename std::vector<typename image::Image<ImageT>::xy_locator>::iterator citeri = convolvedLocatorList.begin();
            typename std::vector<typename image::Image<ImageT>::xy_locator>::iterator citerE = convolvedLocatorList.end();
            for (int kidxi = 0; citeri != citerE; ++citeri, ++kidxi) {
                ImageT const cdImagei = **citeri;
                
                // kernel index j
                typename std::vector<typename image::Image<ImageT>::xy_locator>::iterator citerj = citeri;
                for (int kidxj = kidxi; citerj != citerE; ++citerj, ++kidxj) {
                    ImageT const cdImagej  = **citerj;
                    M(kidxi, kidxj) += cdImagei*cdImagej*iVariance;
                } 
                
                B(kidxi) += ncImage*cdImagei*iVariance;
                
                // Constant background term; effectively j = kidxj + 1 */
                M(kidxi, nParameters-1) += cdImagei*iVariance;
            } 
            
            // Background term; effectively i = kidxi + 1 
            B(nParameters-1)                += ncImage*iVariance;
            M(nParameters-1, nParameters-1) += 1.0*iVariance;
            
            // Step each accessor in column
            ++imageToConvolveLocator.x();
            ++imageToNotConvolveLocator.x();
            ++varianceLocator.x();
            for (int ki = 0; ki < nKernelParameters; ++ki) {
                ++convolvedLocatorList[ki].x();
            }             
            
        } // col
        
        // Get to next row, first col
        imageToConvolveLocator    += rowStep;
        imageToNotConvolveLocator += rowStep;
        varianceLocator += rowStep;

        for (int ki = 0; ki < nKernelParameters; ++ki) {
            convolvedLocatorList[ki] += rowStep;
        }
        
    } // row
    
    /** @note If we are going to regularize the solution to M, this is the place
     * to do it 
     */
    
    // Fill in rest of M
    for (int kidxi=0; kidxi < nParameters; ++kidxi) {
        for (int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) {
            M(kidxj, kidxi) = M(kidxi, kidxj);
        }
    }
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                       "Total compute time to step through pixels : %.2f s", time);
    t.restart();

    //std::cout << "B eigen : " << B << std::endl;

    // To use Cholesky decomposition, the matrix needs to be symmetric (M is, by
    // design) and positive definite.  
    //
    // Eventually put a check in here to make sure its positive definite
    //
    Eigen::VectorXd Soln = Eigen::VectorXd::Zero(nParameters);;
    if (!( M.ldlt().solve(B, &Soln) )) {
        logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                           "Unable to determine kernel via Cholesky LDL^T");
        if (!( M.llt().solve(B, &Soln) )) {
            logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                               "Unable to determine kernel via Cholesky LL^T");
            if (!( M.lu().solve(B, &Soln) )) {
                logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                                   "Unable to determine kernel via LU");
                // LAST RESORT
                try {
                    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eVecValues(M);
                    Eigen::MatrixXd const& R = eVecValues.eigenvectors();
                    Eigen::VectorXd eValues  = eVecValues.eigenvalues();
                    
                    for (int i = 0; i != eValues.rows(); ++i) {
                        if (eValues(i) != 0.0) {
                            eValues(i) = 1.0/eValues(i);
                        }
                    }
                    
                    Soln = R*eValues.asDiagonal()*R.transpose()*B;
                } catch (exceptions::Exception& e) {
                    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                                       "Unable to determine kernel via eigen-values");
                    
                    throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel solution in PsfMatchingFunctor::apply");
                }
            }
        }
    }
    //std::cout << "Soln eigen : " << Soln << std::endl;
    //return;

    // Estimate of parameter uncertainties comes from the inverse of the
    // covariance matrix (noise spectrum).  
    // N.R. 15.4.8 to 15.4.15
    // 
    // Since this is a linear problem no need to use Fisher matrix
    // N.R. 15.5.8

    // Although I might be able to take advantage of the solution above.
    // Since this now works and is not the rate limiting step, keep as-is for DC3a.

    // Use Cholesky decomposition again.
    // Cholkesy:
    // Cov       =  L L^t
    // Cov^(-1)  = (L L^t)^(-1)
    //           = (L^T)^-1 L^(-1)
    Eigen::MatrixXd             Cov    = M.transpose() * M;
    Eigen::LLT<Eigen::MatrixXd> llt    = Cov.llt();
    Eigen::MatrixXd             Error2 = llt.matrixL().transpose().inverse() * llt.matrixL().inverse();
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                       "Total compute time to do matrix math : %.2f s", time);
    
    // Translate from Eigen vectors into LSST classes
    int const kCols = policy.getInt("kernelCols");
    int const kRows = policy.getInt("kernelRows");
    std::vector<double> kValues(kCols*kRows);
    std::vector<double> kErrValues(kCols*kRows);
    for (int row = 0, idx = 0; row < kRows; row++) {
        for (int col = 0; col < kCols; col++, idx++) {
            
            // Insanity checking
            if (std::isnan( Soln(idx) )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel solution (nan)");
            }
            if (std::isnan( Error2(idx, idx) )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
            }
            if (Error2(idx, idx) < 0.0) {
                throw LSST_EXCEPT(exceptions::Exception,
                                  str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                                      Error2(idx, idx)
                                      ));
            }
            
            kValues[idx]    = Soln(idx);
            kErrValues[idx] = sqrt(Error2(idx, idx));
        }
    }
    this->_kernel.reset( new math::LinearCombinationKernel(this->_basisList, kValues) );
    this->_kernelError.reset( new math::LinearCombinationKernel(this->_basisList, kErrValues) );
    
    // Estimate of Background and Background Error */
    if (std::isnan( Error2(nParameters-1, nParameters-1) )) {
        throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
    }
    if (Error2(nParameters-1, nParameters-1) < 0.0) {
        throw LSST_EXCEPT(exceptions::Exception, 
                          str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                              Error2(nParameters-1, nParameters-1) 
                              ));
    }
    this->_background      = Soln(nParameters-1);
    this->_backgroundError = sqrt(Error2(nParameters-1, nParameters-1));
}


template <typename ImageT, typename VarT>
void diffim::PsfMatchingFunctorGsl<ImageT, VarT>::apply(
    image::MaskedImage<ImageT> const& imageToConvolve,
    image::MaskedImage<ImageT> const& imageToNotConvolve,
    image::Image<VarT>         const& varianceEstimate,
    lsst::pex::policy::Policy             const& policy
    ) {
    
    // Make sure you do not overwrite anyone else's kernels
    this->reset();

    // grab mask bits from the image to convolve, since that is what we'll be operating on
    int edgeMaskBit = imageToConvolve.getMask()->getMaskPlane("EDGE");
    
    int nKernelParameters     = this->_basisList.size();
    int nBackgroundParameters = 1;
    int nParameters           = nKernelParameters + nBackgroundParameters;
    
    boost::timer t;
    double time;
    t.restart();
    
    gsl_vector *B = gsl_vector_alloc (nParameters);
    gsl_matrix *M = gsl_matrix_alloc (nParameters, nParameters);
    gsl_vector_set_zero(B);
    gsl_matrix_set_zero(M);
    
    std::vector<boost::shared_ptr<image::MaskedImage<ImageT> > > convolvedImageList(nKernelParameters);
    typename std::vector<boost::shared_ptr<image::MaskedImage<ImageT> > >::iterator 
        citer = convolvedImageList.begin();
    std::vector<boost::shared_ptr<math::Kernel> >::const_iterator 
        kiter = this->_basisList.begin();
    
    // Create C_ij in the formalism of Alard & Lupton */
    for (; kiter != this->_basisList.end(); ++kiter, ++citer) {
        
        /* NOTE : we could also *precompute* the entire template image convolved with these functions */
        /*        and save them somewhere to avoid this step each time.  however, our paradigm is to */
        /*        compute whatever is needed on the fly.  hence this step here. */
        image::MaskedImage<ImageT> image(imageToConvolve.getDimensions());
        math::convolve(image,
                       imageToConvolve,
                       **kiter,
                       false,
                       edgeMaskBit);
        boost::shared_ptr<image::MaskedImage<ImageT> > imagePtr( new image::MaskedImage<ImageT>(image) );
        *citer = imagePtr;
    } 

    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctorGsl.apply", 
                       "Total compute time to do basis convolutions : %.2f s", time);
    t.restart();
    
    kiter = this->_basisList.begin();
    citer = convolvedImageList.begin();

    // Ignore buffers around edge of convolved images :
    //
    // If the kernel has width 5, it has center pixel 2.  The first good pixel
    // is the (5-2)=3rd pixel, which is array index 2, and ends up being the
    // index of the central pixel.
    //
    // You also have a buffer of unusable pixels on the other side, numbered
    // width-center-1.  The last good usable pixel is N-width+center+1.

    // Example : the kernel is width = 5, center = 2
    //
    //     ---|---|-c-|---|---|
    //          
    //           the image is width = N
    //           convolve this with the kernel, and you get
    //
    //    |-x-|-x-|-g-|---|---| ... |---|---|-g-|-x-|-x-|
    //
    //           g = first/last good pixel
    //           x = bad
    // 
    //           the first good pixel is the array index that has the value "center", 2
    //           the last good pixel has array index N-(5-2)+1
    //           eg. if N = 100, you want to use up to index 97
    //               100-3+1 = 98, and the loops use i < 98, meaning the last
    //               index you address is 97.
   
    unsigned int startCol = (*kiter)->getCtrX();
    unsigned int startRow = (*kiter)->getCtrY();
    unsigned int endCol   = (*citer)->getWidth()  - ((*kiter)->getWidth()  - (*kiter)->getCtrX()) + 1;
    unsigned int endRow   = (*citer)->getHeight() - ((*kiter)->getHeight() - (*kiter)->getCtrY()) + 1;

    std::vector<xy_locator> convolvedLocatorList;
    for (citer = convolvedImageList.begin(); citer != convolvedImageList.end(); ++citer) {
        convolvedLocatorList.push_back( (**citer).xy_at(startCol,startRow) );
    }
    xy_locator  imageToConvolveLocator    = imageToConvolve.xy_at(startCol, startRow);
    xy_locator  imageToNotConvolveLocator = imageToNotConvolve.xy_at(startCol, startRow);
    xyi_locator varianceLocator           = varianceEstimate.xy_at(startCol, startRow);

    // Unit test ImageSubtract_1.py should show
    // Image range : 9 9 -> 31 31 : 2804.000000 2798.191162
    logging::TTrace<8>("lsst.ip.diffim.PsfMatchingFunctor.applyGsl",
                       "Image range : %d %d -> %d %d : %f %f",
                       startCol, startRow, endCol, endRow, 
                       imageToConvolveLocator.image(), 
                       imageToNotConvolveLocator.image());

    std::pair<int, int> rowStep = std::make_pair(static_cast<int>(-(endCol-startCol)), 1);
    for (unsigned int row = startRow; row < endRow; ++row) {
        
        for (unsigned int col = startCol; col < endCol; ++col) {
            
            ImageT const ncImage          = imageToNotConvolveLocator.image();
            ImageT const ncVariance       = imageToNotConvolveLocator.variance();
            image::MaskPixel const ncMask = imageToNotConvolveLocator.mask();
            double const iVariance        = 1.0 / *varianceLocator;
            
            // Unit test ImageSubtract_1.py should show
            // Accessing image row 9 col 9  : 2798.191 23.426 0 1792.511475
            // Accessing image row 9 col 10 : 2805.171 23.459 0 1774.878662
            // ...
            // Accessing image row 9 col 30 : 2793.281 23.359 0 1779.194946
            // Accessing image row 10 col 9 : 2802.968 23.464 0 1770.467163
            // ...
            logging::TTrace<8>("lsst.ip.diffim.PsfMatchingFunctor.applyGsl",
                               "Accessing image row %d col %d : %.3f %.3f %d %f",
                               row, col, ncImage, ncVariance, ncMask, 1.0 * *varianceLocator);
            
            // kernel index i
            typename std::vector<xy_locator>::iterator 
                citeri = convolvedLocatorList.begin();
            typename std::vector<xy_locator>::iterator 
                citerE = convolvedLocatorList.end();
            for (int kidxi = 0; citeri != citerE; ++citeri, ++kidxi) {
                ImageT           const cdImagei = (*citeri).image();
                image::MaskPixel const cdMaski  = (*citeri).mask();
                if (cdMaski != 0) {
                    throw LSST_EXCEPT(exceptions::Exception, 
                                      str(boost::format("Accessing invalid pixel (%d) in PsfMatchingFunctor::applyGsl") % 
                                          kidxi));
                }                
                
                // kernel index j
                typename std::vector<xy_locator>::iterator 
                    citerj = citeri;
                for (int kidxj = kidxi; citerj != citerE; ++citerj, ++kidxj) {
                    ImageT const cdImagej  = (*citerj).image();
                    *gsl_matrix_ptr(M, kidxi, kidxj) += cdImagei * cdImagej * iVariance;
                } 
                
                *gsl_vector_ptr(B, kidxi) += ncImage * cdImagei * iVariance;
                
                // Constant background term; effectively j=kidxj+1 */
                *gsl_matrix_ptr(M, kidxi, nParameters-1) += cdImagei * iVariance;
            } 
            
            // Background term; effectively i=kidxi+1 
            *gsl_vector_ptr(B, nParameters-1)                += ncImage * iVariance;
            *gsl_matrix_ptr(M, nParameters-1, nParameters-1) += 1.0 * iVariance;
            
            // Step each accessor in column
            ++imageToConvolveLocator.x();
            ++imageToNotConvolveLocator.x();
            ++varianceLocator.x();
            for (int ki = 0; ki < nKernelParameters; ++ki) {
                ++convolvedLocatorList[ki].x();
            }             
            
        } // col
        
        // Get to next row, first col
        imageToConvolveLocator    += rowStep;
        imageToNotConvolveLocator += rowStep;

        // HACK UNTIL Ticket #647 FIXED
        varianceLocator            = varianceEstimate.xy_at(startCol, row+1);

        for (int ki = 0; ki < nKernelParameters; ++ki) {
            convolvedLocatorList[ki] += rowStep;
        }
        
    } // row

    
    /** @note If we are going to regularize the solution to M, this is the place
     * to do it 
     */
    
    // Fill in rest of M
    for (int kidxi=0; kidxi < nParameters; ++kidxi) 
        for (int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) 
            gsl_matrix_set(M, kidxj, kidxi,
                           gsl_matrix_get(M, kidxi, kidxj));
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctorGsl.apply", 
                       "Total compute time to step through pixels : %.2f s", time);
    t.restart();

    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (nParameters, nParameters);
    gsl_vector *Soln                    = gsl_vector_alloc (nParameters);
    gsl_matrix *Error2                  = gsl_matrix_alloc (nParameters, nParameters);
    double chi2;
    size_t rank;
    gsl_multifit_linear_svd(M, B, GSL_DBL_EPSILON, &rank, Soln, Error2, &chi2, work);

    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctorGsl.apply", 
                       "Total compute time to do matrix math : %.2f s", time);
    
    // Translate from Gsl vectors into LSST classes
    unsigned int kCols = policy.getInt("kernelCols");
    unsigned int kRows = policy.getInt("kernelRows");
    std::vector<double> kValues(kCols*kRows);
    std::vector<double> kErrValues(kCols*kRows);
    for (unsigned int row = 0, idx = 0; row < kRows; row++) {
        for (unsigned int col = 0; col < kCols; col++, idx++) {
            
            // Insanity checking
            if (std::isnan( gsl_vector_get(Soln,idx) )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel solution (nan)");
            }
            if (std::isnan( gsl_matrix_get(Error2, idx, idx) )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
            }
            if (gsl_matrix_get(Error2, idx, idx) < 0.0) {
                throw LSST_EXCEPT(exceptions::Exception,
                                  str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                                      gsl_matrix_get(Error2, idx, idx)
                                      ));
            }

            kValues[idx]    = gsl_vector_get(Soln, idx);
            kErrValues[idx] = sqrt(gsl_matrix_get(Error2, idx, idx));
        }
    }
    this->_kernel.reset( new math::LinearCombinationKernel(this->_basisList, kValues) );
    this->_kernelError.reset( new math::LinearCombinationKernel(this->_basisList, kErrValues) );
    
    // Estimate of Background and Background Error */
    if (std::isnan( gsl_matrix_get(Error2, nParameters-1, nParameters-1) )) {
        throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
    }
    if (gsl_matrix_get(Error2, nParameters-1, nParameters-1) < 0.0) {
        throw LSST_EXCEPT(exceptions::Exception, 
                          str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                              gsl_matrix_get(Error2, nParameters-1, nParameters-1) 
                              ));
    }
    this->_background      = gsl_vector_get(Soln, nParameters-1);
    this->_backgroundError = sqrt(gsl_matrix_get(Error2, nParameters-1, nParameters-1));
}

#if USE_VW
template <typename ImageT, typename VarT>
void diffim::PsfMatchingFunctorVw<ImageT, VarT>::apply(
    image::MaskedImage<ImageT> const& imageToConvolve,
    image::MaskedImage<ImageT> const& imageToNotConvolve,
    image::Image<VarT>         const& varianceEstimate,
    lsst::pex::policy::Policy             const& policy
    ) {
    
    // Make sure you do not overwrite anyone else's kernels
    this->reset();

    // grab mask bits from the image to convolve, since that is what we'll be operating on
    int edgeMaskBit = imageToConvolve.getMask()->getMaskPlane("EDGE");
    
    int nKernelParameters     = this->_basisList.size();
    int nBackgroundParameters = 1;
    int nParameters           = nKernelParameters + nBackgroundParameters;
    
    boost::timer t;
    double time;
    t.restart();
    
    vw::math::Vector<double> B(nParameters);
    vw::math::Matrix<double> M(nParameters, nParameters);
    for (unsigned int i = nParameters; i--;) {
        B(i) = 0;
        for (unsigned int j = nParameters; j--;) {
            M(i,j) = 0;
        }
    }
    
    std::vector<boost::shared_ptr<image::MaskedImage<ImageT> > > convolvedImageList(nKernelParameters);
    typename std::vector<boost::shared_ptr<image::MaskedImage<ImageT> > >::iterator 
        citer = convolvedImageList.begin();
    std::vector<boost::shared_ptr<math::Kernel> >::const_iterator 
        kiter = this->_basisList.begin();
    
    // Create C_ij in the formalism of Alard & Lupton */
    for (; kiter != this->_basisList.end(); ++kiter, ++citer) {
        
        /* NOTE : we could also *precompute* the entire template image convolved with these functions */
        /*        and save them somewhere to avoid this step each time.  however, our paradigm is to */
        /*        compute whatever is needed on the fly.  hence this step here. */
        image::MaskedImage<ImageT> image(imageToConvolve.getDimensions());
        math::convolve(image,
                       imageToConvolve,
                       **kiter,
                       false,
                       edgeMaskBit);
        boost::shared_ptr<image::MaskedImage<ImageT> > imagePtr( new image::MaskedImage<ImageT>(image) );
        *citer = imagePtr;
    } 
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctorVw.apply", 
                       "Total compute time to do basis convolutions : %.2f s", time);
    t.restart();

    kiter = this->_basisList.begin();
    citer = convolvedImageList.begin();

    // Ignore buffers around edge of convolved images :
    //
    // If the kernel has width 5, it has center pixel 2.  The first good pixel
    // is the (5-2)=3rd pixel, which is array index 2, and ends up being the
    // index of the central pixel.
    //
    // You also have a buffer of unusable pixels on the other side, numbered
    // width-center-1.  The last good usable pixel is N-width+center+1.

    // Example : the kernel is width = 5, center = 2
    //
    //     ---|---|-c-|---|---|
    //          
    //           the image is width = N
    //           convolve this with the kernel, and you get
    //
    //    |-x-|-x-|-g-|---|---| ... |---|---|-g-|-x-|-x-|
    //
    //           g = first/last good pixel
    //           x = bad
    // 
    //           the first good pixel is the array index that has the value "center", 2
    //           the last good pixel has array index N-(5-2)+1
    //           eg. if N = 100, you want to use up to index 97
    //               100-3+1 = 98, and the loops use i < 98, meaning the last
    //               index you address is 97.
   
    unsigned int startCol = (*kiter)->getCtrX();
    unsigned int startRow = (*kiter)->getCtrY();
    unsigned int endCol   = (*citer)->getWidth()  - ((*kiter)->getWidth()  - (*kiter)->getCtrX()) + 1;
    unsigned int endRow   = (*citer)->getHeight() - ((*kiter)->getHeight() - (*kiter)->getCtrY()) + 1;

    std::vector<xy_locator> convolvedLocatorList;
    for (citer = convolvedImageList.begin(); citer != convolvedImageList.end(); ++citer) {
        convolvedLocatorList.push_back( (**citer).xy_at(startCol,startRow) );
    }
    xy_locator  imageToConvolveLocator    = imageToConvolve.xy_at(startCol, startRow);
    xy_locator  imageToNotConvolveLocator = imageToNotConvolve.xy_at(startCol, startRow);
    xyi_locator varianceLocator           = varianceEstimate.xy_at(startCol, startRow);

    // Unit test ImageSubtract_1.py should show
    // Image range : 9 9 -> 31 31 : 2804.000000 2798.191162
    logging::TTrace<8>("lsst.ip.diffim.PsfMatchingFunctor.applyVW",
                       "Image range : %d %d -> %d %d : %f %f",
                       startCol, startRow, endCol, endRow, 
                       imageToConvolveLocator.image(), 
                       imageToNotConvolveLocator.image());

    std::pair<int, int> rowStep = std::make_pair(static_cast<int>(-(endCol-startCol)), 1);
    for (unsigned int row = startRow; row < endRow; ++row) {
        
        for (unsigned int col = startCol; col < endCol; ++col) {
            
            ImageT const ncImage          = imageToNotConvolveLocator.image();
            ImageT const ncVariance       = imageToNotConvolveLocator.variance();
            image::MaskPixel const ncMask = imageToNotConvolveLocator.mask();
            double const iVariance        = 1.0 / *varianceLocator;
            
            // Unit test ImageSubtract_1.py should show
            // Accessing image row 9 col 9  : 2798.191 23.426 0 1792.511475
            // Accessing image row 9 col 10 : 2805.171 23.459 0 1774.878662
            // ...
            // Accessing image row 9 col 30 : 2793.281 23.359 0 1779.194946
            // Accessing image row 10 col 9 : 2802.968 23.464 0 1770.467163
            // ...
            logging::TTrace<8>("lsst.ip.diffim.PsfMatchingFunctor.applyVW",
                               "Accessing image row %d col %d : %.3f %.3f %d %f",
                               row, col, ncImage, ncVariance, ncMask, 1.0 * *varianceLocator);
            
            // kernel index i
            typename std::vector<xy_locator>::iterator 
                citeri = convolvedLocatorList.begin();
            typename std::vector<xy_locator>::iterator 
                citerE = convolvedLocatorList.end();
            for (int kidxi = 0; citeri != citerE; ++citeri, ++kidxi) {
                ImageT           const cdImagei = (*citeri).image();
                image::MaskPixel const cdMaski  = (*citeri).mask();
                if (cdMaski != 0) {
                    throw LSST_EXCEPT(exceptions::Exception, 
                                      str(boost::format("Accessing invalid pixel (%d) in PsfMatchingFunctor::applyVW") % 
                                          kidxi));
                }                
                
                // kernel index j
                typename std::vector<xy_locator>::iterator 
                    citerj = citeri;
                for (int kidxj = kidxi; citerj != citerE; ++citerj, ++kidxj) {
                    ImageT const cdImagej  = (*citerj).image();
                    M[kidxi][kidxj]   += cdImagei * cdImagej * iVariance;
                } 
                
                B[kidxi] += ncImage * cdImagei * iVariance;
                
                // Constant background term; effectively j=kidxj+1 */
                M[kidxi][nParameters-1] += cdImagei * iVariance;
            } 
            
            // Background term; effectively i=kidxi+1 
            B[nParameters-1]                += ncImage * iVariance;
            M[nParameters-1][nParameters-1] += 1.0 * iVariance;
            
            // Step each accessor in column
            ++imageToConvolveLocator.x();
            ++imageToNotConvolveLocator.x();
            ++varianceLocator.x();
            for (int ki = 0; ki < nKernelParameters; ++ki) {
                ++convolvedLocatorList[ki].x();
            }             
            
        } // col
        
        // Get to next row, first col
        imageToConvolveLocator    += rowStep;
        imageToNotConvolveLocator += rowStep;

        // HACK UNTIL Ticket #647 FIXED
        varianceLocator            = varianceEstimate.xy_at(startCol, row+1);

        for (int ki = 0; ki < nKernelParameters; ++ki) {
            convolvedLocatorList[ki] += rowStep;
        }
        
    } // row

    
    /** @note If we are going to regularize the solution to M, this is the place
     * to do it 
     */
    
    // Fill in rest of M
    for (int kidxi=0; kidxi < nParameters; ++kidxi) 
        for (int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) 
            M[kidxj][kidxi] = M[kidxi][kidxj];
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctorVw.apply", 
                       "Total compute time to step through pixels : %.2f s", time);
    t.restart();

    // Invert using VW's internal method
    vw::math::Vector<double> Soln      = vw::math::least_squares(M, B);

    // Additional gymnastics to get the parameter uncertainties
    vw::math::Matrix<double> Mt        = vw::math::transpose(M);
    vw::math::Matrix<double> MtM       = Mt * M;
    vw::math::Matrix<double> Error2    = vw::math::pseudoinverse(MtM);

    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctorVw.apply", 
                       "Total compute time to do matrix math : %.2f s", time);
    
    // Translate from VW vectors into LSST classes
    unsigned int kCols = policy.getInt("kernelCols");
    unsigned int kRows = policy.getInt("kernelRows");
    std::vector<double> kValues(kCols*kRows);
    std::vector<double> kErrValues(kCols*kRows);
    for (unsigned int row = 0, idx = 0; row < kRows; row++) {
        for (unsigned int col = 0; col < kCols; col++, idx++) {
            
            // Insanity checking
            if (std::isnan( Soln[idx] )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel solution (nan)");
            }
            if (std::isnan( Error2[idx][idx] )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
            }
            if (Error2[idx][idx] < 0.0) {
                throw LSST_EXCEPT(exceptions::Exception,
                                  str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                                      Error2[idx][idx]
                                      ));
            }
            
            kValues[idx]    = Soln[idx];
            kErrValues[idx] = sqrt(Error2[idx][idx]);
        }
    }
    this->_kernel.reset( new math::LinearCombinationKernel(this->_basisList, kValues) );
    this->_kernelError.reset( new math::LinearCombinationKernel(this->_basisList, kErrValues) );
    
    // Estimate of Background and Background Error */
    if (std::isnan( Error2[nParameters-1][nParameters-1] )) {
        throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
    }
    if (Error2[nParameters-1][nParameters-1] < 0.0) {
        throw LSST_EXCEPT(exceptions::Exception, 
                          str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                              Error2[nParameters-1][nParameters-1] 
                              ));
    }
    this->_background      = Soln[nParameters-1];
    this->_backgroundError = sqrt(Error2[nParameters-1][nParameters-1]);
}
#endif

//
// Subroutines
//

/** 
 * @brief Generate a basis set of delta function Kernels.
 *
 * Generates a vector of Kernels sized nCols * nRows, where each Kernel has
 * a unique pixel set to value 1.0 with the other pixels valued 0.0.  This
 * is the "delta function" basis set.
 * 
 * @return Vector of orthonormal delta function Kernels.
 *
 * @throw lsst::pex::exceptions::DomainError if nRows or nCols not positive
 *
 * @ingroup diffim
 */
math::KernelList<math::Kernel>
diffim::generateDeltaFunctionKernelSet(
    unsigned int width,                 ///< number of columns in the set
    unsigned int height                 ///< number of rows in the set
    ) {
    if ((width < 1) || (height < 1)) {
        throw LSST_EXCEPT(exceptions::Exception, "nRows and nCols must be positive");
    }
    const int signedWidth = static_cast<int>(width);
    const int signedHeight = static_cast<int>(height);
    math::KernelList<math::Kernel> kernelBasisList;
    for (int row = 0; row < signedHeight; ++row) {
        for (int col = 0; col < signedWidth; ++col) {
            boost::shared_ptr<math::Kernel> 
                kernelPtr( new math::DeltaFunctionKernel(width, height, image::PointI(col,row) ) );
            kernelBasisList.push_back(kernelPtr);
        }
    }
    return kernelBasisList;
}

/** 
 * @brief Generate an Alard-Lupton basis set of Kernels.
 *
 * Not implemented.
 * 
 * @return Vector of Alard-Lupton Kernels.
 *
 * @throw lsst::pex::exceptions::DomainError if nRows or nCols not positive
 * @throw lsst::pex::exceptions::DomainError until implemented
 *
 * @ingroup diffim
 */
math::KernelList<math::Kernel>
diffim::generateAlardLuptonKernelSet(
    unsigned int nRows, 
    unsigned int nCols, 
    std::vector<double> const& sigGauss, 
    std::vector<double> const& degGauss  
    ) {
    if ((nCols < 1) || (nRows < 1)) {
        throw LSST_EXCEPT(exceptions::Exception, "nRows and nCols must be positive");
    }
    throw LSST_EXCEPT(exceptions::Exception, "Not implemented");
    
    math::KernelList<math::Kernel> kernelBasisList;
    return kernelBasisList;
}

/************************************************************************************************************/
/*
 * Adds a Function to an Image
 *
 * @note MAJOR NOTE; I need to check if my scaling of the image range from -1 to
 * 1 gets messed up here.  ACB.
 *
 * @note This routine assumes that the pixel coordinates start at (0, 0) which is
 * in general not true
 *
 * @node this function was renamed from addFunctionToImage to addSomethingToImage to allow generic programming
 */
namespace {
    template <typename ImageT, typename FunctionT>
    void addSomethingToImage(ImageT &image,
                             FunctionT const& function
                            ) {

        // Set the pixels row by row, to avoid repeated checks for end-of-row
        for (int y = 0; y != image.getHeight(); ++y) {
            double yPos = image::positionToIndex(y);
        
            double xPos = image::positionToIndex(0);
            for (typename ImageT::x_iterator ptr = image.row_begin(y), end = image.row_end(y);
                 ptr != end; ++ptr, ++xPos) {            
                *ptr += function(xPos, yPos);
            }
        }
    }
    //
    // Add a scalar.
    //
    template <typename ImageT>
    void addSomethingToImage(image::Image<ImageT> &image,
                             double value
                            ) {
        if (value != 0.0) {
            image += value;
        }
    }
}

/** 
 * @brief Implement fundamental difference imaging step of convolution and
 * subtraction : D = I - (K*T + bg) where * denotes convolution
 * 
 * @note If you convolve the science image, D = (K*I + bg) - T, set invert=False
 *
 * @note The template is taken to be an MaskedImage; this takes c 1.6 times as long
 * as using an Image
 *
 * @return Difference image
 *
 * @ingroup diffim
 */
template <typename ImageT, typename BackgroundT>
image::MaskedImage<ImageT> diffim::convolveAndSubtract(
    image::MaskedImage<ImageT> const& imageToConvolve,    ///< Image T to convolve with Kernel
    image::MaskedImage<ImageT> const& imageToNotConvolve, ///< Image I to subtract convolved template from
    math::Kernel const& convolutionKernel,                ///< PSF-matching Kernel used for convolution
    BackgroundT background,                               ///< Differential background function or scalar
    bool invert                                           ///< Invert the output difference image
    ) {
    
    logging::TTrace<8>("lsst.ip.diffim.convolveAndSubtract", "Convolving using convolve");
    
    int edgeMaskBit = imageToNotConvolve.getMask()->getMaskPlane("EDGE");
    image::MaskedImage<ImageT> convolvedMaskedImage(imageToConvolve.getDimensions());
    convolvedMaskedImage.setXY0(imageToConvolve.getXY0());
    
    math::convolve(convolvedMaskedImage, imageToConvolve, convolutionKernel, false, edgeMaskBit);
    
    /* Add in background */
    addSomethingToImage(*(convolvedMaskedImage.getImage()), background);
    
    /* Do actual subtraction */
    convolvedMaskedImage -= imageToNotConvolve;

    /* Invert */
    if (invert) {
        convolvedMaskedImage *= -1.0;
    }
    
    return convolvedMaskedImage;
}

/** 
 * @brief Implement fundamental difference imaging step of convolution and
 * subtraction : D = I - (K.x.T + bg)
 *
 * @note The template is taken to be an Image, not a MaskedImage; it therefore
 * has neither variance nor bad pixels
 *
 * @note If you convolve the science image, D = (K*I + bg) - T, set invert=False
 * 
 * @return Difference image
 *
 * @ingroup diffim
 */
template <typename ImageT, typename BackgroundT>
image::MaskedImage<ImageT> diffim::convolveAndSubtract(
    image::Image<ImageT> const& imageToConvolve,          ///< Image T to convolve with Kernel
    image::MaskedImage<ImageT> const& imageToNotConvolve, ///< Image I to subtract convolved template from
    math::Kernel const& convolutionKernel,                ///< PSF-matching Kernel used for convolution
    BackgroundT background,                               ///< Differential background function or scalar
    bool invert                                           ///< Invert the output difference image
    ) {
    
    logging::TTrace<8>("lsst.ip.diffim.convolveAndSubtract", "Convolving using convolve");
    
    int edgeMaskBit = imageToNotConvolve.getMask()->getMaskPlane("EDGE");
    image::MaskedImage<ImageT> convolvedMaskedImage(imageToConvolve.getDimensions());
    convolvedMaskedImage.setXY0(imageToConvolve.getXY0());
    
    math::convolve(*convolvedMaskedImage.getImage(), imageToConvolve, convolutionKernel, false, edgeMaskBit);
    
    /* Add in background */
    addSomethingToImage(*convolvedMaskedImage.getImage(), background);
    
    /* Do actual subtraction */
    *convolvedMaskedImage.getImage() -= *imageToNotConvolve.getImage();

    /* Invert */
    if (invert) {
        *convolvedMaskedImage.getImage() *= -1.0;
    }
    *convolvedMaskedImage.getMask() <<= *imageToNotConvolve.getMask();
    *convolvedMaskedImage.getVariance() <<= *imageToNotConvolve.getVariance();
    
    return convolvedMaskedImage;
}

/** 
 * @brief Runs Detection on a single image for significant peaks, and checks
 * returned Footprints for Masked pixels.
 *
 * Accepts two MaskedImages, one of which is to be convolved to match the
 * other.  The Detection package is run on the image to be convolved
 * (assumed to be higher S/N than the other image).  The subimages
 * associated with each returned Footprint in both images are checked for
 * Masked pixels; Footprints containing Masked pixels are rejected.  The
 * Footprints are grown by an amount specified in the Policy.  The
 * acceptible Footprints are returned in a vector.
 *
 * @return Vector of "clean" Footprints around which Image Subtraction
 * Kernels will be built.
 *
 * @ingroup diffim
 */
template <typename ImageT>
std::vector<lsst::afw::detection::Footprint::Ptr> diffim::getCollectionOfFootprintsForPsfMatching(
    image::MaskedImage<ImageT> const& imageToConvolve,    
    image::MaskedImage<ImageT> const& imageToNotConvolve, 
    lsst::pex::policy::Policy  const& policy                                       
    ) {
    
    // Parse the Policy
    unsigned int fpNpixMin      = policy.getInt("fpNpixMin");
    unsigned int fpNpixMax      = policy.getInt("fpNpixMax");
    unsigned int fpGrowPix      = policy.getInt("fpGrowPix");
    int minCleanFp              = policy.getInt("minCleanFp");
    double detThreshold         = policy.getDouble("detThresholdSigma");
    double detThresholdScaling  = policy.getDouble("detThresholdScaling");
    double detThresholdMin      = policy.getDouble("detThresholdMin");
    std::string detThresholdType = policy.getString("detThresholdType");

    // Grab mask bits from the image to convolve, since that is what we'll be operating on
    // Overridden now that we use the FootprintFunctor to look for any masked pixels
    // int badMaskBit = imageToConvolve.getMask()->getMaskPlane("BAD");
    // image::MaskPixel badPixelMask = (badMaskBit < 0) ? 0 : (1 << badMaskBit);
    
    // List of Footprints
    std::vector<lsst::afw::detection::Footprint::Ptr> footprintListIn;
    std::vector<lsst::afw::detection::Footprint::Ptr> footprintListOut;

    // Functors to search through the images for bad pixels within candidate footprints
    diffim::FindSetBits<image::Mask<image::MaskPixel> > itcFunctor(*(imageToConvolve.getMask())); 
    diffim::FindSetBits<image::Mask<image::MaskPixel> > itncFunctor(*(imageToNotConvolve.getMask())); 
 
    int nCleanFp = 0;
    while ( (nCleanFp < minCleanFp) and (detThreshold > detThresholdMin) ) {
        footprintListIn.clear();
        footprintListOut.clear();
        
        // Find detections
        lsst::afw::detection::Threshold threshold = 
                lsst::afw::detection::createThreshold(detThreshold, detThresholdType);
        lsst::afw::detection::DetectionSet<ImageT> detectionSet(
                imageToConvolve, 
                threshold,
                "",
                fpNpixMin);
        
        // Get the associated footprints
        footprintListIn = detectionSet.getFootprints();
        logging::TTrace<4>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                           "Found %d total footprints above threshold %.3f",
                           footprintListIn.size(), detThreshold);

        // Iterate over footprints, look for "good" ones
        nCleanFp = 0;
        for (std::vector<lsst::afw::detection::Footprint::Ptr>::iterator i = footprintListIn.begin(); i != footprintListIn.end(); ++i) {
            // footprint has too many pixels
            if (static_cast<unsigned int>((*i)->getNpix()) > fpNpixMax) {
                logging::TTrace<5>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                               "Footprint has too many pix: %d (max =%d)", 
                               (*i)->getNpix(), fpNpixMax);
                continue;
            } 
            
            logging::TTrace<8>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                               "Footprint in : %d,%d -> %d,%d",
                               (*i)->getBBox().getX0(), (*i)->getBBox().getX1(), 
                               (*i)->getBBox().getY0(), (*i)->getBBox().getY1());

            logging::TTrace<8>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                               "Grow by : %d pixels", fpGrowPix);

            // Grow the footprint
            // true = isotropic grow = slow
            // false = 'manhattan grow' = fast
            lsst::afw::detection::Footprint::Ptr fpGrow = 
                lsst::afw::detection::growFootprint(*i, fpGrowPix, false);
            
            logging::TTrace<6>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                               "Footprint out : %d,%d -> %d,%d (center %d,%d)",
                               (*fpGrow).getBBox().getX0(), (*fpGrow).getBBox().getY0(),
			       (*fpGrow).getBBox().getX1(), (*fpGrow).getBBox().getY1(),
			       int( 0.5 * ((*i)->getBBox().getX0()+(*i)->getBBox().getX1()) ),
			       int( 0.5 * ((*i)->getBBox().getY0()+(*i)->getBBox().getY1()) ) );


            // Grab a subimage; there is an exception if it's e.g. too close to the image */
            try {
                image::BBox fpBBox = (*fpGrow).getBBox();
                fpBBox.shift(-imageToConvolve.getX0(), -imageToConvolve.getY0());
                
                image::MaskedImage<ImageT> subImageToConvolve(imageToConvolve, fpBBox);
                image::MaskedImage<ImageT> subImageToNotConvolve(imageToNotConvolve, fpBBox);
            } catch (exceptions::Exception& e) {
                logging::TTrace<4>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching",
                                   "Exception caught extracting Footprint");
                logging::TTrace<5>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching",
                                   e.what());
                continue;
            }

            // Search for bad pixels within the footprint
            itcFunctor.apply(*fpGrow);
            if (itcFunctor.getBits() > 0) {
                logging::TTrace<5>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                               "Footprint has bad pix in image to convolve"); 
                continue;
            }

            itncFunctor.apply(*fpGrow);
            if (itncFunctor.getBits() > 0) {
                logging::TTrace<5>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                               "Footprint has bad pix in image not to convolve");
                continue;
            }

            // If we get this far, we have a clean footprint
            footprintListOut.push_back(fpGrow);
            nCleanFp += 1;
        }
        
        detThreshold *= detThresholdScaling;
    }
    if (footprintListOut.size() == 0) {
      throw LSST_EXCEPT(exceptions::Exception, 
			"Unable to find any footprints for Psf matching");
    }

    logging::TTrace<3>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                       "Found %d clean footprints above threshold %.3f",
                       footprintListOut.size(), detThreshold/detThresholdScaling);
    
    return footprintListOut;
}

// Explicit instantiations

template class diffim::PsfMatchingFunctor<float, float>;
template class diffim::PsfMatchingFunctor<double, float>;
template class diffim::PsfMatchingFunctorGsl<float, float>;
template class diffim::PsfMatchingFunctorGsl<double, float>;
#if USE_VW
template class diffim::PsfMatchingFunctorVw<float, float>;
template class diffim::PsfMatchingFunctorVw<double, float>;
#endif

template class diffim::FindSetBits<image::Mask<> >;

template class diffim::FindCounts<float>;
template class diffim::FindCounts<double>;

template class diffim::ImageStatistics<float>;
template class diffim::ImageStatistics<double>;

/* */

#define p_INSTANTIATE_convolveAndSubtract(TEMPLATE_IMAGE_T, TYPE)     \
    template \
    image::MaskedImage<TYPE> diffim::convolveAndSubtract( \
        image::TEMPLATE_IMAGE_T<TYPE> const& imageToConvolve, \
        image::MaskedImage<TYPE> const& imageToNotConvolve, \
        math::Kernel const& convolutionKernel, \
        double background, \
        bool invert);      \
    \
    template \
    image::MaskedImage<TYPE> diffim::convolveAndSubtract( \
        image::TEMPLATE_IMAGE_T<TYPE> const& imageToConvolve, \
        image::MaskedImage<TYPE> const& imageToNotConvolve, \
        math::Kernel const& convolutionKernel, \
        math::Function2<double> const& backgroundFunction, \
        bool invert); \

#define INSTANTIATE_convolveAndSubtract(TYPE) \
p_INSTANTIATE_convolveAndSubtract(Image, TYPE) \
p_INSTANTIATE_convolveAndSubtract(MaskedImage, TYPE)
/*
 * Here are the instantiations.
 *
 * Do we really need double diffim code?  It isn't sufficient to remove it here; you'll have to also remove at
 * least SpatialModelKernel<double> and swig instantiations thereof
 */
INSTANTIATE_convolveAndSubtract(float);
INSTANTIATE_convolveAndSubtract(double);

/* */

template
std::vector<lsst::afw::detection::Footprint::Ptr> diffim::getCollectionOfFootprintsForPsfMatching(
    image::MaskedImage<float> const& imageToConvolve,
    image::MaskedImage<float> const& imageToNotConvolve,
    lsst::pex::policy::Policy const& policy);

template
std::vector<lsst::afw::detection::Footprint::Ptr> diffim::getCollectionOfFootprintsForPsfMatching(
    image::MaskedImage<double> const& imageToConvolve,
    image::MaskedImage<double> const& imageToNotConvolve,
    lsst::pex::policy::Policy  const& policy);

#if false
/** 
 * @brief Computes a single Kernel (Model 1) around a single subimage.
 *
 * Accepts two MaskedImages, generally subimages of a larger image, one of which
 * is to be convolved to match the other.  The output Kernel is generated using
 * an input list of basis Kernels by finding the coefficients in front of each
 * basis.  This version accepts an input variance image, and uses GSL for the
 * matrices.
 *
 * @return Vector of coefficients representing the relative contribution of
 * each input basis function.
 *
 * @return Differential background offset between the two images
 *
 * @ingroup diffim
 */
template <typename ImageT, typename VarT>
void diffim::computePsfMatchingKernelForFootprint(
    double                          &background,
    double                          &backgroundError,
    boost::shared_ptr<math::Kernel> &kernelPtr,
    boost::shared_ptr<math::Kernel> &kernelErrorPtr,
    image::MaskedImage<ImageT>         const& imageToConvolve,    
    image::MaskedImage<ImageT>         const& imageToNotConvolve, 
    image::Image<VarT>                 const& varianceImage,      
    math::KernelList<math::Kernel>     const& kernelInBasisList,  
    lsst::pex::policy::Policy          const& policy
    ) { 

    typedef typename image::MaskedImage<ImageT>::xy_locator xy_locator;
    typedef typename image::Image<VarT>::xy_locator xyi_locator;

    // grab mask bits from the image to convolve, since that is what we'll be operating on
    int edgeMaskBit = imageToConvolve.getMask()->getMaskPlane("EDGE");
    
    int nKernelParameters = 0;
    int nBackgroundParameters = 0;
    int nParameters = 0;
    
    boost::timer t;
    double time;
    t.restart();
    
    nKernelParameters     = kernelInBasisList.size();
    nBackgroundParameters = 1;
    nParameters           = nKernelParameters + nBackgroundParameters;
    
    gsl_vector *B = gsl_vector_alloc (nParameters);
    gsl_matrix *M = gsl_matrix_alloc (nParameters, nParameters);
    
    gsl_vector_set_zero(B);
    gsl_matrix_set_zero(M);
    
    std::vector<boost::shared_ptr<image::MaskedImage<ImageT> > > convolvedImageList(nKernelParameters);
    typename std::vector<boost::shared_ptr<image::MaskedImage<ImageT> > >::iterator 
        citer = convolvedImageList.begin();
    std::vector<boost::shared_ptr<math::Kernel> >::const_iterator 
        kiter = kernelInBasisList.begin();
    
    // Create C_ij in the formalism of Alard & Lupton */
    for (; kiter != kernelInBasisList.end(); ++kiter, ++citer) {
        
        /* NOTE : we could also *precompute* the entire template image convolved with these functions */
        /*        and save them somewhere to avoid this step each time.  however, our paradigm is to */
        /*        compute whatever is needed on the fly.  hence this step here. */
        image::MaskedImage<ImageT> image(imageToConvolve.getDimensions());
        math::convolve(image,
                       imageToConvolve,
                       **kiter,
                       false,
                       edgeMaskBit);
        boost::shared_ptr<image::MaskedImage<ImageT> > imagePtr( new image::MaskedImage<ImageT>(image) );
        *citer = imagePtr;
    } 
    
    kiter = kernelInBasisList.begin();
    citer = convolvedImageList.begin();

    // Ignore buffers around edge of convolved images :
    //
    // If the kernel has width 5, it has center pixel 2.  The first good pixel
    // is the (5-2)=3rd pixel, which is array index 2, and ends up being the
    // index of the central pixel.
    //
    // You also have a buffer of unusable pixels on the other side, numbered
    // width-center-1.  The last good usable pixel is N-width+center+1.

    // Example : the kernel is width = 5, center = 2
    //
    //     ---|---|-c-|---|---|
    //          
    //           the image is width = N
    //           convolve this with the kernel, and you get
    //
    //    |-x-|-x-|-g-|---|---| ... |---|---|-g-|-x-|-x-|
    //
    //           g = first/last good pixel
    //           x = bad
    // 
    //           the first good pixel is the array index that has the value "center", 2
    //           the last good pixel has array index N-(5-2)+1
    //           eg. if N = 100, you want to use up to index 97
    //               100-3+1 = 98, and the loops use i < 98, meaning the last
    //               index you address is 97.
   
    unsigned int startCol = (*kiter)->getCtrX();
    unsigned int startRow = (*kiter)->getCtrY();
    unsigned int endCol   = (*citer)->getWidth()  - ((*kiter)->getWidth()  - (*kiter)->getCtrX()) + 1;
    unsigned int endRow   = (*citer)->getHeight() - ((*kiter)->getHeight() - (*kiter)->getCtrY()) + 1;

    std::vector<xy_locator> convolvedLocatorList;
    for (citer = convolvedImageList.begin(); citer != convolvedImageList.end(); ++citer) {
        convolvedLocatorList.push_back( (**citer).xy_at(startCol,startRow) );
    }
    xy_locator  imageToConvolveLocator    = imageToConvolve.xy_at(startCol, startRow);
    xy_locator  imageToNotConvolveLocator = imageToNotConvolve.xy_at(startCol, startRow);
    xyi_locator varianceLocator           = varianceImage.xy_at(startCol, startRow);

    // Unit test ImageSubtract_1.py should show
    // Image range : 9 9 -> 31 31 : 2804.000000 2798.191162
    logging::TTrace<8>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                       "Image range : %d %d -> %d %d : %f %f",
                       startCol, startRow, endCol, endRow, 
                       imageToConvolveLocator.image(), 
                       imageToNotConvolveLocator.image());

    std::pair<int, int> rowStep = std::make_pair(static_cast<int>(-(endCol-startCol)), 1);
    for (unsigned int row = startRow; row < endRow; ++row) {
        
        for (unsigned int col = startCol; col < endCol; ++col) {
            
            ImageT const ncImage          = imageToNotConvolveLocator.image();
            ImageT const ncVariance       = imageToNotConvolveLocator.variance();
            image::MaskPixel const ncMask = imageToNotConvolveLocator.mask();
            double const iVariance        = 1.0 / *varianceLocator;
            
            // Unit test ImageSubtract_1.py should show
            // Accessing image row 9 col 9  : 2798.191 23.426 0 1792.511475
            // Accessing image row 9 col 10 : 2805.171 23.459 0 1774.878662
            // ...
            // Accessing image row 9 col 30 : 2793.281 23.359 0 1779.194946
            // Accessing image row 10 col 9 : 2802.968 23.464 0 1770.467163
            // ...
            logging::TTrace<8>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                               "Accessing image row %d col %d : %.3f %.3f %d %f",
                               row, col, ncImage, ncVariance, ncMask, 1.0 * *varianceLocator);
            
            // kernel index i
            typename std::vector<xy_locator>::iterator 
                citeri = convolvedLocatorList.begin();
            typename std::vector<xy_locator>::iterator 
                citerE = convolvedLocatorList.end();

            for (int kidxi = 0; citeri != citerE; ++citeri, ++kidxi) {
                ImageT           const cdImagei = (*citeri).image();
                image::MaskPixel const cdMaski  = (*citeri).mask();
                if (cdMaski != 0) {
                    throw LSST_EXCEPT(exceptions::Exception, 
                                      str(boost::format("Accessing invalid pixel (%d) in computePsfMatchingKernelForFootprint") % 
                                          kidxi));
                }                
                
                // kernel index j
                typename std::vector<xy_locator>::iterator 
                    citerj = citeri;

                for (int kidxj = kidxi; citerj != citerE; ++citerj, ++kidxj) {
                    ImageT const cdImagej = (*citerj).image();
                    
                    *gsl_matrix_ptr(M, kidxi, kidxj) += cdImagei * cdImagej * iVariance;
                } 
                
                *gsl_vector_ptr(B, kidxi) += ncImage * cdImagei * iVariance;
                
                // Constant background term; effectively j=kidxj+1 */
                *gsl_matrix_ptr(M, kidxi, nParameters-1) += cdImagei * iVariance;
            } 
            
            // Background term; effectively i=kidxi+1 
            *gsl_vector_ptr(B, nParameters-1)                += ncImage * iVariance;
            *gsl_matrix_ptr(M, nParameters-1, nParameters-1) += 1.0 * iVariance;
            
            // Step each accessor in column
            ++imageToConvolveLocator.x();
            ++imageToNotConvolveLocator.x();
            ++varianceLocator.x();
            for (int ki = 0; ki < nKernelParameters; ++ki) {
                ++convolvedLocatorList[ki].x();
            }             
            
        } // col
        
        // Get to next row, first col
        imageToConvolveLocator    += rowStep;
        imageToNotConvolveLocator += rowStep;

        // HACK UNTIL Ticket #647 FIXED
        varianceLocator            = varianceImage.xy_at(startCol, row+1);

        for (int ki = 0; ki < nKernelParameters; ++ki) {
            convolvedLocatorList[ki] += rowStep;
        }
        
    } // row

    
    /** @note If we are going to regularize the solution to M, this is the place
     * to do it 
     */
    
    // Fill in rest of M
    for (int kidxi=0; kidxi < nParameters; ++kidxi) 
        for (int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) 
            gsl_matrix_set(M, kidxj, kidxi,
                           gsl_matrix_get(M, kidxi, kidxj));
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                       "Total compute time before matrix inversions : %.2f s", time);
    
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (nParameters, nParameters);
    gsl_vector *X                       = gsl_vector_alloc (nParameters);
    gsl_matrix *Cov                     = gsl_matrix_alloc (nParameters, nParameters);
    double chi2;
    size_t rank;
    gsl_multifit_linear_svd(M, B, GSL_DBL_EPSILON, &rank, X, Cov, &chi2, work);

    //std::cout << "Soln GSL ";
    //gsl_vector_fprintf(stdout, X, "%f");
    //return;

    /* guts of gsl_multifit_linear_svd
    gsl_matrix *A   = work->A;
    gsl_matrix *Q   = work->Q;
    gsl_matrix *QSI = work->QSI;
    gsl_vector *S   = work->S;
    gsl_vector *xt  = work->xt;
    gsl_vector *D   = work->D;
    gsl_matrix_memcpy (A, M);
    gsl_linalg_balance_columns (A, D);
    gsl_linalg_SV_decomp_mod (A, QSI, Q, S, xt);
    gsl_blas_dgemv (CblasTrans, 1.0, A, B, 0.0, xt);
    gsl_matrix_memcpy (QSI, Q);
    {
        double alpha0 = gsl_vector_get (S, 0);
        size_t p_eff = 0;
        
        const size_t p = M->size2;
        
        for (size_t j = 0; j < p; j++)
            {
                gsl_vector_view column = gsl_matrix_column (QSI, j);
                double alpha = gsl_vector_get (S, j);
                
                if (alpha <= GSL_DBL_EPSILON * alpha0) {
                    alpha = 0.0;
                } else {
                    alpha = 1.0 / alpha;
                    p_eff++;
                }
                
                gsl_vector_scale (&column.vector, alpha);
            }
        
        rank = p_eff;
    }
    gsl_vector_set_zero (X);
    gsl_blas_dgemv (CblasNoTrans, 1.0, QSI, xt, 0.0, X);
    gsl_vector_div (X, D);
    gsl_multifit_linear_free(work);
    */
    //gsl_matrix_fprintf(stdout,M,"%f");
    //gsl_vector_fprintf(stdout,B,"%f");
    //gsl_vector_fprintf(stdout,X,"%f");

    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                       "Total compute time after matrix inversions : %.2f s", time);
    
    // Translate from GSL vectors into LSST classes
    unsigned int kCols = policy.getInt("kernelCols");
    unsigned int kRows = policy.getInt("kernelRows");
    std::vector<double> kValues(kCols*kRows);
    std::vector<double> kErrValues(kCols*kRows);
    for (unsigned int row = 0, idx = 0; row < kRows; row++) {
        for (unsigned int col = 0; col < kCols; col++, idx++) {
            
            // Insanity checking
            if (std::isnan( gsl_vector_get(X,idx) )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel solution (nan)");
            }
            if (std::isnan( gsl_matrix_get(Cov, idx, idx) )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
            }
            if (gsl_matrix_get(Cov, idx, idx) < 0.0) {
                throw LSST_EXCEPT(exceptions::Exception,
                                  str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                                      gsl_matrix_get(Cov, idx, idx)
                                      ));
            }
            
            kValues[idx]    = gsl_vector_get(X, idx);
            kErrValues[idx] = sqrt(gsl_matrix_get(Cov, idx, idx));
        }
    }
    kernelPtr = boost::shared_ptr<math::Kernel> 
        (new math::LinearCombinationKernel(kernelInBasisList, kValues));
    kernelErrorPtr = boost::shared_ptr<math::Kernel> 
        (new math::LinearCombinationKernel(kernelInBasisList, kErrValues));
    //kernelPtr.reset( new math::LinearCombinationKernel(kernelInBasisList, kValues) );
    //kernelErrorPtr.reset( new math::LinearCombinationKernel(kernelInBasisList, kErrValues) );
    
    // Estimate of Background and Background Error */
    if (std::isnan( gsl_matrix_get(Cov, nParameters-1, nParameters-1) )) {
        throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
    }
    if (gsl_matrix_get(Cov, nParameters-1, nParameters-1) < 0.0) {
        throw LSST_EXCEPT(exceptions::Exception, 
                          str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                              gsl_matrix_get(Cov, nParameters-1, nParameters-1) 
                              ));
    }
    background      = gsl_vector_get(X, nParameters-1);
    backgroundError = sqrt(gsl_matrix_get(Cov, nParameters-1, nParameters-1));
    
}


/** 
 * @brief Computes a single Kernel (Model 1) around a single subimage.
 *
 * Accepts two MaskedImages, generally subimages of a larger image, one of which
 * is to be convolved to match the other.  The output Kernel is generated using
 * an input list of basis Kernels by finding the coefficients in front of each
 * basis.  This version accepts an input variance image, and uses Eigen for the
 * matrices.
 *
 * @return Vector of coefficients representing the relative contribution of
 * each input basis function.
 *
 * @return Differential background offset between the two images
 *
 * @ingroup diffim
 */
template <typename ImageT, typename VarT>
void diffim::computePsfMatchingKernelForFootprintEigen(
    double                          &background,
    double                          &backgroundError,
    boost::shared_ptr<math::Kernel> &kernelPtr,
    boost::shared_ptr<math::Kernel> &kernelErrorPtr,
    image::MaskedImage<ImageT>         const& imageToConvolve,    
    image::MaskedImage<ImageT>         const& imageToNotConvolve, 
    image::Image<VarT>                 const& varianceImage,      
    math::KernelList<math::Kernel>     const& kernelInBasisList,  
    lsst::pex::policy::Policy          const& policy
    ) { 

    typedef typename image::MaskedImage<ImageT>::xy_locator xy_locator;
    typedef typename image::Image<VarT>::xy_locator xyi_locator;

    // grab mask bits from the image to convolve, since that is what we'll be operating on
    int edgeMaskBit = imageToConvolve.getMask()->getMaskPlane("EDGE");
    
    int nKernelParameters = 0;
    int nBackgroundParameters = 0;
    int nParameters = 0;
    
    boost::timer t;
    double time;
    t.restart();
    
    nKernelParameters     = kernelInBasisList.size();
    nBackgroundParameters = 1;
    nParameters           = nKernelParameters + nBackgroundParameters;

    Eigen::MatrixXd M     = Eigen::MatrixXd::Zero(nParameters, nParameters);
    Eigen::VectorXd B     = Eigen::VectorXd::Zero(nParameters);
    
    std::vector<boost::shared_ptr<image::MaskedImage<ImageT> > > convolvedImageList(nKernelParameters);
    typename std::vector<boost::shared_ptr<image::MaskedImage<ImageT> > >::iterator 
        citer = convolvedImageList.begin();
    std::vector<boost::shared_ptr<math::Kernel> >::const_iterator 
        kiter = kernelInBasisList.begin();
    
    // Create C_ij in the formalism of Alard & Lupton */
    for (; kiter != kernelInBasisList.end(); ++kiter, ++citer) {
        
        /* NOTE : we could also *precompute* the entire template image convolved with these functions */
        /*        and save them somewhere to avoid this step each time.  however, our paradigm is to */
        /*        compute whatever is needed on the fly.  hence this step here. */
        image::MaskedImage<ImageT> image(imageToConvolve.getDimensions());
        math::convolve(image,
                       imageToConvolve,
                       **kiter,
                       false,
                       edgeMaskBit);
        boost::shared_ptr<image::MaskedImage<ImageT> > imagePtr( new image::MaskedImage<ImageT>(image) );
        *citer = imagePtr;
    } 
    
    kiter = kernelInBasisList.begin();
    citer = convolvedImageList.begin();

    // Ignore buffers around edge of convolved images :
    //
    // If the kernel has width 5, it has center pixel 2.  The first good pixel
    // is the (5-2)=3rd pixel, which is array index 2, and ends up being the
    // index of the central pixel.
    //
    // You also have a buffer of unusable pixels on the other side, numbered
    // width-center-1.  The last good usable pixel is N-width+center+1.

    // Example : the kernel is width = 5, center = 2
    //
    //     ---|---|-c-|---|---|
    //          
    //           the image is width = N
    //           convolve this with the kernel, and you get
    //
    //    |-x-|-x-|-g-|---|---| ... |---|---|-g-|-x-|-x-|
    //
    //           g = first/last good pixel
    //           x = bad
    // 
    //           the first good pixel is the array index that has the value "center", 2
    //           the last good pixel has array index N-(5-2)+1
    //           eg. if N = 100, you want to use up to index 97
    //               100-3+1 = 98, and the loops use i < 98, meaning the last
    //               index you address is 97.
   
    unsigned int startCol = (*kiter)->getCtrX();
    unsigned int startRow = (*kiter)->getCtrY();
    unsigned int endCol   = (*citer)->getWidth()  - ((*kiter)->getWidth()  - (*kiter)->getCtrX()) + 1;
    unsigned int endRow   = (*citer)->getHeight() - ((*kiter)->getHeight() - (*kiter)->getCtrY()) + 1;

    std::vector<xy_locator> convolvedLocatorList;
    for (citer = convolvedImageList.begin(); citer != convolvedImageList.end(); ++citer) {
        convolvedLocatorList.push_back( (**citer).xy_at(startCol,startRow) );
    }
    xy_locator  imageToConvolveLocator    = imageToConvolve.xy_at(startCol, startRow);
    xy_locator  imageToNotConvolveLocator = imageToNotConvolve.xy_at(startCol, startRow);
    xyi_locator varianceLocator           = varianceImage.xy_at(startCol, startRow);

    // Unit test ImageSubtract_1.py should show
    // Image range : 9 9 -> 31 31 : 2804.000000 2798.191162
    logging::TTrace<8>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                       "Image range : %d %d -> %d %d : %f %f",
                       startCol, startRow, endCol, endRow, 
                       imageToConvolveLocator.image(), 
                       imageToNotConvolveLocator.image());

    std::pair<int, int> rowStep = std::make_pair(static_cast<int>(-(endCol-startCol)), 1);
    for (unsigned int row = startRow; row < endRow; ++row) {
        
        for (unsigned int col = startCol; col < endCol; ++col) {
            
            ImageT const ncImage          = imageToNotConvolveLocator.image();
            ImageT const ncVariance       = imageToNotConvolveLocator.variance();
            image::MaskPixel const ncMask = imageToNotConvolveLocator.mask();
            double const iVariance        = 1.0 / *varianceLocator;
            
            // Unit test ImageSubtract_1.py should show
            // Accessing image row 9 col 9  : 2798.191 23.426 0 1792.511475
            // Accessing image row 9 col 10 : 2805.171 23.459 0 1774.878662
            // ...
            // Accessing image row 9 col 30 : 2793.281 23.359 0 1779.194946
            // Accessing image row 10 col 9 : 2802.968 23.464 0 1770.467163
            // ...
            logging::TTrace<8>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                               "Accessing image row %d col %d : %.3f %.3f %d %f",
                               row, col, ncImage, ncVariance, ncMask, 1.0 * *varianceLocator);
            
            // kernel index i
            typename std::vector<xy_locator>::iterator 
                citeri = convolvedLocatorList.begin();
            typename std::vector<xy_locator>::iterator 
                citerE = convolvedLocatorList.end();

            for (int kidxi = 0; citeri != citerE; ++citeri, ++kidxi) {
                ImageT           const cdImagei = (*citeri).image();
                image::MaskPixel const cdMaski  = (*citeri).mask();
                if (cdMaski != 0) {
                    throw LSST_EXCEPT(exceptions::Exception, 
                                      str(boost::format("Accessing invalid pixel (%d) in computePsfMatchingKernelForFootprint") % 
                                          kidxi));
                }                
                
                // kernel index j
                typename std::vector<xy_locator>::iterator 
                    citerj = citeri;

                for (int kidxj = kidxi; citerj != citerE; ++citerj, ++kidxj) {
                    ImageT const cdImagej = (*citerj).image();
                    
                    M(kidxi, kidxj) += cdImagei * cdImagej * iVariance;
                } 
                
                B(kidxi) += ncImage * cdImagei * iVariance;
                
                // Constant background term; effectively j=kidxj+1 */
                M(kidxi, nParameters-1) += cdImagei * iVariance;
            } 
            
            // Background term; effectively i=kidxi+1 
            B(nParameters-1)                += ncImage * iVariance;
            M(nParameters-1, nParameters-1) += 1.0 * iVariance;
            
            // Step each accessor in column
            ++imageToConvolveLocator.x();
            ++imageToNotConvolveLocator.x();
            ++varianceLocator.x();
            for (int ki = 0; ki < nKernelParameters; ++ki) {
                ++convolvedLocatorList[ki].x();
            }             
            
        } // col
        
        // Get to next row, first col
        imageToConvolveLocator    += rowStep;
        imageToNotConvolveLocator += rowStep;

        // HACK UNTIL Ticket #647 FIXED
        varianceLocator            = varianceImage.xy_at(startCol, row+1);

        for (int ki = 0; ki < nKernelParameters; ++ki) {
            convolvedLocatorList[ki] += rowStep;
        }
        
    } // row

    
    /** @note If we are going to regularize the solution to M, this is the place
     * to do it 
     */
    
    // Fill in rest of M
    for (int kidxi=0; kidxi < nParameters; ++kidxi) 
        for (int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) 
            M(kidxj, kidxi) = M(kidxi, kidxj);
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                       "Total compute time before matrix inversions : %.2f s", time);

    //std::cout << "B eigen : " << B << std::endl;

    // To use Cholesky decomposition, the matrix needs to be symmetric (M is, by
    // design) and positive definite.  
    //
    // Eventually put a check in here to make sure its positive definite
    //
    Eigen::VectorXd Soln = Eigen::VectorXd::Zero(nParameters);;
    if (!( M.ldlt().solve(B, &Soln) )) {
        logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                           "Unable to determine kernel via Cholesky LDL^T");
        if (!( M.llt().solve(B, &Soln) )) {
            logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                               "Unable to determine kernel via Cholesky LL^T");
            if (!( M.lu().solve(B, &Soln) )) {
                logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                                   "Unable to determine kernel via LU");
                // LAST RESORT
                try {
                    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eVecValues(M);
                    Eigen::MatrixXd const& R = eVecValues.eigenvectors();
                    Eigen::VectorXd eValues  = eVecValues.eigenvalues();
                    
                    for (int i = 0; i != eValues.rows(); ++i) {
                        if (eValues(i) != 0.0) {
                            eValues(i) = 1.0/eValues(i);
                        }
                    }
                    
                    Soln = R*eValues.asDiagonal()*R.transpose()*B;
                } catch (exceptions::Exception& e) {
                    logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                                       "Unable to determine kernel via eigen-values");
                    
                    throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel solution");
                }
            }
        }
    }
    //std::cout << "Soln eigen : " << Soln << std::endl;
    //return;

    // Estimate of parameter uncertainties comes from the inverse of the
    // covariance matrix (noise spectrum).  Use Cholesky decomposition again.
    // Cholkesy:
    // Cov       =  L L^t
    // Cov^(-1)  = (L L^t)^(-1)
    //           = (L^T)^-1 L^(-1)
    Eigen::MatrixXd             Cov    = M.transpose() * M;
    Eigen::LLT<Eigen::MatrixXd> llt    = Cov.llt();
    Eigen::MatrixXd             Error2 = llt.matrixL().transpose().inverse() * llt.matrixL().inverse();
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                       "Total compute time after matrix inversions : %.2f s", time);
    
    // Translate from Eigen vectors into LSST classes
    unsigned int kCols = policy.getInt("kernelCols");
    unsigned int kRows = policy.getInt("kernelRows");
    std::vector<double> kValues(kCols*kRows);
    std::vector<double> kErrValues(kCols*kRows);
    for (unsigned int row = 0, idx = 0; row < kRows; row++) {
        for (unsigned int col = 0; col < kCols; col++, idx++) {
            
            // Insanity checking
            if (std::isnan( Soln(idx) )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel solution (nan)");
            }
            if (std::isnan( Error2(idx, idx) )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
            }
            if (Error2(idx, idx) < 0.0) {
                throw LSST_EXCEPT(exceptions::Exception,
                                  str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                                      Error2(idx, idx)
                                      ));
            }
            
            kValues[idx]    = Soln(idx);
            kErrValues[idx] = sqrt(Error2(idx, idx));
        }
    }
    //kernelPtr = boost::shared_ptr<math::Kernel> 
    //(new math::LinearCombinationKernel(kernelInBasisList, kValues));
    //kernelErrorPtr = boost::shared_ptr<math::Kernel> 
    //(new math::LinearCombinationKernel(kernelInBasisList, kErrValues));
    kernelPtr.reset( new math::LinearCombinationKernel(kernelInBasisList, kValues) );
    kernelErrorPtr.reset( new math::LinearCombinationKernel(kernelInBasisList, kErrValues) );
    
    // Estimate of Background and Background Error */
    if (std::isnan( Error2(nParameters-1, nParameters-1) )) {
        throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
    }
    if (Error2(nParameters-1, nParameters-1) < 0.0) {
        throw LSST_EXCEPT(exceptions::Exception, 
                          str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                              Error2(nParameters-1, nParameters-1) 
                              ));
    }
    background      = Soln(nParameters-1);
    backgroundError = sqrt(Error2(nParameters-1, nParameters-1));
}

#if USE_VW
/** 
 * @brief Computes a single Kernel (Model 1) around a single subimage.
 *
 * Accepts two MaskedImages, generally subimages of a larger image, one of which
 * is to be convolved to match the other.  The output Kernel is generated using
 * an input list of basis Kernels by finding the coefficients in front of each
 * basis.  This version accepts an input variance image, and uses VW for the
 * matrices.  
 *
 * @return Vector of coefficients representing the relative contribution of
 * each input basis function.
 *
 * @return Differential background offset between the two images
 *
 * @ingroup diffim
 */
template <typename ImageT, typename VarT>
void diffim::computePsfMatchingKernelForFootprintVW(
    double                          &background,
    double                          &backgroundError,
    boost::shared_ptr<math::Kernel> &kernelPtr,
    boost::shared_ptr<math::Kernel> &kernelErrorPtr,
    image::MaskedImage<ImageT>         const& imageToConvolve,    
    image::MaskedImage<ImageT>         const& imageToNotConvolve, 
    image::Image<VarT>                 const& varianceImage,      
    math::KernelList<math::Kernel>     const& kernelInBasisList,  
    lsst::pex::policy::Policy          const& policy
    ) { 

    typedef typename image::MaskedImage<ImageT>::xy_locator xy_locator;
    typedef typename image::Image<VarT>::xy_locator xyi_locator;

    // grab mask bits from the image to convolve, since that is what we'll be operating on
    int edgeMaskBit = imageToConvolve.getMask()->getMaskPlane("EDGE");
    
    int nKernelParameters = 0;
    int nBackgroundParameters = 0;
    int nParameters = 0;
    
    boost::timer t;
    double time;
    t.restart();

    nKernelParameters     = kernelInBasisList.size();
    nBackgroundParameters = 1;
    nParameters           = nKernelParameters + nBackgroundParameters;

    vw::math::Vector<double> B(nParameters);
    vw::math::Matrix<double> M(nParameters, nParameters);
    for (unsigned int i = nParameters; i--;) {
        B(i) = 0;
        for (unsigned int j = nParameters; j--;) {
            M(i,j) = 0;
        }
    }
    
    std::vector<boost::shared_ptr<image::MaskedImage<ImageT> > > convolvedImageList(nKernelParameters);
    typename std::vector<boost::shared_ptr<image::MaskedImage<ImageT> > >::iterator 
        citer = convolvedImageList.begin();
    std::vector<boost::shared_ptr<math::Kernel> >::const_iterator 
        kiter = kernelInBasisList.begin();
    
    // Create C_ij in the formalism of Alard & Lupton */
    for (int idx=0; kiter != kernelInBasisList.end(); ++kiter, ++citer, ++idx) {
        
        /* NOTE : we could also *precompute* the entire template image convolved with these functions */
        /*        and save them somewhere to avoid this step each time.  however, our paradigm is to */
        /*        compute whatever is needed on the fly.  hence this step here. */
        image::MaskedImage<ImageT> image(imageToConvolve.getDimensions());
        math::convolve(image,
                       imageToConvolve,
                       **kiter,
                       false,
                       edgeMaskBit);
        // image.writeFits(str(boost::format("c%d") % idx));
        boost::shared_ptr<image::MaskedImage<ImageT> > imagePtr( new image::MaskedImage<ImageT>(image) );
        *citer = imagePtr;
    } 
    
    kiter = kernelInBasisList.begin();
    citer = convolvedImageList.begin();

    // Ignore buffers around edge of convolved images :
    //
    // If the kernel has width 5, it has center pixel 2.  The first good pixel
    // is the (5-2)=3rd pixel, which is array index 2, and ends up being the
    // index of the central pixel.
    //
    // You also have a buffer of unusable pixels on the other side, numbered
    // width-center-1.  The last good usable pixel is N-width+center+1.

    // Example : the kernel is width = 5, center = 2
    //
    //     ---|---|-c-|---|---|
    //          
    //           the image is width = N
    //           convolve this with the kernel, and you get
    //
    //    |-x-|-x-|-g-|---|---| ... |---|---|-g-|-x-|-x-|
    //
    //           g = first/last good pixel
    //           x = bad
    // 
    //           the first good pixel is the array index that has the value "center", 2
    //           the last good pixel has array index N-(5-2)+1
    //           eg. if N = 100, you want to use up to index 97
    //               100-3+1 = 98, and the loops use i < 98, meaning the last
    //               index you address is 97.
   
    unsigned int startCol = (*kiter)->getCtrX();
    unsigned int startRow = (*kiter)->getCtrY();
    unsigned int endCol   = (*citer)->getWidth()  - ((*kiter)->getWidth()  - (*kiter)->getCtrX()) + 1;
    unsigned int endRow   = (*citer)->getHeight() - ((*kiter)->getHeight() - (*kiter)->getCtrY()) + 1;

    std::vector<xy_locator> convolvedLocatorList;
    for (citer = convolvedImageList.begin(); citer != convolvedImageList.end(); ++citer) {
        convolvedLocatorList.push_back( (**citer).xy_at(startCol,startRow) );
    }
    xy_locator  imageToConvolveLocator    = imageToConvolve.xy_at(startCol, startRow);
    xy_locator  imageToNotConvolveLocator = imageToNotConvolve.xy_at(startCol, startRow);
    xyi_locator varianceLocator           = varianceImage.xy_at(startCol, startRow);

    // Unit test ImageSubtract_1.py should show
    // Image range : 9 9 -> 31 31 : 2804.000000 2798.191162
    logging::TTrace<8>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                       "Image range : %d %d -> %d %d : %f %f",
                       startCol, startRow, endCol, endRow, 
                       imageToConvolveLocator.image(), 
                       imageToNotConvolveLocator.image());

    std::pair<int, int> rowStep = std::make_pair(static_cast<int>(-(endCol-startCol)), 1);
    for (unsigned int row = startRow; row < endRow; ++row) {
        
        for (unsigned int col = startCol; col < endCol; ++col) {
            
            ImageT const ncImage          = imageToNotConvolveLocator.image();
            ImageT const ncVariance       = imageToNotConvolveLocator.variance();
            image::MaskPixel const ncMask = imageToNotConvolveLocator.mask();
            double const iVariance        = 1.0 / *varianceLocator;
            
            // Unit test ImageSubtract_1.py should show
            // Accessing image row 9 col 9  : 2798.191 23.426 0 1792.511475
            // Accessing image row 9 col 10 : 2805.171 23.459 0 1774.878662
            // ...
            // Accessing image row 9 col 30 : 2793.281 23.359 0 1779.194946
            // Accessing image row 10 col 9 : 2802.968 23.464 0 1770.467163
            // ...
            logging::TTrace<8>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                               "Accessing image row %d col %d : %.3f %.3f %d %f",
                               row, col, ncImage, ncVariance, ncMask, 1.0 * *varianceLocator);
            
            // kernel index i
            typename std::vector<xy_locator>::iterator 
                citeri = convolvedLocatorList.begin();
            typename std::vector<xy_locator>::iterator 
                citerE = convolvedLocatorList.end();

            for (int kidxi = 0; citeri != citerE; ++citeri, ++kidxi) {
                ImageT           const cdImagei = (*citeri).image();
                image::MaskPixel const cdMaski  = (*citeri).mask();
                if (cdMaski != 0) {
                    throw LSST_EXCEPT(exceptions::Exception, 
                                      str(boost::format("Accessing invalid pixel (%d) in computePsfMatchingKernelForFootprint") % 
                                          kidxi));
                }                
                
                // kernel index j
                typename std::vector<xy_locator>::iterator 
                    citerj = citeri;

                for (int kidxj = kidxi; citerj != citerE; ++citerj, ++kidxj) {
                    ImageT const cdImagej    = (*citerj).image();
                    M[kidxi][kidxj]   += cdImagei * cdImagej * iVariance;
                } 
                
                B[kidxi] += ncImage * cdImagei * iVariance;
                
                // Constant background term; effectively j=kidxj+1 */
                M[kidxi][nParameters-1] += cdImagei * iVariance;
            } 
            
            // Background term; effectively i=kidxi+1 
            B[nParameters-1]                += ncImage * iVariance;
            M[nParameters-1][nParameters-1] += 1.0 * iVariance;
            
            // Step each accessor in column
            ++imageToConvolveLocator.x();
            ++imageToNotConvolveLocator.x();
            ++varianceLocator.x();
            for (int ki = 0; ki < nKernelParameters; ++ki) {
                ++convolvedLocatorList[ki].x();
            }             
            
        } // col

        // Get to next row, first col
        imageToConvolveLocator    += rowStep;
        imageToNotConvolveLocator += rowStep;

        // HACK UNTIL Ticket #647 FIXED
        varianceLocator            = varianceImage.xy_at(startCol, row+1);

        for (int ki = 0; ki < nKernelParameters; ++ki) {
            convolvedLocatorList[ki] += rowStep;
        }
        
    } // row

    
    /** @note If we are going to regularize the solution to M, this is the place
     * to do it 
     */
    
    // Fill in rest of M
    for (int kidxi=0; kidxi < nParameters; ++kidxi) 
        for (int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) 
            M[kidxj][kidxi] = M[kidxi][kidxj];
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                       "Total compute time before matrix inversions : %.2f s", time);

    // Invert using VW's internal method
    vw::math::Vector<double> Soln      = vw::math::least_squares(M, B);

    //std::cout << "M vw : " << M << std::endl;
    //std::cout << "B vw : " << B << std::endl;
    //std::cout << "Soln vw : " << Soln << std::endl;
    //return;
    
    // Additional gymnastics to get the parameter uncertainties
    vw::math::Matrix<double> Mt        = vw::math::transpose(M);
    vw::math::Matrix<double> MtM       = Mt * M;
    vw::math::Matrix<double> Error2    = vw::math::pseudoinverse(MtM);
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                       "Total compute time after matrix inversions : %.2f s", time);
    
    // Translate from VW vectors into LSST classes
    unsigned int kCols = policy.getInt("kernelCols");
    unsigned int kRows = policy.getInt("kernelRows");
    std::vector<double> kValues(kCols*kRows);
    std::vector<double> kErrValues(kCols*kRows);
    for (unsigned int row = 0, idx = 0; row < kRows; row++) {
        for (unsigned int col = 0; col < kCols; col++, idx++) {
            
            // Insanity checking
            if (std::isnan( Soln[idx] )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel solution (nan)");
            }
            if (std::isnan( Error2[idx][idx] )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
            }
            if (Error2[idx][idx] < 0.0) {
                throw LSST_EXCEPT(exceptions::Exception,
                                  str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                                      Error2[idx][idx]
                                      ));
            }
            
            kValues[idx]    = Soln[idx];
            kErrValues[idx] = sqrt(Error2[idx][idx]);
        }
    }
    //kernelPtr = boost::shared_ptr<math::Kernel> 
    //(new math::LinearCombinationKernel(kernelInBasisList, kValues));
    //kernelErrorPtr = boost::shared_ptr<math::Kernel> 
    //(new math::LinearCombinationKernel(kernelInBasisList, kErrValues));
    kernelPtr.reset( new math::LinearCombinationKernel(kernelInBasisList, kValues) );
    kernelErrorPtr.reset( new math::LinearCombinationKernel(kernelInBasisList, kErrValues) );
    
    
    // Estimate of Background and Background Error */
    if (std::isnan( Error2[nParameters-1][nParameters-1] )) {
        throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
    }
    if (Error2[nParameters-1][nParameters-1] < 0.0) {
        throw LSST_EXCEPT(exceptions::Exception, 
                          str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                              Error2[nParameters-1][nParameters-1]
                              ));
    }
    background      = Soln[nParameters-1];
    backgroundError = sqrt(Error2[nParameters-1][nParameters-1]);
}
#endif

/** 
 * @brief Computes a single Kernel (Model 1) around a single subimage.
 *
 * Accepts two MaskedImages, generally subimages of a larger image, one of which
 * is to be convolved to match the other.  The output Kernel is generated using
 * an input list of basis Kernels by finding the coefficients in front of each
 * basis.  This version accepts an input variance image.
 *
 * @return Vector of coefficients representing the relative contribution of
 * each input basis function.
 *
 * @return Differential background offset between the two images
 *
 * @ingroup diffim
 */
template <typename ImageT, typename MaskT>
void diffim::computePsfMatchingKernelForFootprint(
    image::MaskedImage<ImageT, MaskT>         const &imageToConvolve,    
    image::MaskedImage<ImageT, MaskT>         const &imageToNotConvolve, 
    image::Image<ImageT, MaskT>               const &varianceImage,      
    math::KernelList<math::Kernel> const &kernelInBasisList,  
    lsst::pex::policy::Policy                  &policy,
    boost::shared_ptr<math::Kernel> &kernelPtr,
    boost::shared_ptr<math::Kernel> &kernelErrorPtr,
    double                                     &background,
    double                                     &backgroundError
    ) { 
    
    // grab mask bits from the image to convolve, since that is what we'll be operating on
    int edgeMaskBit = imageToConvolve.getMask()->getMaskPlane("EDGE");
    
    int nKernelParameters = 0;
    int nBackgroundParameters = 0;
    int nParameters = 0;
    
    boost::timer t;
    double time;
    t.restart();
    
    nKernelParameters     = kernelInBasisList.size();
    nBackgroundParameters = 1;
    nParameters           = nKernelParameters + nBackgroundParameters;
    
    vw::math::Vector<double> B(nParameters);
    vw::math::Matrix<double> M(nParameters, nParameters);
    for (unsigned int i = nParameters; i--;) {
        B(i) = 0;
        for (unsigned int j = nParameters; j--;) {
            M(i,j) = 0;
        }
    }
    
    std::vector<boost::shared_ptr<image::MaskedImage<ImageT, MaskT> > > convolvedImageList(nKernelParameters);
    typename std::vector<boost::shared_ptr<image::MaskedImage<ImageT, MaskT> > >::iterator 
        citer = convolvedImageList.begin();
    std::vector<boost::shared_ptr<math::Kernel> >::const_iterator 
        kiter = kernelInBasisList.begin();
    
    // Create C_ij in the formalism of Alard & Lupton
    for (; kiter != kernelInBasisList.end(); ++kiter, ++citer) {
        
        /* NOTE : we could also *precompute* the entire template image convolved with these functions */
        /*        and save them somewhere to avoid this step each time.  however, our paradigm is to */
        /*        compute whatever is needed on the fly.  hence this step here. */
        boost::shared_ptr<image::MaskedImage<ImageT, MaskT> > imagePtr(
                                                                       new image::MaskedImage<ImageT, MaskT>
                                                                       (math::convolveNew(imageToConvolve, **kiter, edgeMaskBit, false))
                                                                       );
        
        *citer = imagePtr;
    } 
    
    kiter = kernelInBasisList.begin();
    citer = convolvedImageList.begin();
    unsigned int startCol = (*kiter)->getCtrCol();
    unsigned int startRow = (*kiter)->getCtrRow();
    
    unsigned int endCol   = (*citer)->getCols() - ((*kiter)->getCols() - (*kiter)->getCtrCol()) + 1;
    unsigned int endRow   = (*citer)->getRows() - ((*kiter)->getRows() - (*kiter)->getCtrRow()) + 1;
    
    std::vector<image::MaskedPixelAccessor<ImageT, MaskT> > convolvedAccessorRowList;
    for (citer = convolvedImageList.begin(); citer != convolvedImageList.end(); ++citer) {
        convolvedAccessorRowList.push_back(image::MaskedPixelAccessor<ImageT, MaskT>(**citer));
    }
    
    // An accessor for each input image; address rows and cols separately
    image::MaskedPixelAccessor<ImageT, MaskT> imageToConvolveRow(imageToConvolve);
    image::MaskedPixelAccessor<ImageT, MaskT> imageToNotConvolveRow(imageToNotConvolve);
    image::MaskedPixelAccessor<ImageT, MaskT> varianceRow(varianceImage);
    
    // Address input images
    imageToConvolveRow.advance(startCol, startRow);
    imageToNotConvolveRow.advance(startCol, startRow);
    varianceRow.advance(startCol, startRow);
    // Address kernel images
    for (int ki = 0; ki < nKernelParameters; ++ki) {
        convolvedAccessorRowList[ki].advance(startCol, startRow);
    }
    
    for (unsigned int row = startRow; row < endRow; ++row) {
        // An accessor for each convolution plane 
        std::vector<image::MaskedPixelAccessor<ImageT, MaskT> > convolvedAccessorColList = convolvedAccessorRowList;
        
        // An accessor for each input image
        image::MaskedPixelAccessor<ImageT, MaskT> imageToConvolveCol = imageToConvolveRow;
        image::MaskedPixelAccessor<ImageT, MaskT> imageToNotConvolveCol = imageToNotConvolveRow;
        image::MaskedPixelAccessor<ImageT, MaskT> varianceCol = varianceRow;
        
        for (unsigned int col = startCol; col < endCol; ++col) {
            
            ImageT ncImage    = *imageToNotConvolveCol.image;
            ImageT ncVariance = *imageToNotConvolveCol.variance;
            MaskT  ncMask     = *imageToNotConvolveCol.mask;
            
            ImageT iVariance  = 1.0 / *varianceCol.variance;
            
            logging::TTrace<8>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                               "Accessing image row %d col %d : %.3f %.3f %d",
                               row, col, ncImage, ncVariance, ncMask);
            
            // kernel index i
            typename std::vector<image::MaskedPixelAccessor<ImageT, MaskT> >::iterator
                kiteri   = convolvedAccessorColList.begin();
            typename std::vector<image::MaskedPixelAccessor<ImageT, MaskT> >::iterator
                kiterEnd = convolvedAccessorColList.end();
            
            for (int kidxi = 0; kiteri != kiterEnd; ++kiteri, ++kidxi) {
                ImageT cdImagei   = *kiteri->image;
                
                /** @note Commenting in these additional pixel accesses yields
                 * an additional second of run-time per kernel with opt=1 at 2.8
                 * GHz
                 */
                
                /* ignore unnecessary pixel accesses 
                   ImageT cdVariancei = *kiteri->variance;
                   MaskT  cdMaski     = *kiteri->mask;
                   logging::TTrace<8>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                   "Accessing convolved image %d : %.3f %.3f %d",
                   kidxi, cdImagei, cdVariancei, cdMaski);
                */
                
                // kernel index j
                typename std::vector<image::MaskedPixelAccessor<ImageT, MaskT> >::iterator kiterj = kiteri;
                
                for (int kidxj = kidxi; kiterj != kiterEnd; ++kiterj, ++kidxj) {
                    ImageT cdImagej   = *kiterj->image;
                    M[kidxi][kidxj]   += cdImagei * cdImagej * iVariance;
                } 
                
                B[kidxi] += ncImage * cdImagei * iVariance;
                
                // Constant background term; effectively j=kidxj+1 */
                M[kidxi][nParameters-1] += cdImagei * iVariance;
            } 
            
            // Background term; effectively i=kidxi+1 */
            B[nParameters-1]                += ncImage * iVariance;
            M[nParameters-1][nParameters-1] += 1.0 * iVariance;
            
            logging::TTrace<8>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                               "Background terms : %.3f %.3f",
                               B[nParameters-1], M[nParameters-1][nParameters-1]);
            
            // Step each accessor in column 
            imageToConvolveCol.nextCol();
            imageToNotConvolveCol.nextCol();
            varianceCol.nextCol();
            for (int ki = 0; ki < nKernelParameters; ++ki) {
                convolvedAccessorColList[ki].nextCol();
            }             
            
        } // col
        
        // Step each accessor in row
        imageToConvolveRow.nextRow();
        imageToNotConvolveRow.nextRow();
        varianceRow.nextRow();
        for (int ki = 0; ki < nKernelParameters; ++ki) {
            convolvedAccessorRowList[ki].nextRow();
        }
        
    } // row
    
    /** @note If we are going to regularize the solution to M, this is the place
     * to do it 
     */
    
    // Fill in rest of M
    for (int kidxi=0; kidxi < nParameters; ++kidxi) 
        for (int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) 
            M[kidxj][kidxi] = M[kidxi][kidxj];
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                       "Total compute time before matrix inversions : %.2f s", time);
    
    /* Invert using SVD and Pseudoinverse : This is a full second slower per
     * kernel than vw::math::least_squares, compiled at opt=1 at 2.8 GHz.
     
     vw::math::Matrix<double> Minv = vw::math::pseudoinverse(M);
     vw::math::Vector<double> Soln = Minv * B;
    */     
    
    // Invert using VW's internal method
    vw::math::Vector<double> kSolution = vw::math::least_squares(M, B);
    
    // Additional gymnastics to get the parameter uncertainties
    vw::math::Matrix<double> Mt        = vw::math::transpose(M);
    vw::math::Matrix<double> MtM       = Mt * M;
    vw::math::Matrix<double> kError    = vw::math::pseudoinverse(MtM);
    
    /*
      NOTE : for any real kernels I have looked at, these solutions have agreed
      exactly.  However, when designing the testDeconvolve unit test with
      hand-built gaussians as objects and non-realistic noise, the solutions did
      *NOT* agree.
      
      std::cout << "Soln : " << Soln << std::endl;
      std::cout << "Soln2 : " << Soln2 << std::endl;
    */
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                       "Total compute time after matrix inversions : %.2f s", time);
    
    // Translate from VW vectors into LSST classes
    unsigned int kCols = policy.getInt("kernelCols");
    unsigned int kRows = policy.getInt("kernelRows");
    vector<double> kValues(kCols*kRows);
    vector<double> kErrValues(kCols*kRows);
    for (unsigned int row = 0, idx = 0; row < kRows; row++) {
        for (unsigned int col = 0; col < kCols; col++, idx++) {
            
            // Insanity checking
            if (std::isnan(kSolution[idx])) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel solution (nan)");
            }
            if (std::isnan(kError[idx][idx])) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
            }
            if (kError[idx][idx] < 0.0) {
                throw LSST_EXCEPT(exceptions::Exception
                                  str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                                      kError[idx][idx]
                                      ));
            }
            
            kValues[idx]    = kSolution[idx];
            kErrValues[idx] = sqrt(kError[idx][idx]);
        }
    }
    kernelPtr = boost::shared_ptr<math::Kernel> (
                                                 new math::LinearCombinationKernel(kernelInBasisList, kValues)
                                                 );
    kernelErrorPtr = boost::shared_ptr<math::Kernel> (
                                                      new math::LinearCombinationKernel(kernelInBasisList, kErrValues)
                                                      );
    
    // Estimate of Background and Background Error */
    if (std::isnan(kError[nParameters-1][nParameters-1])) {
        throw LSST_EXCEPT(exceptions::Exception, "Unable to determine background uncertainty (nan)");
    }
    if (kError[nParameters-1][nParameters-1] < 0.0) {
        throw LSST_EXCEPT(exceptions::Exception,
                          str(boost::format("Unable to determine background uncertainty, negative variance (%.3e)") % 
                              kError[nParameters-1][nParameters-1]
                              ));
    }
    background      = kSolution[nParameters-1];
    backgroundError = sqrt(kError[nParameters-1][nParameters-1]);
}


/** 
 * @brief Computes a single Kernel (Model 1) around a single subimage.
 *
 * Accepts two MaskedImages, generally subimages of a larger image, one of which
 * is to be convolved to match the other.  The output Kernel is generated using
 * an input list of basis Kernels by finding the coefficients in front of each
 * basis.  An old version with comments and considerations for posterity.
 *
 * @return Vector of coefficients representing the relative contribution of
 * each input basis function.
 *
 * @return Differential background offset between the two images
 *
 * @ingroup diffim
 */
template <typename ImageT, typename MaskT>
std::vector<double> diffim::computePsfMatchingKernelForFootprint_Legacy(
    double &background,                                                            
    image::MaskedImage<ImageT, MaskT> const &imageToConvolve,           
    image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve,        
    math::KernelList<math::Kernel> const &kernelInBasisList, 
    lsst::pex::policy::Policy &policy                                              
    ) { 
    
    /* grab mask bits from the image to convolve, since that is what we'll be operating on */
    int edgeMaskBit = imageToConvolve.getMask()->getMaskPlane("EDGE");
    
    int nKernelParameters = 0;
    int nBackgroundParameters = 0;
    int nParameters = 0;
    
    boost::timer t;
    double time;
    t.restart();
    
    logging::TTrace<6>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                       "Entering subroutine computePsfMatchingKernelForFootprint");
    
    /* We assume that each kernel in the Set has 1 parameter you fit for */
    nKernelParameters = kernelInBasisList.size();
    /* Or, we just assume that across a single kernel, background 0th order.  This quite makes sense. */
    nBackgroundParameters = 1;
    /* Total number of parameters */
    nParameters = nKernelParameters + nBackgroundParameters;
    
    vw::math::Vector<double> B(nParameters);
    vw::math::Matrix<double> M(nParameters, nParameters);
    for (unsigned int i = nParameters; i--;) {
        B(i) = 0;
        for (unsigned int j = nParameters; j--;) {
            M(i,j) = 0;
        }
    }
    
    /* convolve creates a MaskedImage, push it onto the back of the Vector */
    /* need to use shared pointers because MaskedImage copy does not work */
    std::vector<boost::shared_ptr<image::MaskedImage<ImageT, MaskT> > > convolvedImageList(nKernelParameters);
    /* and an iterator over this */
    typename std::vector<boost::shared_ptr<image::MaskedImage<ImageT, MaskT> > >::iterator citer = convolvedImageList.begin();
    
    /* Iterator for input kernel basis */
    std::vector<boost::shared_ptr<math::Kernel> >::const_iterator kiter = kernelInBasisList.begin();
    /* Create C_ij in the formalism of Alard & Lupton */
    for (; kiter != kernelInBasisList.end(); ++kiter, ++citer) {
        
        logging::TTrace<7>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                           "Convolving an Object with Basis");
        
        /* NOTE : we could also *precompute* the entire template image convolved with these functions */
        /*        and save them somewhere to avoid this step each time.  however, our paradigm is to */
        /*        compute whatever is needed on the fly.  hence this step here. */
        boost::shared_ptr<image::MaskedImage<ImageT, MaskT> > imagePtr(
                                                                       new image::MaskedImage<ImageT, MaskT>
                                                                       (math::convolveNew(imageToConvolve, **kiter, edgeMaskBit, false))
                                                                       );
        
        logging::TTrace<7>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                           "Convolved an Object with Basis");
        
        *citer = imagePtr;
        
    } 
    
    /* NOTE ABOUT CONVOLUTION : */
    /* getCtrCol:getCtrRow pixels are masked on the left:bottom side */
    /* getCols()-getCtrCol():getRows()-getCtrRow() masked on right/top side */
    /* */
    /* The convolved image and the input image are by default the same size, so */
    /* we offset our initial pixel references by the same amount */
    kiter = kernelInBasisList.begin();
    citer = convolvedImageList.begin();
    unsigned int startCol = (*kiter)->getCtrCol();
    unsigned int startRow = (*kiter)->getCtrRow();
    /* NOTE - I determined I needed this +1 by eye */
    unsigned int endCol   = (*citer)->getCols() - ((*kiter)->getCols() - (*kiter)->getCtrCol()) + 1;
    unsigned int endRow   = (*citer)->getRows() - ((*kiter)->getRows() - (*kiter)->getCtrRow()) + 1;
    /* NOTE - we need to enforce that the input images are large enough */
    /* How about some multiple of the PSF FWHM?  Or second moments? */
    
    /* An accessor for each convolution plane */
    /* NOTE : MaskedPixelAccessor has no empty constructor, therefore we need to push_back() */
    std::vector<image::MaskedPixelAccessor<ImageT, MaskT> > convolvedAccessorRowList;
    for (citer = convolvedImageList.begin(); citer != convolvedImageList.end(); ++citer) {
        convolvedAccessorRowList.push_back(image::MaskedPixelAccessor<ImageT, MaskT>(**citer));
    }
    
    /* An accessor for each input image; address rows and cols separately */
    image::MaskedPixelAccessor<ImageT, MaskT> imageToConvolveRow(imageToConvolve);
    image::MaskedPixelAccessor<ImageT, MaskT> imageToNotConvolveRow(imageToNotConvolve);
    
    /* Take into account buffer for kernel images */
    imageToConvolveRow.advance(startCol, startRow);
    imageToNotConvolveRow.advance(startCol, startRow);
    for (int ki = 0; ki < nKernelParameters; ++ki) {
        convolvedAccessorRowList[ki].advance(startCol, startRow);
    }
    
    for (unsigned int row = startRow; row < endRow; ++row) {
        
        /* An accessor for each convolution plane */
        std::vector<image::MaskedPixelAccessor<ImageT, MaskT> > convolvedAccessorColList = convolvedAccessorRowList;
        
        /* An accessor for each input image; places the col accessor at the correct row */
        image::MaskedPixelAccessor<ImageT, MaskT> imageToConvolveCol = imageToConvolveRow;
        image::MaskedPixelAccessor<ImageT, MaskT> imageToNotConvolveCol = imageToNotConvolveRow;
        
        for (unsigned int col = startCol; col < endCol; ++col) {
            
            ImageT ncImage   = *imageToNotConvolveCol.image;
            ImageT ncVariance = *imageToNotConvolveCol.variance;
            MaskT  ncMask     = *imageToNotConvolveCol.mask;
            
            ImageT cVariance  = *imageToConvolveCol.variance;
            
            /* Variance for a particlar pixel; do we use this variance of the */
            /* input data, or include the variance after its been convolved with */
            /* the basis?  For now, use the average of the input varianes. */
            ImageT iVariance  = 1.0 / (cVariance + ncVariance);
            
            logging::TTrace<8>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                               "Accessing image row %d col %d : %.3f %.3f %d",
                               row, col, ncImage, ncVariance, ncMask);
            
            /* kernel index i */
            typename std::vector<image::MaskedPixelAccessor<ImageT, MaskT> >::iterator
                kiteri   = convolvedAccessorColList.begin();
            typename std::vector<image::MaskedPixelAccessor<ImageT, MaskT> >::iterator
                kiterEnd = convolvedAccessorColList.end();
            
            
            for (int kidxi = 0; kiteri != kiterEnd; ++kiteri, ++kidxi) {
                ImageT cdImagei   = *kiteri->image;
                
                ImageT cdVariancei = *kiteri->variance;
                MaskT  cdMaski     = *kiteri->mask;
                logging::TTrace<8>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                                   "Accessing convolved image %d : %.3f %.3f %d",
                                   kidxi, cdImagei, cdVariancei, cdMaski);
                
                /* kernel index j  */
                typename std::vector<image::MaskedPixelAccessor<ImageT, MaskT> >::iterator kiterj = kiteri;
                for (int kidxj = kidxi; kiterj != kiterEnd; ++kiterj, ++kidxj) {
                    ImageT cdImagej   = *kiterj->image;
                    
                    /* NOTE - These inner trace statements can ENTIRELY kill the run time */
                    ImageT cdVariancej = *kiterj->variance;
                    MaskT  cdMaskj     = *kiterj->mask;
                    logging::TTrace<8>("lsst.ip.diffim.computePsfMatchingKernelForFootprint",
                                       "Accessing convolved image %d : %.3f %.3f %d",
                                       kidxj, cdImagej, cdVariancej, cdMaskj);
                    
                    M[kidxi][kidxj] += cdImagei * cdImagej * iVariance;
                } 
                
                B[kidxi] += ncImage * cdImagei * iVariance;
                
                /* Constant background term; effectively j=kidxj+1 */
                M[kidxi][nParameters-1] += cdImagei * iVariance;
            } 
            
            /* Background term; effectively i=kidxi+1 */
            B[nParameters-1] += ncImage * iVariance;
            M[nParameters-1][nParameters-1] += 1.0 * iVariance;
            
            logging::TTrace<7>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                               "Background terms : %.3f %.3f",
                               B[nParameters-1], M[nParameters-1][nParameters-1]);
            
            /* Step each accessor in column */
            imageToConvolveCol.nextCol();
            imageToNotConvolveCol.nextCol();
            for (int ki = 0; ki < nKernelParameters; ++ki) {
                convolvedAccessorColList[ki].nextCol();
            }             
            
        } /* col */
        
        /* Step each accessor in row */
        imageToConvolveRow.nextRow();
        imageToNotConvolveRow.nextRow();
        for (int ki = 0; ki < nKernelParameters; ++ki) {
            convolvedAccessorRowList[ki].nextRow();
        }
        
    } /* row */
    
    /* NOTE: If we are going to regularize the solution to M, this is the place to do it */
    
    /* Fill in rest of M */
    for (int kidxi=0; kidxi < nParameters; ++kidxi) 
        for (int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) 
            M[kidxj][kidxi] = M[kidxi][kidxj];
    
#if DEBUG_MATRIX
    std::cout << "B : " << B << std::endl;
    std::cout << "M : " << M << std::endl;
#endif
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                       "Total compute time before matrix inversions : %.2f s", time);
    
    /* Invert using SVD and Pseudoinverse */
    vw::math::Matrix<double> Minv;
    Minv = vw::math::pseudoinverse(M);
    /*Minv = vw::math::inverse(M); */
    
#if DEBUG_MATRIX
    std::cout << "Minv : " << Minv << std::endl;
#endif
    
    /* Solve for x in Mx = B */
    vw::math::Vector<double> Soln = Minv * B;
    
#if DEBUG_MATRIX
    std::cout << "Solution : " << Soln << std::endl;
#endif
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                       "Total compute time after matrix inversions : %.2f s", time);
    
    /* Translate from VW std::vectors to std std::vectors */
    std::vector<double> kernelCoeffs(kernelInBasisList.size());
    for (int ki = 0; ki < nKernelParameters; ++ki) {
        kernelCoeffs[ki] = Soln[ki];
    }
    background = Soln[nParameters-1];
    
    logging::TTrace<6>("lsst.ip.diffim.computePsfMatchingKernelForFootprint", 
                       "Leaving subroutine computePsfMatchingKernelForFootprint");
    
    return kernelCoeffs;
}
#endif

