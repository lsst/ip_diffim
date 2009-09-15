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

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>

// NOTE -  trace statements >= 6 can ENTIRELY kill the run time
// #define LSST_MAX_TRACE 5

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
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> const& basisList
    ) :
    _basisList(basisList),
    _H(Eigen::MatrixXd()),
    _background(0.),
    _backgroundError(0.),
    _kernel(boost::shared_ptr<lsst::afw::math::Kernel>()),
    _kernelError(boost::shared_ptr<lsst::afw::math::Kernel>())
{;}

template <typename ImageT, typename VarT>
diffim::PsfMatchingFunctor<ImageT, VarT>::PsfMatchingFunctor(
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> const& basisList,
    Eigen::MatrixXd const H
    ) :
    _basisList(basisList),
    _H(H),
    _background(0.),
    _backgroundError(0.),
    _kernel(boost::shared_ptr<lsst::afw::math::Kernel>()),
    _kernelError(boost::shared_ptr<lsst::afw::math::Kernel>())
{;}

//
// Public Member Functions
//

/** Create PSF matching kernel
 */
template <typename ImageT, typename VarT>
void diffim::PsfMatchingFunctor<ImageT, VarT>::apply(
    lsst::afw::image::Image<ImageT> const& imageToConvolve,    ///< Image to apply kernel to
    lsst::afw::image::Image<ImageT> const& imageToNotConvolve, ///< Image whose PSF you want to match to
    lsst::afw::image::Image<VarT>   const& varianceEstimate,   ///< Estimate of the variance per pixel
    lsst::pex::policy::Policy       const& policy              ///< Policy file
    ) {
    
    unsigned int const nKernelParameters     = _basisList.size();
    unsigned int const nBackgroundParameters = 1;
    unsigned int const nParameters           = nKernelParameters + nBackgroundParameters;
    
    boost::timer t;
    t.restart();
    
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(nParameters, nParameters);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(nParameters);
    
    std::vector<boost::shared_ptr<image::Image<ImageT> > > convolvedImageList(nKernelParameters);
    typename std::vector<boost::shared_ptr<image::Image<ImageT> > >::iterator citer = convolvedImageList.begin();
    std::vector<boost::shared_ptr<math::Kernel> >::const_iterator kiter = _basisList.begin();
    
    // Create C_ij in the formalism of Alard & Lupton */
    for (; kiter != _basisList.end(); ++kiter, ++citer) {
        /*
         * NOTE : we could also *precompute* the entire template image convolved with these functions
         *        and save them somewhere to avoid this step each time.  however, our paradigm is to
         *        compute whatever is needed on the fly.  hence this step here.
         */
        *citer = typename image::Image<ImageT>::Ptr(new image::Image<ImageT>(imageToConvolve.getDimensions()));
        math::convolve(**citer, imageToConvolve, **kiter, false);
    } 

    double time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                       "Total compute time to do basis convolutions : %.2f s", time);
    t.restart();
     
    kiter = _basisList.begin();
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
    typename image::Image<ImageT>::xy_locator imageToConvolveLocator    = imageToConvolve.xy_at(startCol, startRow);
    typename image::Image<ImageT>::xy_locator imageToNotConvolveLocator = imageToNotConvolve.xy_at(startCol, startRow);
    xyi_locator varianceLocator                                         = varianceEstimate.xy_at(startCol, startRow);

    // Unit test ImageSubtract_1.py should show
    // Image range : 9 9 -> 31 31 : 2804.000000 2798.191162
    logging::TTrace<8>("lsst.ip.diffim.PsfMatchingFunctor.apply",
                       "Image range : %d %d -> %d %d : %f %f %f",
                       startCol, startRow, endCol, endRow, 
                       0 + *imageToConvolveLocator, 0 + *imageToNotConvolveLocator, 0 + *varianceLocator);

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
		    /*
		    logging::TTrace<8>("lsst.ip.diffim.PsfMatchingFunctor.apply",
				       "%f %f %f",
				       ncImage, cdImagei, cdImagej, 1./iVariance);
		    */
                } 
                B(kidxi) += ncImage*cdImagei*iVariance;

                // Background term; effectively j = kidxj + 1 */
                M(kidxi, nParameters-1) += cdImagei*iVariance;
            } 
            
            // Background term; effectively i = kidxi + 1, j = kidxj + 1
            B(nParameters-1)                += ncImage*iVariance;
            M(nParameters-1, nParameters-1) += 1.0*iVariance;
            
            // Step each accessor in column
            ++imageToConvolveLocator.x();
            ++imageToNotConvolveLocator.x();
            ++varianceLocator.x();
            for (unsigned int ki = 0; ki < nKernelParameters; ++ki) {
                ++convolvedLocatorList[ki].x();
            }             
            
        } // col
        
        // Get to next row, first col
        imageToConvolveLocator    += rowStep;
        imageToNotConvolveLocator += rowStep;
        varianceLocator           += rowStep;
        for (unsigned int ki = 0; ki < nKernelParameters; ++ki) {
            convolvedLocatorList[ki] += rowStep;
        }
        
    } // row
    
    // Fill in rest of M
    for (unsigned int kidxi=0; kidxi < nParameters; ++kidxi) {
        for (unsigned int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) {
            M(kidxj, kidxi) = M(kidxi, kidxj);
        }
    }
    
    std::cout << "M1 " << std::endl;
    std::cout << M << std::endl;
    std::cout << "B1 " << std::endl;
    std::cout << B << std::endl;

    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                       "Total compute time to step through pixels : %.2f s", time);
    t.restart();

    /* If the regularization matrix is here and not null, we use it by default */
    if ( (_H.rows() != 0) && (_H.cols() != 0) ) {
        double regularizationScaling = policy.getDouble("regularizationScaling");        
        /* 
           See N.R. 18.5 equation 18.5.8 for the solution to the regularized
           normal equations.  For M.x = b, and solving for x, 

           M -> (Mt.M + lambda*H)
           B -> (Mt.B)

           An estimate of lambda is NR 18.5.16

           lambda = Trace(Mt.M) / Tr(H)

         */

        Eigen::MatrixXd Mt = M.transpose();
        M = Mt * M;

        double lambda = M.trace() / _H.trace();
        lambda *= regularizationScaling;

        M = M + lambda * _H;
        B = Mt * B;
        logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                           "Applying kernel regularization with lambda = %.2e", lambda);
    }

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

    /* Fill in the kernel results */
    std::vector<double> kValues(nKernelParameters);
    for (unsigned int idx = 0; idx < nKernelParameters; idx++) {
        if (std::isnan( Soln(idx) )) {
            throw LSST_EXCEPT(exceptions::Exception, 
                              str(boost::format("Unable to determine kernel solution %d (nan)") % idx));
        }
        kValues[idx] = Soln(idx);
    }
    _kernel.reset( new math::LinearCombinationKernel(_basisList, kValues) );

    if (std::isnan( Soln(nParameters-1) )) {
        throw LSST_EXCEPT(exceptions::Exception, 
                          str(boost::format("Unable to determine background solution %d (nan)") % (nParameters-1)));
    }
    _background = Soln(nParameters-1);


    /* Calculate the kernel & background uncertainties only if requested */
    bool calculateKernelUncertainty = policy.getBool("calculateKernelUncertainty");            
    if (calculateKernelUncertainty) {
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
        
    
        std::vector<double> kErrValues(nKernelParameters);
        for (unsigned int idx = 0; idx < nKernelParameters; idx++) {
            // Insanity checking
            if (std::isnan( Error2(idx, idx) )) {
                throw LSST_EXCEPT(exceptions::Exception, 
                                  str(boost::format("Unable to determine kernel uncertainty %d (nan)") % idx));
            }
            if (Error2(idx, idx) < 0.0) {
                throw LSST_EXCEPT(exceptions::Exception,
                                  str(boost::format("Unable to determine kernel uncertainty, negative variance %d (%.3e)") % 
                                      idx % Error2(idx, idx)));
            }
            kErrValues[idx] = sqrt(Error2(idx, idx));
        }
        _kernelError.reset( new math::LinearCombinationKernel(_basisList, kErrValues) );
    
        // Estimate of Background and Background Error */
        if (std::isnan( Error2(nParameters-1, nParameters-1) )) {
            throw LSST_EXCEPT(exceptions::Exception, "Unable to determine background uncertainty (nan)");
        }
        if (Error2(nParameters-1, nParameters-1) < 0.0) {
            throw LSST_EXCEPT(exceptions::Exception, 
                              str(boost::format("Unable to determine background uncertainty, negative variance (%.3e)") % 
                                  Error2(nParameters-1, nParameters-1) 
                                  ));
        }
        _backgroundError = sqrt(Error2(nParameters-1, nParameters-1));
    }

    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                       "Total compute time to do matrix math : %.2f s", time);
    
}


template <typename ImageT>
Eigen::MatrixXd diffim::imageToEigenMatrix(
    lsst::afw::image::Image<ImageT> const& img
    ) {
    unsigned int rows = img.getHeight();
    unsigned int cols = img.getWidth();
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(rows, cols);

    for (int y = 0; y != img.getHeight(); ++y) {
        int x = 0;
        for (typename lsst::afw::image::Image<ImageT>::x_iterator ptr = img.row_begin(y); ptr != img.row_end(y); ++ptr, ++x) {
            M(x,y) = *ptr;
        }
    }
    return M;
}
    

template <typename ImageT, typename VarT>
void diffim::PsfMatchingFunctor<ImageT, VarT>::apply2(
    lsst::afw::image::Image<ImageT> const& imageToConvolve,    ///< Image to apply kernel to
    lsst::afw::image::Image<ImageT> const& imageToNotConvolve, ///< Image whose PSF you want to match to
    lsst::afw::image::Image<VarT>   const& varianceEstimate,   ///< Estimate of the variance per pixel
    lsst::pex::policy::Policy       const& policy              ///< Policy file
    ) {
    
    unsigned int const nKernelParameters     = _basisList.size();
    unsigned int const nBackgroundParameters = 1;
    unsigned int const nParameters           = nKernelParameters + nBackgroundParameters;
    std::vector<boost::shared_ptr<math::Kernel> >::const_iterator kiter = _basisList.begin();
    
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
    unsigned int const startCol = (*kiter)->getCtrX();
    unsigned int const startRow = (*kiter)->getCtrY();
    unsigned int const endCol   = imageToConvolve.getWidth()  - ((*kiter)->getWidth()  - (*kiter)->getCtrX()) + 1;
    unsigned int const endRow   = imageToConvolve.getHeight() - ((*kiter)->getHeight() - (*kiter)->getCtrY()) + 1;
    
    boost::timer t;
    t.restart();
    
    /* least squares matrices */
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(nParameters, nParameters);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(nParameters);
    /* eigen representation of input images; only the pixels that are unconvolved in cimage below */
    Eigen::MatrixXd eigenToConvolveM    = diffim::imageToEigenMatrix(imageToConvolve).block(startRow, startCol, endRow-startRow, endCol-startCol);
    Eigen::MatrixXd eigenToNotConvolveM = diffim::imageToEigenMatrix(imageToNotConvolve).block(startRow, startCol, endRow-startRow, endCol-startCol);
    Eigen::MatrixXd eigeniVarianceM     = diffim::imageToEigenMatrix(varianceEstimate).block(startRow, startCol, endRow-startRow, endCol-startCol).cwise().inverse();
    /* 
       Finally, turn them into vectors for quick matrix updating.  We can technically do everything we want below using Matrices, 
       since we are doing cwise() coefficient-wise multiplication, but doing this as Vectors seems to provide a slight speed up 
       in the multiplication stages, 5% or so, at the expense of a bit more complexity in the code.
    */
    eigenToConvolveM.resize(eigenToConvolveM.rows()*eigenToConvolveM.cols(), 1);
    eigenToNotConvolveM.resize(eigenToNotConvolveM.rows()*eigenToNotConvolveM.cols(), 1);
    eigeniVarianceM.resize(eigeniVarianceM.rows()*eigeniVarianceM.cols(), 1);
    Eigen::VectorXd eigenToConvolveV     = eigenToConvolveM.col(0);
    Eigen::VectorXd eigenToNotConvolveV  = eigenToNotConvolveM.col(0);
    Eigen::VectorXd eigeniVarianceV      = eigeniVarianceM.col(0);
    
    /* holds image convolved with basis function */
    image::Image<ImageT> cimage(imageToConvolve.getDimensions());
    
    /* holds eigen representation of image convolved with all basis functions */
    std::vector<boost::shared_ptr<Eigen::VectorXd> > convolvedEigenList(nKernelParameters);
    
    /* iterators over convolved image list and basis list */
    typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiter = convolvedEigenList.begin();
    /* create C_i in the formalism of Alard & Lupton */
    for (; kiter != _basisList.end(); ++kiter, ++eiter) {
        math::convolve(cimage, imageToConvolve, **kiter, false); /* cimage stores convolved image */
        Eigen::MatrixXd cmat = diffim::imageToEigenMatrix(cimage).block(startRow, startCol, endRow-startRow, endCol-startCol);
	cmat.resize(cmat.rows()*cmat.cols(), 1);
	boost::shared_ptr<Eigen::VectorXd> vmat (new Eigen::VectorXd(cmat.col(0)));
	*eiter = vmat;
    } 
    
    double time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                       "Total compute time to do basis convolutions : %.2f s", time);
    t.restart();
    
    typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiteri = convolvedEigenList.begin();
    typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiterE = convolvedEigenList.end();
    for (unsigned int kidxi = 0; eiteri != eiterE; eiteri++, kidxi++) {
        
        typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiterj = eiteri;
        for (unsigned int kidxj = kidxi; eiterj != eiterE; eiterj++, kidxj++) {
            M(kidxi, kidxj) = ( ( (*eiteri)->cwise() * (**eiterj)).cwise() * eigeniVarianceV).sum();
            /* Equivalent to :
               Eigen::VectorXd mij = (*eiteri)->cwise() * (**eiterj);
               mij.cwise()        *= eigeniVarianceV;
               M(kidxi, kidxj)     = mij.sum();
            */
        }
	B(kidxi)                 = ( ( (eigenToNotConvolveV.cwise()) * (**eiteri)).cwise() * eigeniVarianceV).sum();
	M(kidxi, nParameters-1)  = ((*eiteri)->cwise() * eigeniVarianceV).sum();
    }
    /* background term */
    B(nParameters-1)                = (eigenToNotConvolveV.cwise() * eigeniVarianceV).sum();
    M(nParameters-1, nParameters-1) = eigeniVarianceV.sum();

    // Fill in rest of M
    for (unsigned int kidxi=0; kidxi < nParameters; ++kidxi) {
        for (unsigned int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) {
            M(kidxj, kidxi) = M(kidxi, kidxj);
        }
    }
    
    std::cout << "M2 " << std::endl;
    std::cout << M << std::endl;
    std::cout << "B2 " << std::endl;
    std::cout << B << std::endl;

    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                       "Total compute time to step through pixels : %.2f s", time);
    t.restart();

    /* If the regularization matrix is here and not null, we use it by default */
    if ( (_H.rows() != 0) && (_H.cols() != 0) ) {
        double regularizationScaling = policy.getDouble("regularizationScaling");        
        /* 
           See N.R. 18.5 equation 18.5.8 for the solution to the regularized
           normal equations.  For M.x = b, and solving for x, 

           M -> (Mt.M + lambda*H)
           B -> (Mt.B)

           An estimate of lambda is NR 18.5.16

           lambda = Trace(Mt.M) / Tr(H)

         */

        Eigen::MatrixXd Mt = M.transpose();
        M = Mt * M;

        double lambda = M.trace() / _H.trace();
        lambda *= regularizationScaling;

        M = M + lambda * _H;
        B = Mt * B;
        logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                           "Applying kernel regularization with lambda = %.2e", lambda);
    }
    

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

    /* Fill in the kernel results */
    std::vector<double> kValues(nKernelParameters);
    for (unsigned int idx = 0; idx < nKernelParameters; idx++) {
        if (std::isnan( Soln(idx) )) {
            throw LSST_EXCEPT(exceptions::Exception, 
                              str(boost::format("Unable to determine kernel solution %d (nan)") % idx));
        }
        kValues[idx] = Soln(idx);
    }
    _kernel.reset( new math::LinearCombinationKernel(_basisList, kValues) );

    if (std::isnan( Soln(nParameters-1) )) {
        throw LSST_EXCEPT(exceptions::Exception, 
                          str(boost::format("Unable to determine background solution %d (nan)") % (nParameters-1)));
    }
    _background = Soln(nParameters-1);


    /* Calculate the kernel & background uncertainties only if requested */
    bool calculateKernelUncertainty = policy.getBool("calculateKernelUncertainty");            
    if (calculateKernelUncertainty) {
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
        
    
        std::vector<double> kErrValues(nKernelParameters);
        for (unsigned int idx = 0; idx < nKernelParameters; idx++) {
            // Insanity checking
            if (std::isnan( Error2(idx, idx) )) {
                throw LSST_EXCEPT(exceptions::Exception, 
                                  str(boost::format("Unable to determine kernel uncertainty %d (nan)") % idx));
            }
            if (Error2(idx, idx) < 0.0) {
                throw LSST_EXCEPT(exceptions::Exception,
                                  str(boost::format("Unable to determine kernel uncertainty, negative variance %d (%.3e)") % 
                                      idx % Error2(idx, idx)));
            }
            kErrValues[idx] = sqrt(Error2(idx, idx));
        }
        _kernelError.reset( new math::LinearCombinationKernel(_basisList, kErrValues) );
    
        // Estimate of Background and Background Error */
        if (std::isnan( Error2(nParameters-1, nParameters-1) )) {
            throw LSST_EXCEPT(exceptions::Exception, "Unable to determine background uncertainty (nan)");
        }
        if (Error2(nParameters-1, nParameters-1) < 0.0) {
            throw LSST_EXCEPT(exceptions::Exception, 
                              str(boost::format("Unable to determine background uncertainty, negative variance (%.3e)") % 
                                  Error2(nParameters-1, nParameters-1) 
                                  ));
        }
        _backgroundError = sqrt(Error2(nParameters-1, nParameters-1));
    }

    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                       "Total compute time to do matrix math : %.2f s", time);
    

}

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

Eigen::MatrixXd
diffim::generateFiniteDifferenceRegularization(
    unsigned int width,
    unsigned int height,
    unsigned int order,
    unsigned int boundary_style,  // 0 = unwrapped, 1 = wrapped, 2 = order-tappered ('order' is highest used)
    unsigned int difference_style, // 0 = forward, 1 = central
    bool printB // a debug flag ... remove when done.
					    ) {

    if ( (order < 0) || (order > 2) ) throw LSST_EXCEPT(exceptions::Exception, "Only orders 0..2 allowed");
    if ( (width < 0) )  throw LSST_EXCEPT(exceptions::Exception, "Width < 0");
    if ( (height < 0) ) throw LSST_EXCEPT(exceptions::Exception, "Height < 0");

    if ( (boundary_style < 0) || (boundary_style > 2) ) { 
	throw LSST_EXCEPT(exceptions::Exception, "Boundary styles 0..2 defined");
    }
    if ( (difference_style < 0) || (difference_style > 1) ) {
	throw LSST_EXCEPT(exceptions::Exception, "Only forward (0), and central (1) difference styles defined.");
    }

    /* what works, and what doesn't */
    // == good job == 
    // - order 0, wrapped, forward
    // - order 1, wrapped or tapered, central or forward
    // - order 2, wrapped or tapered, forward
    // == bad job (usually diagongal stripes) ==
    // - all others


    /* 
       Instead of Taylor expanding the forward difference approximation of
       derivatives (N.R. section 18.5) lets just hard code in the expansion of
       the 1st through 3rd derivatives, which will try and enforce smoothness of
       0 through 2nd derivatives.

       A property of the basic "finite difference regularization" is that their
       rows (column multipliers) sum to 0.

       Another consideration is to use *multiple* finite difference operators as
       a constraint.

    */


    // ===========================================================================
    // Get the coeffs for the finite differencing
    // note: The coeffs are stored 2D although they are essentially 1D entities.
    //       The 2d design was chosen to allow cross-terms to be included,
    //         though they are not yet implemented.
    //
    std::vector<std::vector<std::vector<float> > > 
	coeffs(3, std::vector<std::vector<float> >(5, std::vector<float>(5,0)));
    unsigned int x_cen = 0,  y_cen = 0;  // center of reqested order coeffs
    unsigned int x_cen1 = 0, y_cen1 = 0; // center of order 1 coeffs
    unsigned int x_cen2 = 0, y_cen2 = 0; // center of order 2 coeffs
    unsigned int x_size = 0, y_size = 0;

    // forward difference coefficients
    if ( difference_style == 0 ) {
	
	y_cen  = x_cen  = 0;
	x_cen1 = y_cen1 = 0;
	x_cen2 = y_cen2 = 0;

	x_size = y_size = order + 2;

	// default forward difference suggested in NR chap 18
	// 0th order
	coeffs[0][0][0] = -2; coeffs[0][0][1] = 1; 
	coeffs[0][1][0] = 1;  coeffs[0][1][1] = 0;

	// 1st 2
	coeffs[1][0][0] = -2; coeffs[1][0][1] = 2;  coeffs[1][0][2] = -1; 
	coeffs[1][1][0] = 2;  coeffs[1][1][1] = 0;  coeffs[1][1][2] =  0; 
	coeffs[1][2][0] = -1; coeffs[1][2][1] = 0;  coeffs[1][2][2] =  0; 

	// 2nd 2
	coeffs[2][0][0] = -2; coeffs[2][0][1] = 3;  coeffs[2][0][2] = -3; coeffs[2][0][3] = 1; 
	coeffs[2][1][0] = 3;  coeffs[2][1][1] = 0;  coeffs[2][1][2] =  0; coeffs[2][1][3] = 0; 
	coeffs[2][2][0] = -3; coeffs[2][2][1] = 0;  coeffs[2][2][2] =  0; coeffs[2][2][3] = 0; 
	coeffs[2][3][0] = 1;  coeffs[2][3][1] = 0;  coeffs[2][3][2] =  0; coeffs[2][3][3] = 0; 

    }

    // central difference coefficients
    if ( difference_style == 1 ) {

	// this is asymmetric and produces diagonal banding in the kernel
	// from: http://www.holoborodko.com/pavel/?page_id=239
	if ( order == 0 ) { 
	    y_cen = x_cen = 1;
	    x_size = y_size = 3;
	}
	coeffs[0][0][0] =  0; coeffs[0][0][1] = -1;  coeffs[0][0][2] =  0; 
	coeffs[0][1][0] = -1; coeffs[0][1][1] =  0;  coeffs[0][1][2] =  1; 
	coeffs[0][2][0] =  0; coeffs[0][2][1] =  1;  coeffs[0][2][2] =  0; 

	// this works well and is largely the same as order=1 forward-diff.
	// from: http://www.holoborodko.com/pavel/?page_id=239
	if ( order == 1 ) { 
	    y_cen = x_cen = 1;
	    x_size = y_size = 3;
	}
	y_cen1 = x_cen1 = 1;
	coeffs[1][0][0] =  0; coeffs[1][0][1] =  1;  coeffs[1][0][2] =  0;  
	coeffs[1][1][0] =  1; coeffs[1][1][1] = -4;  coeffs[1][1][2] =  1; 
	coeffs[1][2][0] =  0; coeffs[1][2][1] =  1;  coeffs[1][2][2] =  0;  

	// asymmetric and produces diagonal banding in the kernel
	// from http://www.holoborodko.com/pavel/?page_id=239
	if ( order == 2 ) { 
	    y_cen = x_cen = 2;
	    x_size = y_size = 5;
	}
	y_cen2 = x_cen2 = 2;
	coeffs[2][0][0] =  0; coeffs[2][0][1] =  0;  coeffs[2][0][2] = -1; coeffs[2][0][3] =  0; coeffs[2][0][4] =  0; 
	coeffs[2][1][0] =  0; coeffs[2][1][1] =  0;  coeffs[2][1][2] =  2; coeffs[2][1][3] =  0; coeffs[2][1][4] =  0; 
	coeffs[2][2][0] = -1; coeffs[2][2][1] =  2;  coeffs[2][2][2] =  0; coeffs[2][2][3] = -2; coeffs[2][2][4] =  1; 
	coeffs[2][3][0] =  0; coeffs[2][3][1] =  0;  coeffs[2][3][2] = -2; coeffs[2][3][3] =  0; coeffs[2][3][4] =  0; 
	coeffs[2][4][0] =  0; coeffs[2][4][1] =  0;  coeffs[2][4][2] =  1; coeffs[2][4][3] =  0; coeffs[2][4][4] =  0; 
	
    }


    /* Note we have to add 1 extra (empty) term here because of the differential
     * background fitting */
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(width*height+1, width*height+1);

    /* Forward difference approximation */
    for (unsigned int i = 0; i < width*height; i++) {

	unsigned int const x0 = i % width;  // the x coord in the kernel image
	unsigned int const y0 = i / width;  // the y coord in the kernel image
	
	unsigned int x_edge_distance = (x0 > (width - x0 - 1))  ? width - x0 - 1  : x0;
	unsigned int y_edge_distance = (y0 > (height - y0 - 1)) ? height - y0 - 1 : y0;
	unsigned int edge_distance = (x_edge_distance < y_edge_distance) ? x_edge_distance : y_edge_distance;

        for (unsigned int dx = 0; dx < x_size; dx++) {
	    for (unsigned int dy = 0; dy < y_size; dy++) {

		// determine where to put this coeff

		// handle the boundary condition
		// note: adding width and height in the sum prevents negatives
		unsigned int x = 0;
		unsigned int y = 0; 
		double this_coeff = 0;

		// no-wrapping at edges
		if ( boundary_style == 0 ) {
		    x = x0 + dx - x_cen;
		    y = y0 + dy - y_cen;
		    if ( y < 0 || y > height - 1 || x < 0 || x > width - 1 ) { continue; }
		    this_coeff = coeffs[order][dx][dy];

		// wrapping at edges
		} else if ( boundary_style == 1) {
		    x = (width  + x0 + dx - x_cen) % width;
		    y = (height + y0 + dy - y_cen) % height;
		    this_coeff = coeffs[order][dx][dy];

		// order tapering to the edge (just clone wrapping for now)
		// - use the lowest order possible
		} else if ( boundary_style == 2) {

		    // edge rows and columns ... set to constant
		    if ( edge_distance == 0 ) {
			x = x0;
			y = y0;
			this_coeff = 1;
		    }
		    // in one from edge, use 1st order
		    else if ( edge_distance == 1 && order > 0 ) {
			x = (width  + x0 + dx - x_cen1) % width;
			y = (height + y0 + dy - y_cen1) % height;
			if ( (dx < 3) && (dy < 3) ) { this_coeff = coeffs[1][dx][dy]; } 
		    }
		    // in two from edge, use 2st order if order > 1
		    else if ( edge_distance == 2 && order > 1 ){
			x = (width  + x0 + dx - x_cen2) % width;
			y = (height + y0 + dy - y_cen2) % height;
			if ( (dx < 5) && (dy < 5) ) { this_coeff = coeffs[2][dx][dy]; } 
		    } 
		    // if we're somewhere in the middle
		    else if ( edge_distance > order ) {
			x = (width  + x0 + dx - x_cen) % width;
			y = (height + y0 + dy - y_cen) % height;
		    	this_coeff = coeffs[order][dx][dy];
		    }

		} 

		B(i, y*width + x) = this_coeff;
		
	    }

        }

    }

    if ( printB )  {
	std::cout << B << std::endl;
    }
    
    Eigen::MatrixXd H = B.transpose() * B;
    return H;
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
    unsigned int halfWidth,                ///< size is 2*N + 1
    unsigned int nGauss,                   ///< number of gaussians
    std::vector<double> const& sigGauss,   ///< width of the gaussians
    std::vector<int>    const& degGauss    ///< local spatial variation of gaussians
    ) {
    typedef lsst::afw::math::Kernel::PixelT PixelT;
    typedef image::Image<double> ImageT;

    if (halfWidth < 1) {
        throw LSST_EXCEPT(exceptions::Exception, "halfWidth must be positive");
    }
    if (nGauss != sigGauss.size()) {
        throw LSST_EXCEPT(exceptions::Exception, "sigGauss does not have enough entries");
    }
    if (nGauss != degGauss.size()) {
        throw LSST_EXCEPT(exceptions::Exception, "degGauss does not have enough entries");
    }
    int fullWidth = 2 * halfWidth + 1;
    ImageT image0(fullWidth, fullWidth);
    ImageT image(fullWidth, fullWidth);
    
    math::KernelList<math::Kernel> kernelBasisList;
    for (unsigned int i = 0; i < nGauss; i++) {
        /* 
           sigma = FWHM / ( 2 * sqrt(2 * ln(2)) )
        */
        double sig = sigGauss[i];
        unsigned int deg  = degGauss[i];
        //unsigned int nPar = (deg + 1) * (deg + 2) / 2;

        math::GaussianFunction2<PixelT> gaussian(sig, sig);
        math::AnalyticKernel kernel(fullWidth, fullWidth, gaussian);
        math::PolynomialFunction2<PixelT> polynomial(deg);

        /* 

        We want all the bases except for the first to sum to 0.0.  This allows
        us to achieve kernel flux conservation (Ksum) across the image since all
        the power will be in the first term, which will not vary spatially.

           K(x,y) = Ksum * B_0 + Sum_i : a(x,y) * B_i

        To do this, normalize all Kernels to sum = 1. and subtract B_0 from all
        subsequent kenrels.  

        To get an idea of the relative contribution of each of these basis
        functions later on down the line, lets also normalize them such that 

           Sum(B_i)  == 0.0   *and*
           B_i * B_i == 1.0

        For completeness 

           Sum(B_0)  == 1.0
           B_0 * B_0 != 1.0

        */


        for (unsigned int j = 0, n = 0; j <= deg; j++) {
            for (unsigned int k = 0; k <= (deg - j); k++, n++) {
                /* gaussian to be modified by this term in the polynomial */
                polynomial.setParameter(n, 1.);

                if ( (n == 0) && (i == 0) ) {
                    /* Very first kernel */
                    (void)kernel.computeImage(image0, true);                    
                    boost::shared_ptr<math::Kernel> 
                        kernelPtr( new math::FixedKernel(image0) );
                    kernelBasisList.push_back(kernelPtr);
                    polynomial.setParameter(n, 0.);
                    continue;
                }

                (void)kernel.computeImage(image, false); /* no need to normalize as its done below */
                double ksum = 0.;
                for (int y = 0, v = -halfWidth; y < image.getHeight(); y++, v++) {
                    int u = -halfWidth;
                    for (ImageT::xy_locator ptr = image.xy_at(0, y), end = image.xy_at(image.getWidth(), y); ptr != end; ++ptr.x(), u++) {
                        /* Evaluate from -1 to 1 */
                        *ptr  = *ptr * polynomial(u/static_cast<double>(halfWidth), v/static_cast<double>(halfWidth));
                        ksum += *ptr;
                    }
                }

                /* Normalize resulting basis to have sum = 1 */
                image /= ksum;
                /* All bases except the very first have zero total flux */
                image -= image0;

                /* Lets divide by a factor to make its dot product with itself 1.0 */
                double ksum2 = 0.;
                for (int y = 0; y < image.getHeight(); y++) {
                    for (ImageT::xy_locator ptr = image.xy_at(0, y), end = image.xy_at(image.getWidth(), y); ptr != end; ++ptr.x()) {
                        ksum2 += *ptr * *ptr;
                    }
                }
                image /= sqrt(ksum2);

                boost::shared_ptr<math::Kernel> 
                    kernelPtr( new math::FixedKernel(image) );
                kernelBasisList.push_back(kernelPtr);
                polynomial.setParameter(n, 0.);
            }
        }
    }
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
    lsst::afw::image::MaskedImage<ImageT> const& imageToConvolve,    ///< Image T to convolve with Kernel
    lsst::afw::image::MaskedImage<ImageT> const& imageToNotConvolve, ///< Image I to subtract convolved template from
    lsst::afw::math::Kernel const& convolutionKernel,                ///< PSF-matching Kernel used for convolution
    BackgroundT background,                               ///< Differential background function or scalar
    bool invert                                           ///< Invert the output difference image
    ) {

    logging::TTrace<8>("lsst.ip.diffim.convolveAndSubtract", "Convolving using convolve");
    
    image::MaskedImage<ImageT> convolvedMaskedImage(imageToConvolve.getDimensions());
    convolvedMaskedImage.setXY0(imageToConvolve.getXY0());
    math::convolve(convolvedMaskedImage, imageToConvolve, convolutionKernel, false);
    
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
    lsst::afw::image::Image<ImageT> const& imageToConvolve,          ///< Image T to convolve with Kernel
    lsst::afw::image::MaskedImage<ImageT> const& imageToNotConvolve, ///< Image I to subtract convolved template from
    lsst::afw::math::Kernel const& convolutionKernel,                ///< PSF-matching Kernel used for convolution
    BackgroundT background,                                          ///< Differential background function or scalar
    bool invert                                                      ///< Invert the output difference image
    ) {
    
    logging::TTrace<8>("lsst.ip.diffim.convolveAndSubtract", "Convolving using convolve");
    
    image::MaskedImage<ImageT> convolvedMaskedImage(imageToConvolve.getDimensions());
    convolvedMaskedImage.setXY0(imageToConvolve.getXY0());
    
    math::convolve(*convolvedMaskedImage.getImage(), imageToConvolve, convolutionKernel, false);
    
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
    lsst::afw::image::MaskedImage<ImageT> const& imageToConvolve,    
    lsst::afw::image::MaskedImage<ImageT> const& imageToNotConvolve, 
    lsst::pex::policy::Policy  const& policy                                       
    ) {
    
    // Parse the Policy
    unsigned int fpNpixMin      = policy.getInt("fpNpixMin");
    unsigned int fpNpixMax      = policy.getInt("fpNpixMax");

    int const kCols             = policy.getInt("kernelCols");
    int const kRows             = policy.getInt("kernelRows");
    double fpGrowKsize          = policy.getDouble("fpGrowKsize");

    int minCleanFp              = policy.getInt("minCleanFp");
    double detThreshold         = policy.getDouble("detThreshold");
    double detThresholdScaling  = policy.getDouble("detThresholdScaling");
    double detThresholdMin      = policy.getDouble("detThresholdMin");
    std::string detThresholdType = policy.getString("detThresholdType");

    // Number of pixels to grow each Footprint, based upon the Kernel size
    int fpGrowPix = int(fpGrowKsize * ( (kCols > kRows) ? kCols : kRows ));

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
                logging::TTrace<6>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
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
                logging::TTrace<5>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching",
                                   "Exception caught extracting Footprint");
                logging::TTrace<6>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching",
                                   e.what());
                continue;
            }

            // Search for bad pixels within the footprint
            itcFunctor.apply(*fpGrow);
            if (itcFunctor.getBits() > 0) {
                logging::TTrace<6>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                               "Footprint has bad pix in image to convolve"); 
                continue;
            }

            itncFunctor.apply(*fpGrow);
            if (itncFunctor.getBits() > 0) {
                logging::TTrace<6>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
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

