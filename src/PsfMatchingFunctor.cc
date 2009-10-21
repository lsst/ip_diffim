// -*- lsst-c++ -*-
/**
 * @file PsfMatchingFunctor.cc
 *
 * @brief Implementation of image subtraction functions declared in PsfMatchingFunctor.h
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */
#include <cmath>

#include "boost/timer.hpp" 

#include "Eigen/Core"
#include "Eigen/QR"
#include "Eigen/LU"

#include "lsst/pex/exceptions/Exception.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math.h"
#include "lsst/ip/diffim/PsfMatchingFunctor.h"
#include "lsst/ip/diffim/ImageSubtract.h"

#define DEBUG_MATRIX 0

namespace pexExcept  = lsst::pex::exceptions; 
namespace pexLog     = lsst::pex::logging; 
namespace pexPolicy  = lsst::pex::policy; 
namespace afwImage   = lsst::afw::image;
namespace afwMath    = lsst::afw::math;

namespace lsst { 
namespace ip { 
namespace diffim {

//
// Constructors
//
template <typename PixelT, typename VarT>
PsfMatchingFunctor<PixelT, VarT>::PsfMatchingFunctor(
    lsst::afw::math::KernelList const &basisList
    ) :
    _basisList(basisList),
    _mMat(),
    _bVec(),
    _sVec(),
    _hMat(),
    _initialized(false),
    _regularize(false)
{;}

template <typename PixelT, typename VarT>
PsfMatchingFunctor<PixelT, VarT>::PsfMatchingFunctor(
    lsst::afw::math::KernelList const &basisList,
    boost::shared_ptr<Eigen::MatrixXd> const &hMat
    ) :
    _basisList(basisList),
    _mMat(),
    _bVec(),
    _sVec(),
    _hMat(hMat),
    _initialized(false),
    _regularize(true)
{;}

template <typename PixelT, typename VarT>
PsfMatchingFunctor<PixelT, VarT>::PsfMatchingFunctor(
    const PsfMatchingFunctor<PixelT,VarT> &rhs
    ) :
    _basisList(rhs._basisList),
    _mMat(),
    _bVec(),
    _sVec(),
    _hMat(rhs._hMat),
    _initialized(false),
    _regularize(rhs._regularize)
{;}

//
// Public Member Functions
//

template <typename PixelT, typename VarT>
void PsfMatchingFunctor<PixelT, VarT>::apply(
    lsst::afw::image::Image<PixelT> const &imageToConvolve,    ///< Image to apply kernel to
    lsst::afw::image::Image<PixelT> const &imageToNotConvolve, ///< Image whose PSF you want to match to
    lsst::afw::image::Image<VarT>   const &varEstimate,        ///< Estimate of the variance per pixel
    lsst::pex::policy::Policy       const &policy              ///< Policy file
    ) {
    
    unsigned int const nKernelParameters     = _basisList.size();
    unsigned int const nBackgroundParameters = 1;
    unsigned int const nParameters           = nKernelParameters + nBackgroundParameters;
    std::vector<boost::shared_ptr<afwMath::Kernel> >::const_iterator kiter = _basisList.begin();
    
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
    unsigned int const endCol   = imageToConvolve.getWidth()  - 
        ((*kiter)->getWidth()  - (*kiter)->getCtrX()) + 1;
    unsigned int const endRow   = imageToConvolve.getHeight() - 
        ((*kiter)->getHeight() - (*kiter)->getCtrY()) + 1;
    
    boost::timer t;
    t.restart();
    
    /* Least squares matrices */
    Eigen::MatrixXd mMat = Eigen::MatrixXd::Zero(nParameters, nParameters);
    Eigen::VectorXd bVec = Eigen::VectorXd::Zero(nParameters);
    /* Eigen representation of input images; only the pixels that are unconvolved in cimage below */
    Eigen::MatrixXd eigenToConvolve = imageToEigenMatrix(imageToConvolve).block(startRow, startCol, 
                                                                                endRow-startRow, 
                                                                                endCol-startCol);
    Eigen::MatrixXd eigenToNotConvolve = imageToEigenMatrix(imageToNotConvolve).block(startRow, startCol, 
                                                                                      endRow-startRow, 
                                                                                      endCol-startCol);
    Eigen::MatrixXd eigeniVariance = imageToEigenMatrix(varEstimate).block(startRow, startCol, 
                                                                           endRow-startRow, 
                                                                           endCol-startCol).cwise().inverse();
    /* Resize into 1-D for later usage */
    eigenToConvolve.resize(eigenToConvolve.rows()*eigenToConvolve.cols(), 1);
    eigenToNotConvolve.resize(eigenToNotConvolve.rows()*eigenToNotConvolve.cols(), 1);
    eigeniVariance.resize(eigeniVariance.rows()*eigeniVariance.cols(), 1);

    /* Holds image convolved with basis function */
    afwImage::Image<PixelT> cimage(imageToConvolve.getDimensions());
    
    /* Holds eigen representation of image convolved with all basis functions */
    std::vector<boost::shared_ptr<Eigen::MatrixXd> > convolvedEigenList(nKernelParameters);
    
    /* Iterators over convolved image list and basis list */
    typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiter = convolvedEigenList.begin();
    /* Create C_i in the formalism of Alard & Lupton */
    for (; kiter != _basisList.end(); ++kiter, ++eiter) {
        afwMath::convolve(cimage, imageToConvolve, **kiter, false); /* cimage stores convolved image */
        boost::shared_ptr<Eigen::MatrixXd> cMat (
            new Eigen::MatrixXd(imageToEigenMatrix(cimage).block(startRow, 
                                                                 startCol, 
                                                                 endRow-startRow, 
                                                                 endCol-startCol))
            );
        cMat->resize(cMat->rows()*cMat->cols(), 1);
        *eiter = cMat;
    } 
    
    double time = t.elapsed();
    pexLog::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                      "Total compute time to do basis convolutions : %.2f s", time);
    t.restart();
    
    /* 
     * 
     * NOTE - 
     * 
     * Below is the original Eigen representation of the matrix math needed.
     * Its a bit more readable but 5-10% slower than the as-implemented Eigen
     * math.  Left here for reference as it nicely and simply outlines the math
     * that goes into the construction of M and B.
     * 

    typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiteri = convolvedEigenList.begin();
    typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiterE = convolvedEigenList.end();
    for (unsigned int kidxi = 0; eiteri != eiterE; eiteri++, kidxi++) {
        Eigen::VectorXd eiteriDotiVariance = (*eiteri)->cwise() * eigeniVarianceV;

        typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiterj = eiteri;
        for (unsigned int kidxj = kidxi; eiterj != eiterE; eiterj++, kidxj++) {
            M(kidxi, kidxj) = (eiteriDotiVariance.cwise() * (**eiterj)).sum();
            M(kidxj, kidxi) = M(kidxi, kidxj);
        }
        B(kidxi)                 = (eiteriDotiVariance.cwise() * eigenToNotConvolveV).sum();
        M(kidxi, nParameters-1)  = eiteriDotiVariance.sum();
        M(nParameters-1, kidxi)  = M(kidxi, nParameters-1);
    }
    B(nParameters-1)                = (eigenToNotConvolveV.cwise() * eigeniVarianceV).sum();
    M(nParameters-1, nParameters-1) = eigeniVarianceV.sum();

    */
    
    /* 
       Load matrix with all values from convolvedEigenList : all images
       (eigeniVariance, convolvedEigenList) must be the same size
    */
    Eigen::MatrixXd cMat(eigeniVariance.col(0).size(), nParameters);
    typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiterj = convolvedEigenList.begin();
    typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiterE = convolvedEigenList.end();
    for (unsigned int kidxj = 0; eiterj != eiterE; eiterj++, kidxj++) {
        cMat.col(kidxj) = (*eiterj)->col(0);
    }
    /* Treat the last "image" as all 1's to do the background calculation. */
    cMat.col(nParameters-1).fill(1.);
    
    /* Caculate the variance-weighted pixel values */
    Eigen::MatrixXd vcMat = eigeniVariance.col(0).asDiagonal() * cMat;
    
    /* Calculate M as the variance-weighted inner product of C */
    mMat = cMat.transpose() * vcMat;
    bVec = vcMat.transpose() * eigenToNotConvolve.col(0);

    if (DEBUG_MATRIX) {
        std::cout << "M " << std::endl;
        std::cout << mMat << std::endl;
        std::cout << "B " << std::endl;
        std::cout << bVec << std::endl;
    }

    time = t.elapsed();
    pexLog::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                      "Total compute time to step through pixels : %.2f s", time);
    t.restart();

    /* If the regularization matrix is here and not null, we use it by default */
    if (_regularize) {
        double regularizationScaling = policy.getDouble("regularizationScaling");        
        /* 
           See N.R. 18.5 equation 18.5.8 for the solution to the regularized
           normal equations.  For M.x = b, and solving for x, 

           M -> (Mt.M + lambda*H)
           B -> (Mt.B)

           An estimate of lambda is NR 18.5.16

           lambda = Trace(Mt.M) / Tr(H)

         */

        Eigen::MatrixXd mtMat = mMat.transpose();
        mMat = mtMat * mMat;

        double lambda = mMat.trace() / _hMat->trace();
        lambda *= regularizationScaling;

        mMat = mMat + lambda * *_hMat;
        bVec = mtMat * bVec;
        pexLog::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                          "Applying kernel regularization with lambda = %.2e", lambda);
        
    }
    

    // To use Cholesky decomposition, the matrix needs to be symmetric (M is, by
    // design) and positive definite.  
    //
    // Eventually put a check in here to make sure its positive definite
    //
    Eigen::VectorXd sVec = Eigen::VectorXd::Zero(nParameters);
    if (!(mMat.ldlt().solve(bVec, &sVec))) {
        pexLog::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                          "Unable to determine kernel via Cholesky LDL^T");
        if (!(mMat.llt().solve(bVec, &sVec))) {
            pexLog::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                              "Unable to determine kernel via Cholesky LL^T");
            if (!(mMat.lu().solve(bVec, &sVec))) {
                pexLog::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                                  "Unable to determine kernel via LU");
                // LAST RESORT
                try {
                    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eVecValues(mMat);
                    Eigen::MatrixXd const& rMat = eVecValues.eigenvectors();
                    Eigen::VectorXd eValues = eVecValues.eigenvalues();
                    
                    for (int i = 0; i != eValues.rows(); ++i) {
                        if (eValues(i) != 0.0) {
                            eValues(i) = 1.0/eValues(i);
                        }
                    }
                    
                    sVec = rMat*eValues.asDiagonal()*rMat.transpose()*bVec;
                } catch (pexExcept::Exception& e) {
                    pexLog::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                                      "Unable to determine kernel via eigen-values");
                    
                    throw LSST_EXCEPT(pexExcept::Exception, "Unable to determine kernel solution");
                }
            }
        }
    }

    /* Save matrices as they are expensive to calculate.
     * 
     * NOTE : one might consider saving the VC and B vectors instead of M and B;
     * however then we would not be able to maintain the regularization of M
     * even though the stored B would be regularized.
     * 
     * ANOTHER NOTE : we might also consider *not* solving for Soln here, in the
     * case that we don't care about the results of the single-kernel fit.  That
     * is, if we decide to only do sigma clipping on the spatial results.
     * 
     */
    _mMat = boost::shared_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd(mMat));
    _bVec = boost::shared_ptr<Eigen::VectorXd>(new Eigen::VectorXd(bVec));
    _sVec = boost::shared_ptr<Eigen::VectorXd>(new Eigen::VectorXd(sVec));
    _initialized = true;
    time = t.elapsed();
    pexLog::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                      "Total compute time to do matrix math : %.2f s", time);
    
}

template <typename PixelT, typename VarT>
std::pair<boost::shared_ptr<lsst::afw::math::Kernel>, double>
PsfMatchingFunctor<PixelT, VarT>::getSolution() {

    if (!(_initialized)) {
        throw LSST_EXCEPT(pexExcept::Exception, "Kernel not initialized");
    }

    unsigned int const nKernelParameters     = _basisList.size();
    unsigned int const nBackgroundParameters = 1;
    unsigned int const nParameters           = nKernelParameters + nBackgroundParameters;

    /* Fill in the kernel results */
    std::vector<double> kValues(nKernelParameters);
    for (unsigned int idx = 0; idx < nKernelParameters; idx++) {
        if (std::isnan((*_sVec)(idx))) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              str(boost::format("Unable to determine kernel solution %d (nan)") % idx));
        }
        kValues[idx] = (*_sVec)(idx);
    }
    boost::shared_ptr<afwMath::Kernel> kernel( 
        new afwMath::LinearCombinationKernel(_basisList, kValues) 
        );
    
    if (std::isnan((*_sVec)(nParameters-1))) {
        throw LSST_EXCEPT(pexExcept::Exception, 
                          str(boost::format("Unable to determine background solution %d (nan)") % (nParameters-1)));
    }
    double background = (*_sVec)(nParameters-1);

    return std::make_pair(kernel, background);
}

template <typename PixelT, typename VarT>
std::pair<boost::shared_ptr<lsst::afw::math::Kernel>, double>
PsfMatchingFunctor<PixelT, VarT>::getSolutionUncertainty() {

    if (!(_initialized)) {
        throw LSST_EXCEPT(pexExcept::Exception, "Kernel not initialized");
    }

    unsigned int const nKernelParameters     = _basisList.size();
    unsigned int const nBackgroundParameters = 1;
    unsigned int const nParameters           = nKernelParameters + nBackgroundParameters;

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
    Eigen::MatrixXd             Cov    = (*_mMat).transpose() * (*_mMat);
    Eigen::LLT<Eigen::MatrixXd> llt    = Cov.llt();
    Eigen::MatrixXd             Error2 = llt.matrixL().transpose().inverse() * llt.matrixL().inverse();
        
    std::vector<double> kErrValues(nKernelParameters);
    for (unsigned int idx = 0; idx < nKernelParameters; idx++) {
        // Insanity checking
        if (std::isnan(Error2(idx, idx))) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              str(boost::format("Unable to determine kernel err %d (nan)") % idx));
        }
        if (Error2(idx, idx) < 0.0) {
            throw LSST_EXCEPT(pexExcept::Exception,
                              str(boost::format("Unable to determine kernel err, negative var %d (%.3e)") % 
                                  idx % Error2(idx, idx)));
        }
        kErrValues[idx] = std::sqrt(Error2(idx, idx));
    }
    boost::shared_ptr<afwMath::Kernel> kernelErr( 
        new afwMath::LinearCombinationKernel(_basisList, kErrValues) 
        );
 
    // Estimate of Background and Background Error */
    if (std::isnan(Error2(nParameters-1, nParameters-1))) {
        throw LSST_EXCEPT(pexExcept::Exception, "Unable to determine background err (nan)");
    }
    if (Error2(nParameters-1, nParameters-1) < 0.0) {
        throw LSST_EXCEPT(pexExcept::Exception, 
                          str(boost::format("Unable to determine background err, negative var (%.3e)") % 
                              Error2(nParameters-1, nParameters-1) 
                              ));
    }
    double backgroundErr = std::sqrt(Error2(nParameters-1, nParameters-1));

    return std::make_pair(kernelErr, backgroundErr);
}

template <typename PixelT, typename VarT>
std::pair<boost::shared_ptr<Eigen::MatrixXd>, boost::shared_ptr<Eigen::VectorXd> >
PsfMatchingFunctor<PixelT, VarT>::getAndClearMB() {

    boost::shared_ptr<Eigen::MatrixXd> mOut = _mMat;
    boost::shared_ptr<Eigen::VectorXd> bOut = _bVec;
    _mMat.reset();
    _bVec.reset();
    _sVec.reset();
    _initialized=false;
    return std::make_pair(mOut, bOut);
}

template class PsfMatchingFunctor<float, float>;
template class PsfMatchingFunctor<double, float>;

}}} // end of namespace lsst::ip::diffim
