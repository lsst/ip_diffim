// -*- lsst-c++ -*-
/**
 * @file KernelCandidate.cc
 *
 * @brief Implementation of KernelCandidate class
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include "boost/timer.hpp" 

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/ip/diffim/KernelCandidate.h"
#include "lsst/ip/diffim/ImageSubtract.h"
#include "lsst/ip/diffim/KernelSolution.h"

#define DEBUG_MATRIX 0

namespace afwMath        = lsst::afw::math;
namespace afwImage       = lsst::afw::image;
namespace pexLog         = lsst::pex::logging; 
namespace pexExcept      = lsst::pex::exceptions; 

namespace lsst { 
namespace ip { 
namespace diffim {

    template <typename PixelT>
    KernelCandidate<PixelT>::KernelCandidate(
        float const xCenter,
        float const yCenter, 
        MaskedImagePtr const& miToConvolvePtr,
        MaskedImagePtr const& miToNotConvolvePtr,
        lsst::pex::policy::Policy const& policy
        ) :
        lsst::afw::math::SpatialCellImageCandidate<ImageT>(xCenter, yCenter),
        _miToConvolvePtr(miToConvolvePtr),
        _miToNotConvolvePtr(miToNotConvolvePtr),
        _varianceEstimate(),
        _policy(policy),
        _coreFlux(),
        _isInitialized(false),
        _useRegularization(false),
        _fitForBackground(_policy.getBool("fitForBackground")),
        _kernelSolutionOrig(),
        _kernelSolutionPca()
    {
        
        /* Rank by mean core S/N in science image */
        ImageStatistics<PixelT> imstats;
        int candidateCoreRadius = _policy.getInt("candidateCoreRadius");
        imstats.apply(*_miToNotConvolvePtr, candidateCoreRadius);
        _coreFlux = imstats.getMean();

        pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate",
                          "Candidate %d at %.2f %.2f with ranking %.2f", 
                          this->getId(), this->getXCenter(), this->getYCenter(), _coreFlux);
    }
    

    template <typename PixelT>
    void KernelCandidate<PixelT>::build(
        lsst::afw::math::KernelList const& basisList
        ) {
        build(basisList, boost::shared_ptr<Eigen::MatrixXd>());
    }
    
    
    template <typename PixelT>
    void KernelCandidate<PixelT>::build(
        lsst::afw::math::KernelList const& basisList,
        boost::shared_ptr<Eigen::MatrixXd> hMat
        ) {

        /* Examine the policy for control over the variance estimate */
        afwImage::MaskedImage<PixelT> var = 
            afwImage::MaskedImage<PixelT>(*(_miToNotConvolvePtr), true);
        if (_policy.getBool("constantVarianceWeighting")) {
            /* Constant variance weighting */
            *var.getVariance() = 1.;
        }
        else {
            /* Variance estimate comes from the straight difference */
            var -= *(_miToConvolvePtr);
        }
        _varianceEstimate = var.getVariance();

        /* Do we have a regularization matrix?  If so use it */
        if (hMat) {
            _useRegularization = true;
            pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate.build", 
                              "Using kernel regularization");

            if (_isInitialized) {
                _kernelSolutionPca = boost::shared_ptr<StaticKernelSolution2<PixelT> >(
                    new RegularizedKernelSolution<PixelT>(basisList, _fitForBackground, hMat, _policy)
                    );
                _kernelSolutionPca->build(*(_miToConvolvePtr->getImage()),
                                          *(_miToNotConvolvePtr->getImage()),
                                          *_varianceEstimate);
                _kernelSolutionPca->solve();
            }
            else {
                _kernelSolutionOrig = boost::shared_ptr<StaticKernelSolution2<PixelT> >(
                    new RegularizedKernelSolution<PixelT>(basisList, _fitForBackground, hMat, _policy)
                    );
                _kernelSolutionOrig->build(*(_miToConvolvePtr->getImage()),
                                           *(_miToNotConvolvePtr->getImage()),
                                           *_varianceEstimate);
                _kernelSolutionOrig->solve();
            }
        }
        else {
            _useRegularization = false;
            pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate.build", 
                              "Not using kernel regularization");
            if (_isInitialized) {
                _kernelSolutionPca = boost::shared_ptr<StaticKernelSolution2<PixelT> >(
                    new StaticKernelSolution2<PixelT>(basisList, _fitForBackground)
                    );
                _kernelSolutionPca->build(*(_miToConvolvePtr->getImage()),
                                          *(_miToNotConvolvePtr->getImage()),
                                          *_varianceEstimate);
                _kernelSolutionPca->solve();
            }
            else {
                _kernelSolutionOrig = boost::shared_ptr<StaticKernelSolution2<PixelT> >(
                    new StaticKernelSolution2<PixelT>(basisList, _fitForBackground)
                    );
                _kernelSolutionOrig->build(*(_miToConvolvePtr->getImage()),
                                           *(_miToNotConvolvePtr->getImage()),
                                           *_varianceEstimate);
                _kernelSolutionOrig->solve();
            }
        }

        
        /*
        try {
            buildEngine(basisList, hMat);
        } catch (pexExcept::Exception &e) {
            throw e;
        }
        */

        if (_policy.getBool("iterateSingleKernel") && (!(_policy.getBool("constantVarianceWeighting")))) {
            afwImage::MaskedImage<PixelT> diffim = getDifferenceImage(KernelCandidate::RECENT);
            _varianceEstimate = diffim.getVariance();
            //try {
            //buildEngine(basisList, hMat);
            //} catch (pexExcept::Exception &e) {
            //throw e;
            //}
        }

        /* Don't put this in buildEngine() since you might have to iterate the
         * first kernel */
        _isInitialized = true;

    }
            

//    template <typename PixelT>
//    void KernelCandidate<PixelT>::buildEngine(
//        lsst::afw::math::KernelList const& basisList,
//        boost::shared_ptr<Eigen::MatrixXd> hMat
//        ) {
//
//        lsst::afw::image::Image<PixelT> const &imageToConvolve = *(_miToConvolvePtr->getImage());
//        lsst::afw::image::Image<PixelT> const &imageToNotConvolve = *(_miToNotConvolvePtr->getImage());
//        
//        unsigned int const nKernelParameters     = basisList.size();
//        unsigned int const nBackgroundParameters = _fitForBackground ? 1 : 0;
//        unsigned int const nParameters           = nKernelParameters + nBackgroundParameters;
//        std::vector<boost::shared_ptr<afwMath::Kernel> >::const_iterator kiter = basisList.begin();
//        
//        /* Ignore buffers around edge of convolved images :
//         * 
//         * If the kernel has width 5, it has center pixel 2.  The first good pixel
//         * is the (5-2)=3rd pixel, which is array index 2, and ends up being the
//         * index of the central pixel.
//         * 
//         * You also have a buffer of unusable pixels on the other side, numbered
//         * width-center-1.  The last good usable pixel is N-width+center+1.
//         * 
//         * Example : the kernel is width = 5, center = 2
//         * 
//         * ---|---|-c-|---|---|
//         * 
//         * the image is width = N
//         * convolve this with the kernel, and you get
//         * 
//         * |-x-|-x-|-g-|---|---| ... |---|---|-g-|-x-|-x-|
//         * 
//         * g = first/last good pixel
//         * x = bad
//         * 
//         * the first good pixel is the array index that has the value "center", 2
//         * the last good pixel has array index N-(5-2)+1
//         * eg. if N = 100, you want to use up to index 97
//         * 100-3+1 = 98, and the loops use i < 98, meaning the last
//         * index you address is 97.
//         */
//        unsigned int const startCol = (*kiter)->getCtrX();
//        unsigned int const startRow = (*kiter)->getCtrY();
//        unsigned int const endCol   = imageToConvolve.getWidth()  - 
//            ((*kiter)->getWidth()  - (*kiter)->getCtrX()) + 1;
//        unsigned int const endRow   = imageToConvolve.getHeight() - 
//            ((*kiter)->getHeight() - (*kiter)->getCtrY()) + 1;
//        
//        boost::timer t;
//        t.restart();
//        
//        /* Eigen representation of input images; only the pixels that are unconvolved in cimage below */
//        Eigen::MatrixXd eigenToConvolve = imageToEigenMatrix(imageToConvolve).block(startRow, 
//                                                                                    startCol, 
//                                                                                    endRow-startRow, 
//                                                                                    endCol-startCol);
//        Eigen::MatrixXd eigenToNotConvolve = imageToEigenMatrix(imageToNotConvolve).block(startRow, 
//                                                                                          startCol, 
//                                                                                          endRow-startRow, 
//                                                                                          endCol-startCol);
//        Eigen::MatrixXd eigeniVariance = imageToEigenMatrix(*_varianceEstimate).block(
//            startRow, startCol, endRow-startRow, endCol-startCol).cwise().inverse();
//
//        /* Resize into 1-D for later usage */
//        eigenToConvolve.resize(eigenToConvolve.rows()*eigenToConvolve.cols(), 1);
//        eigenToNotConvolve.resize(eigenToNotConvolve.rows()*eigenToNotConvolve.cols(), 1);
//        eigeniVariance.resize(eigeniVariance.rows()*eigeniVariance.cols(), 1);
//        
//        /* Holds image convolved with basis function */
//        afwImage::Image<PixelT> cimage(imageToConvolve.getDimensions());
//        
//        /* Holds eigen representation of image convolved with all basis functions */
//        std::vector<boost::shared_ptr<Eigen::MatrixXd> > convolvedEigenList(nKernelParameters);
//        
//        /* Iterators over convolved image list and basis list */
//        typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiter = 
//            convolvedEigenList.begin();
//        /* Create C_i in the formalism of Alard & Lupton */
//        for (; kiter != basisList.end(); ++kiter, ++eiter) {
//            afwMath::convolve(cimage, imageToConvolve, **kiter, false); /* cimage stores convolved image */
//            boost::shared_ptr<Eigen::MatrixXd> cMat (
//                new Eigen::MatrixXd(imageToEigenMatrix(cimage).block(startRow, 
//                                                                     startCol, 
//                                                                     endRow-startRow, 
//                                                                     endCol-startCol))
//                );
//            cMat->resize(cMat->rows()*cMat->cols(), 1);
//            *eiter = cMat;
//
//        } 
//
//        double time = t.elapsed();
//        pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate.build", 
//                          "Total compute time to do basis convolutions : %.2f s", time);
//        t.restart();
//        
//        /* 
//         * 
//         * NOTE - 
//         * 
//         * Below is the original Eigen representation of the matrix math needed.
//         * Its a bit more readable but 5-10% slower than the as-implemented Eigen
//         * math.  Left here for reference as it nicely and simply outlines the math
//         * that goes into the construction of M and B.
//         * 
//         
//         typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiteri = 
//             convolvedEigenList.begin();
//         typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiterE = 
//             convolvedEigenList.end();
//         for (unsigned int kidxi = 0; eiteri != eiterE; eiteri++, kidxi++) {
//         Eigen::VectorXd eiteriDotiVariance = (*eiteri)->cwise() * eigeniVarianceV;
//         
//         typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiterj = eiteri;
//         for (unsigned int kidxj = kidxi; eiterj != eiterE; eiterj++, kidxj++) {
//         M(kidxi, kidxj) = (eiteriDotiVariance.cwise() * (**eiterj)).sum();
//         M(kidxj, kidxi) = M(kidxi, kidxj);
//         }
//         B(kidxi)                 = (eiteriDotiVariance.cwise() * eigenToNotConvolveV).sum();
//         M(kidxi, nParameters-1)  = eiteriDotiVariance.sum();
//         M(nParameters-1, kidxi)  = M(kidxi, nParameters-1);
//         }
//         B(nParameters-1)                = (eigenToNotConvolveV.cwise() * eigeniVarianceV).sum();
//         M(nParameters-1, nParameters-1) = eigeniVarianceV.sum();
//         
//        */
//        
//        /* 
//           Load matrix with all values from convolvedEigenList : all images
//           (eigeniVariance, convolvedEigenList) must be the same size
//        */
//        Eigen::MatrixXd cMat(eigeniVariance.col(0).size(), nParameters);
//        typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiterj = 
//            convolvedEigenList.begin();
//        typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiterE = 
//            convolvedEigenList.end();
//        for (unsigned int kidxj = 0; eiterj != eiterE; eiterj++, kidxj++) {
//            cMat.col(kidxj) = (*eiterj)->col(0);
//        }
//        /* Treat the last "image" as all 1's to do the background calculation. */
//        if (_fitForBackground)
//            cMat.col(nParameters-1).fill(1.);
//        
//        /* Caculate the variance-weighted pixel values */
//        Eigen::MatrixXd vcMat = eigeniVariance.col(0).asDiagonal() * cMat;
//        
//
//        /* Calculate M as the variance-weighted inner product of C */
//        boost::shared_ptr<Eigen::MatrixXd> mMat (
//            new Eigen::MatrixXd(cMat.transpose() * vcMat)
//            );
//        boost::shared_ptr<Eigen::VectorXd> bVec (
//            new Eigen::VectorXd(vcMat.transpose() * eigenToNotConvolve.col(0))
//            );
//        
//        if (DEBUG_MATRIX) {
//            std::cout << "M "    << std::endl;
//            std::cout << (*mMat) << std::endl;
//            std::cout << "B "    << std::endl;
//            std::cout << (*bVec) << std::endl;
//        }
//        
//        time = t.elapsed();
//        pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate.build", 
//                          "Total compute time to step through pixels : %.2f s", time);
//        t.restart();
//        
//        /* If the regularization matrix is here and not null, we use it by default */
//        if (_useRegularization) {
//            std::string lambdaType = _policy.getString("lambdaType");        
//            double lambdaValue  = _policy.getDouble("lambdaValue");
//            
//            /* See N.R. 18.5 */
    //
//            /* 
//               We regularize the least squares problem:
//               
//               (M + lambda H) a = b
//
//               Regularizing the normal equations instead would suggest we make
//               a new problem of the form:
//
//               (Mt M + lambda H) a = Mt b
//
//               which yields a slightly different effective noise weighting in
//               the solution.
//
//            */
//            
//            double lambda;
//            if (lambdaType == "absolute") {
//                lambda = lambdaValue;
//            }
//            else if (lambdaType == "relative") {
//                lambda  = mMat->trace() / hMat->trace();
//                lambda *= lambdaValue;
//            }
//            else {
//                throw LSST_EXCEPT(pexExcept::Exception, "lambdaType in Policy not recognized");
//            }
//            (*mMat) += lambda * (*hMat);
//            
//            pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate.build", 
//                              "Applying kernel regularization with lambda = %.2e", lambda);
//            
//        }
//        
//        if (_isInitialized) {
//            _kernelSolutionPca = boost::shared_ptr<StaticKernelSolution>(
//                new StaticKernelSolution(mMat, bVec, _fitForBackground, basisList)
//                );
//            _kernelSolutionPca->solve(false);
//        }
//        else {
//            _kernelSolutionOrig = boost::shared_ptr<StaticKernelSolution>(
//                new StaticKernelSolution(mMat, bVec, _fitForBackground, basisList)
//                );
//            _kernelSolutionOrig->solve(false);
//        }
//    }
    
    
    
    template <typename PixelT>
    lsst::afw::math::Kernel::Ptr KernelCandidate<PixelT>::getKernel(CandidateSwitch cand) const {
        if (cand == KernelCandidate::ORIG) {
            if (_kernelSolutionOrig) 
                return _kernelSolutionOrig->getKernel();
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Orignal kernel does not exist");
        }
        else if (cand == KernelCandidate::PCA) {
            if (_kernelSolutionPca) 
                return _kernelSolutionPca->getKernel();
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Pca kernel does not exist");
        }
        else if (cand == KernelCandidate::RECENT) {
            if (_kernelSolutionPca) 
                return _kernelSolutionPca->getKernel();
            else if (_kernelSolutionOrig)
                return _kernelSolutionOrig->getKernel();
            else
                throw LSST_EXCEPT(pexExcept::Exception, "No kernels exist");
        }
        else {
            throw LSST_EXCEPT(pexExcept::Exception, "Invalid CandidateSwitch, cannot get kernel");
        }
    }

    template <typename PixelT>
    double KernelCandidate<PixelT>::getBackground(CandidateSwitch cand) const {
        if (cand == KernelCandidate::ORIG) {
            if (_kernelSolutionOrig) 
                return _kernelSolutionOrig->getBackground();
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Orignal kernel does not exist");
        }
        else if (cand == KernelCandidate::PCA) {
            if (_kernelSolutionPca) 
                return _kernelSolutionPca->getBackground();
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Pca kernel does not exist");
        }
        else if (cand == KernelCandidate::RECENT) {
            if (_kernelSolutionPca) 
                return _kernelSolutionPca->getBackground();
            else if (_kernelSolutionOrig)
                return _kernelSolutionOrig->getBackground();
            else
                throw LSST_EXCEPT(pexExcept::Exception, "No kernels exist");
        }
        else {
            throw LSST_EXCEPT(pexExcept::Exception, "Invalid CandidateSwitch, cannot get background");
        }
    }

    template <typename PixelT>
    double KernelCandidate<PixelT>::getKsum(CandidateSwitch cand) const {
        if (cand == KernelCandidate::ORIG) {
            if (_kernelSolutionOrig) 
                return _kernelSolutionOrig->getKsum();
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Orignal kernel does not exist");
        }
        else if (cand == KernelCandidate::PCA) {
            if (_kernelSolutionPca) 
                return _kernelSolutionPca->getKsum();
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Pca kernel does not exist");
        }
        else if (cand == KernelCandidate::RECENT) {
            if (_kernelSolutionPca) 
                return _kernelSolutionPca->getKsum();
            else if (_kernelSolutionOrig)
                return _kernelSolutionOrig->getKsum();
            else
                throw LSST_EXCEPT(pexExcept::Exception, "No kernels exist");
        }
        else {
            throw LSST_EXCEPT(pexExcept::Exception, "Invalid CandidateSwitch, cannot get kSum");
        }
    }

    template <typename PixelT>
    KernelCandidate<PixelT>::ImageT::Ptr KernelCandidate<PixelT>::getKernelImage(
        CandidateSwitch cand) const {
        if (cand == KernelCandidate::ORIG) {
            if (_kernelSolutionOrig) 
                return _kernelSolutionOrig->makeKernelImage();
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Orignal kernel does not exist");
        }
        else if (cand == KernelCandidate::PCA) {
            if (_kernelSolutionPca) 
                return _kernelSolutionPca->makeKernelImage();
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Pca kernel does not exist");
        }
        else if (cand == KernelCandidate::RECENT) {
            if (_kernelSolutionPca) 
                return _kernelSolutionPca->makeKernelImage();
            else if (_kernelSolutionOrig)
                return _kernelSolutionOrig->makeKernelImage();
            else
                throw LSST_EXCEPT(pexExcept::Exception, "No kernels exist");
        }
        else {
            throw LSST_EXCEPT(pexExcept::Exception, "Invalid CandidateSwitch, cannot get kernel image");
        }
    }

    template <typename PixelT>
    KernelCandidate<PixelT>::ImageT::ConstPtr KernelCandidate<PixelT>::getImage() const {
        return getKernelImage(KernelCandidate::ORIG);
    }

    template <typename PixelT>
    boost::shared_ptr<StaticKernelSolution2<PixelT> > KernelCandidate<PixelT>::getKernelSolution(
        CandidateSwitch cand) const {
        if (cand == KernelCandidate::ORIG) {
            if (_kernelSolutionOrig) 
                return _kernelSolutionOrig;
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Orignal kernel does not exist");
        }
        else if (cand == KernelCandidate::PCA) {
            if (_kernelSolutionPca) 
                return _kernelSolutionPca;
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Pca kernel does not exist");
        }
        else if (cand == KernelCandidate::RECENT) {
            if (_kernelSolutionPca) 
                return _kernelSolutionPca;
            else if (_kernelSolutionOrig)
                return _kernelSolutionOrig;
            else
                throw LSST_EXCEPT(pexExcept::Exception, "No kernels exist");
        }
        else {
            throw LSST_EXCEPT(pexExcept::Exception, "Invalid CandidateSwitch, cannot get solution");
        }
    }

    template <typename PixelT>
    lsst::afw::image::MaskedImage<PixelT> KernelCandidate<PixelT>::getDifferenceImage(
        CandidateSwitch cand) {
        if (cand == KernelCandidate::ORIG) {
            if (_kernelSolutionOrig) 
                return getDifferenceImage(_kernelSolutionOrig->getKernel(),
                                          _kernelSolutionOrig->getBackground());
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Orignal kernel does not exist");
        }
        else if (cand == KernelCandidate::PCA) {
            if (_kernelSolutionPca) 
                return getDifferenceImage(_kernelSolutionPca->getKernel(),
                                          _kernelSolutionPca->getBackground());
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Pca kernel does not exist");
        }
        else if (cand == KernelCandidate::RECENT) {
            if (_kernelSolutionPca) 
                return getDifferenceImage(_kernelSolutionPca->getKernel(),
                                          _kernelSolutionPca->getBackground());
            else if (_kernelSolutionOrig)
                return getDifferenceImage(_kernelSolutionOrig->getKernel(),
                                          _kernelSolutionOrig->getBackground());
            else
                throw LSST_EXCEPT(pexExcept::Exception, "No kernels exist");
        }
        else {
            throw LSST_EXCEPT(pexExcept::Exception, "Invalid CandidateSwitch, cannot get diffim");
        }
    }

    template <typename PixelT>
    lsst::afw::image::MaskedImage<PixelT> KernelCandidate<PixelT>::getDifferenceImage(
        lsst::afw::math::Kernel::Ptr kernel,
        double background
        ) {
        /* Make diffim and set chi2 from result */
        afwImage::MaskedImage<PixelT> diffIm = convolveAndSubtract(*_miToConvolvePtr,
                                                                   *_miToNotConvolvePtr,
                                                                   *kernel,
                                                                   background);
        return diffIm;
    }

    
/***********************************************************************************************************/
//
// Explicit instantiations
//
    typedef float PixelT;

    template class KernelCandidate<PixelT>;

}}} // end of namespace lsst::ip::diffim
