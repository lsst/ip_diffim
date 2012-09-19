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
#include "lsst/ip/diffim/ImageStatistics.h"
#include "lsst/ip/diffim/KernelSolution.h"

namespace afwMath        = lsst::afw::math;
namespace afwImage       = lsst::afw::image;
namespace pexLog         = lsst::pex::logging; 
namespace pexExcept      = lsst::pex::exceptions; 
namespace pexLogging     = lsst::pex::logging; 

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
        lsst::afw::math::SpatialCellImageCandidate<lsst::afw::math::Kernel::Pixel>(xCenter, yCenter),
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
        ImageStatistics<PixelT> imstats(_policy);
        int candidateCoreRadius = _policy.getInt("candidateCoreRadius");
        try {
            imstats.apply(*_miToNotConvolvePtr, candidateCoreRadius);
        } catch (pexExcept::Exception& e) {
            pexLogging::TTrace<3>("lsst.ip.diffim.KernelCandidate", 
                                  "Unable to calculate core imstats for ranking Candidate %d", this->getId()); 
            this->setStatus(afwMath::SpatialCellCandidate::BAD);
            return;
        }

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
        afwImage::Image<afwImage::VariancePixel> var = 
            afwImage::Image<afwImage::VariancePixel>(*(_miToNotConvolvePtr->getVariance()), true);
        /* Variance estimate comes from sum of image variances */
        var += (*(_miToConvolvePtr->getVariance()));

        if (_policy.getBool("constantVarianceWeighting")) {
            /* Constant variance weighting */
            afwMath::Statistics varStats = afwMath::makeStatistics(var, afwMath::MEDIAN);
            float varValue;
            if (varStats.getValue(afwMath::MEDIAN) <= 0.0)
                varValue = 1.0;
            else
                varValue = varStats.getValue(afwMath::MEDIAN);
            pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate",
                              "Candidate %d using constant variance of %.2f", varValue);
            var = varValue;

        }

        _varianceEstimate = VariancePtr( new afwImage::Image<afwImage::VariancePixel>(var) );

        try {
            _buildKernelSolution(basisList, hMat);
        } catch (pexExcept::Exception &e) {
            throw e;
        }
        

        if (_policy.getBool("iterateSingleKernel") && (!(_policy.getBool("constantVarianceWeighting")))) {
            afwImage::MaskedImage<PixelT> diffim = getDifferenceImage(KernelCandidate::RECENT);
            _varianceEstimate = diffim.getVariance();

            try {
                _buildKernelSolution(basisList, hMat);
            } catch (pexExcept::Exception &e) {
                throw e;
            }
        }

        _isInitialized = true;

    }

    template <typename PixelT>
    void KernelCandidate<PixelT>::_buildKernelSolution(lsst::afw::math::KernelList const& basisList,
                                                       boost::shared_ptr<Eigen::MatrixXd> hMat) 
    {
        bool checkConditionNumber = _policy.getBool("checkConditionNumber");
        double maxConditionNumber = _policy.getDouble("maxConditionNumber");
        std::string conditionNumberType = _policy.getString("conditionNumberType");
        KernelSolution::ConditionNumberType ctype;
        if (conditionNumberType == "SVD") {
            ctype = KernelSolution::SVD;
        }
        else if (conditionNumberType == "EIGENVALUE") {
            ctype = KernelSolution::EIGENVALUE;
        }
        else {
            throw LSST_EXCEPT(pexExcept::Exception, "conditionNumberType not recognized");            
        }
        
        /* Do we have a regularization matrix?  If so use it */
        if (hMat) {
            _useRegularization = true;
            pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate.build", 
                              "Using kernel regularization");

            if (_isInitialized) {
                _kernelSolutionPca = boost::shared_ptr<StaticKernelSolution<PixelT> >(
                    new RegularizedKernelSolution<PixelT>(basisList, _fitForBackground, hMat, _policy)
                    );
                _kernelSolutionPca->build(*(_miToConvolvePtr->getImage()),
                                          *(_miToNotConvolvePtr->getImage()),
                                          *_varianceEstimate);
                if (checkConditionNumber) {
                    if (_kernelSolutionPca->getConditionNumber(ctype) > maxConditionNumber) {
                        pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate",
                                          "Candidate %d solution has bad condition number", 
                                          this->getId());
                        this->setStatus(afwMath::SpatialCellCandidate::BAD);
                        return;
                    }
                }
                _kernelSolutionPca->solve();
            }
            else {
                _kernelSolutionOrig = boost::shared_ptr<StaticKernelSolution<PixelT> >(
                    new RegularizedKernelSolution<PixelT>(basisList, _fitForBackground, hMat, _policy)
                    );
                _kernelSolutionOrig->build(*(_miToConvolvePtr->getImage()),
                                           *(_miToNotConvolvePtr->getImage()),
                                           *_varianceEstimate);
                if (checkConditionNumber) {
                    if (_kernelSolutionOrig->getConditionNumber(ctype) > maxConditionNumber) {
                        pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate",
                                          "Candidate %d solution has bad condition number", 
                                          this->getId());
                        this->setStatus(afwMath::SpatialCellCandidate::BAD);
                        return;
                    }
                }
                _kernelSolutionOrig->solve();
            }
        }
        else {
            _useRegularization = false;
            pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate.build", 
                              "Not using kernel regularization");
            if (_isInitialized) {
                _kernelSolutionPca = boost::shared_ptr<StaticKernelSolution<PixelT> >(
                    new StaticKernelSolution<PixelT>(basisList, _fitForBackground)
                    );
                _kernelSolutionPca->build(*(_miToConvolvePtr->getImage()),
                                          *(_miToNotConvolvePtr->getImage()),
                                          *_varianceEstimate);
                if (checkConditionNumber) {
                    if (_kernelSolutionPca->getConditionNumber(ctype) > maxConditionNumber) {
                        pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate",
                                          "Candidate %d solution has bad condition number", 
                                          this->getId());
                        this->setStatus(afwMath::SpatialCellCandidate::BAD);
                        return;
                    }
                }
                _kernelSolutionPca->solve();
            }
            else {
                _kernelSolutionOrig = boost::shared_ptr<StaticKernelSolution<PixelT> >(
                    new StaticKernelSolution<PixelT>(basisList, _fitForBackground)
                    );
                _kernelSolutionOrig->build(*(_miToConvolvePtr->getImage()),
                                           *(_miToNotConvolvePtr->getImage()),
                                           *_varianceEstimate);
                if (checkConditionNumber) {
                    if (_kernelSolutionOrig->getConditionNumber(ctype) > maxConditionNumber) {
                        pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate",
                                          "Candidate %d solution has bad condition number", 
                                          this->getId());
                        this->setStatus(afwMath::SpatialCellCandidate::BAD);
                        return;
                    }
                }
                _kernelSolutionOrig->solve();
            }
        }
    }
            
    template <typename PixelT>
    lsst::afw::math::Kernel::Ptr KernelCandidate<PixelT>::getKernel(CandidateSwitch cand) const {
        if (cand == KernelCandidate::ORIG) {
            if (_kernelSolutionOrig) 
                return _kernelSolutionOrig->getKernel();
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Original kernel does not exist");
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
                throw LSST_EXCEPT(pexExcept::Exception, "Original kernel does not exist");
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
                throw LSST_EXCEPT(pexExcept::Exception, "Original kernel does not exist");
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
                throw LSST_EXCEPT(pexExcept::Exception, "Original kernel does not exist");
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
    boost::shared_ptr<StaticKernelSolution<PixelT> > KernelCandidate<PixelT>::getKernelSolution(
        CandidateSwitch cand) const {
        if (cand == KernelCandidate::ORIG) {
            if (_kernelSolutionOrig) 
                return _kernelSolutionOrig;
            else 
                throw LSST_EXCEPT(pexExcept::Exception, "Original kernel does not exist");
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
                throw LSST_EXCEPT(pexExcept::Exception, "Original kernel does not exist");
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
