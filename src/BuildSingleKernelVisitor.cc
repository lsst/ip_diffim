// -*- lsst-c++ -*-
/**
 * @file BuildSingleKernelVisitor.h
 *
 * @brief Implementation of BuildSingleKernelVisitor 
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include "boost/shared_ptr.hpp" 
#include "Eigen/Core"

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"

#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/ip/diffim/ImageSubtract.h"
#include "lsst/ip/diffim/KernelCandidate.h"
#include "lsst/ip/diffim/PsfMatchingFunctor.h"
#include "lsst/ip/diffim/BuildSingleKernelVisitor.h"

#define DEBUG_MATRIX 0

namespace afwMath        = lsst::afw::math;
namespace afwImage       = lsst::afw::image;
namespace pexLogging     = lsst::pex::logging; 
namespace pexPolicy      = lsst::pex::policy; 
namespace pexExcept      = lsst::pex::exceptions; 

namespace lsst { 
namespace ip { 
namespace diffim {
namespace detail {

    /**
     * @class BuildSingleKernelVisitor
     * @ingroup ip_diffim
     *
     * @brief Builds the convolution kernel for a given candidate
     *
     * @code
        Policy::Ptr policy(new Policy);
        policy->set("constantVarianceWeighting", false);
        policy->set("iterateSingleKernel", false);
        policy->set("singleKernelClipping", true);
        policy->set("candidateResidualMeanMax", 0.25);
        policy->set("candidateResidualStdMax", 1.25);
    
        detail::BuildSingleKernelVisitor<PixelT> singleKernelFitter(kFunctor, *policy);
        int nRejected = -1;
        while (nRejected != 0) {
            singleKernelFitter.reset();
            kernelCells.visitCandidates(&singleKernelFitter, nStarPerCell);
            nRejected = singleKernelFitter.getNRejected();
        }
     * @endcode
     *
     * @note Visits each current candidate in a afwMath::SpatialCellSet, and builds
     * its kernel using member object kFunctor.  We don't build the kernel for
     * *every* candidate since this is computationally expensive, only when its the
     * current candidate in the cell.  During the course of building the kernel, it
     * also assesses the quality of the difference image.  If it is determined to be
     * bad (based on the Policy paramters) the candidate is flagged as
     * afwMath::SpatialCellCandidate::BAD; otherwise its marked as
     * afwMath::SpatialCellCandidate::GOOD.  Keeps a running sample of all the new
     * candidates it visited that turned out to be bad.
     *
     * @note Because this visitor does not have access to the next candidate in the
     * cell, it must be called iteratively until no candidates are rejected.  This
     * ensures that the current candidate of every cell has an initialized Kernel.
     * This also requires that this class re-Visit all the cells after any other
     * Visitors with the ability to mark something as BAD.
     *
     * @note Because we are frequently re-Visiting entirely GOOD candidates during
     * these iterations, the option of _skipBuilt=true will enable the user to *not*
     * rebuilt the kernel on every visit.
     *
     * @note For the particular use case of creating a Pca basis from the raw
     * kernels, we want to re-Visit each candidate and re-fit the kernel using this
     * Pca basis.  However, we don't want to overwrite the original raw kernel,
     * since that is what is used to create the Pca basis in the first place.  Thus
     * the user has the option to setCandidateKernel(false), which will not override
     * the candidates original kernel, but will override its _M and _B matrices for
     * use in the spatial modeling.  This also requires the user to
     * setSkipBuilt(false) so that the candidate is reprocessed with this new basis.
     * 
     * @note When sending data to _kFunctor, this class uses an estimate of the
     * variance which is the straight difference of the 2 images.  If requested in
     * the Policy ("iterateSingleKernel"), the kernel will be rebuilt using the
     * variance of the difference image resulting from this first approximate step.
     * This is particularly useful when convolving a single-depth science image; the
     * variance (and thus resulting kernel) generally converges after 1 iteration.
     * If "constantVarianceWeighting" is requested in the Policy, no iterations will
     * be performed even if requested.
     * 
     */
    template<typename PixelT>
    BuildSingleKernelVisitor<PixelT>::BuildSingleKernelVisitor(
        PsfMatchingFunctor<PixelT> &kFunctor,    ///< Functor that builds the kernels
        lsst::pex::policy::Policy const& policy  ///< Policy file directing behavior
        ) :
        afwMath::CandidateVisitor(),
        _kFunctor(kFunctor),
        _policy(policy),
        _imstats(ImageStatistics<PixelT>()),
        _setCandidateKernel(true),
        _skipBuilt(true),
        _nRejected(0)
    {};

    
    template<typename PixelT>
    void BuildSingleKernelVisitor<PixelT>::processCandidate(
        lsst::afw::math::SpatialCellCandidate *candidate
        ) {
        
        KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
        if (kCandidate == NULL) {
            throw LSST_EXCEPT(pexExcept::LogicErrorException,
                              "Failed to cast SpatialCellCandidate to KernelCandidate");
        }
        
        if (_skipBuilt and kCandidate->hasKernel()) {
            return;
        }
        
        pexLogging::TTrace<3>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                              "Processing candidate %d", kCandidate->getId());
        
        /* Estimate of the variance */
        MaskedImageT var = MaskedImageT(*(kCandidate->getMiToNotConvolvePtr()), true);
        if (_policy.getBool("constantVarianceWeighting")) {
            /* Constant variance weighting */
            *var.getVariance() = 1;
        }
        else {
            /* Variance estimate is the straight difference */
            var -= *(kCandidate->getMiToConvolvePtr());
        }
        
        /* Build its kernel here */
        try {
            _kFunctor.apply(*(kCandidate->getMiToConvolvePtr()->getImage()),
                            *(kCandidate->getMiToNotConvolvePtr()->getImage()),
                            *(var.getVariance()),
                            _policy);
        } catch (pexExcept::Exception &e) {
            kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
            pexLogging::TTrace<4>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                  "Unable to process candidate %d; exception caught (%s)", 
                                  kCandidate->getId(),
                                  e.what());
            _nRejected += 1;
            return;
        }
        
        /* If we need to renormalize the kernel and its B matrix, do it here.
           This is particularly relevant when you are building a kernel matching
           an image to a Gaussian model, when each star is a different
           brightness.  We don't want to have to scale each Gaussian model to
           the flux of the star; we just scale all the kernels to have the same
           kernel sum, but we have to also scale the B matrix so that this does
           not go awry in the spatial modeling.
        */
        if (_policy.getBool("psfMatchToGaussian")) {
            _kFunctor.normalizeKernel();
        }
        
        /* 
           Sometimes you do not want to override the kernel; e.g. on a
           second fitting loop after the results of the first fitting loop
           are used to define a PCA basis
        */
        std::pair<boost::shared_ptr<afwMath::Kernel>, double> kb;
        try {
            kb = _kFunctor.getSolution();
        } catch (pexExcept::Exception &e) {
            kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
            pexLogging::TTrace<4>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                  "Unable to process candidate %d; exception caught (%s)", 
                                  kCandidate->getId(),
                                  e.what());
            _nRejected += 1;
            return;
        }
        
        if (_setCandidateKernel) {
            kCandidate->setKernel(kb.first);
            kCandidate->setBackground(kb.second);
        }
        
        /* 
         * However you *always* need to reset M and B since these are used 
         * in the spatial fitting
         */
        std::pair<boost::shared_ptr<Eigen::MatrixXd>, boost::shared_ptr<Eigen::VectorXd> > mb = 
            _kFunctor.getAndClearMB();
        kCandidate->setM(mb.first);
        kCandidate->setB(mb.second);
        
        /* 
         * Make diffim and set chi2 from result.  Note that you need to send
         * the newly-derived kernel and background in the case that
         * _setCandidateKernel = false.
         */
        MaskedImageT diffim = kCandidate->returnDifferenceImage(kb.first, kb.second);
        
        /* 
         * Remake the kernel using the first iteration difference image
         * variance as a better estimate of the true diffim variance.  If
         * you are setting "constantVarianceWeighting" it makes no sense to
         * do this
         */
        if (_policy.getBool("iterateSingleKernel") && (!(_policy.getBool("constantVarianceWeighting")))) {
            try {
                _kFunctor.apply(*(kCandidate->getMiToConvolvePtr()->getImage()),
                                *(kCandidate->getMiToNotConvolvePtr()->getImage()),
                                *(diffim.getVariance()),
                                _policy);
            } catch (pexExcept::Exception &e) {
                LSST_EXCEPT_ADD(e, "Unable to recalculate Kernel");
                throw e;
            }
            
            if (_policy.getBool("psfMatchToGaussian")) {
                _kFunctor.normalizeKernel();
            }
            
            try {
                kb = _kFunctor.getSolution();
            } catch (pexExcept::Exception &e) {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
                pexLogging::TTrace<4>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                      "Unable to process candidate %d; exception caught (%s)", 
                                      kCandidate->getId(), 
                                      e.what());
                _nRejected += 1;
                return;
            }
            
            if (_setCandidateKernel) {
                kCandidate->setKernel(kb.first);
                kCandidate->setBackground(kb.second);
            }
            
            mb = _kFunctor.getAndClearMB();
            kCandidate->setM(mb.first);
            kCandidate->setB(mb.second);
            diffim = kCandidate->returnDifferenceImage(kb.first, kb.second);                
        }
        
        /* Official resids */
        _imstats.apply(diffim);
        kCandidate->setChi2(_imstats.getVariance());
        
        /* When using a Pca basis, we don't reset the kernel or background,
           so we need to evaluate these locally for the Trace */
        afwImage::Image<double> kImage(kb.first->getDimensions());
        double kSum = kb.first->computeImage(kImage, false);
        double background = kb.second;
        
        pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                              "Chi2 = %.2f", kCandidate->getChi2());
        pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate",
                              "X = %.2f Y = %.2f",
                              kCandidate->getXCenter(), 
                              kCandidate->getYCenter());
        pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate",
                              "Kernel Sum = %.3f", kSum);
        pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate",
                              "Background = %.3f", background);
        pexLogging::TTrace<4>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate",
                              "Candidate %d resids = %.2f +/- %.2f sigma (%d pix)",
                              kCandidate->getId(),
                              _imstats.getMean(),
                              _imstats.getRms(),
                              _imstats.getNpix());
        
        bool meanIsNan = std::isnan(_imstats.getMean());
        bool rmsIsNan  = std::isnan(_imstats.getRms());
        if (meanIsNan || rmsIsNan) {
            kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
            pexLogging::TTrace<4>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                  "Rejecting candidate %d, encountered NaN",
                                  kCandidate->getId());
            _nRejected += 1;
            return;
        }
        
        if (_policy.getBool("singleKernelClipping")) {
            if (fabs(_imstats.getMean()) > _policy.getDouble("candidateResidualMeanMax")) {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
                pexLogging::TTrace<4>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                      "Rejecting candidate %d; bad mean residual : |%.2f| > %.2f",
                                      kCandidate->getId(),
                                      _imstats.getMean(),
                                      _policy.getDouble("candidateResidualMeanMax"));
                _nRejected += 1;
            }
            else if (_imstats.getRms() > _policy.getDouble("candidateResidualStdMax")) {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
                pexLogging::TTrace<4>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                      "Rejecting candidate %d; bad residual rms : %.2f > %.2f",
                                      kCandidate->getId(),
                                      _imstats.getRms(),
                                      _policy.getDouble("candidateResidualStdMax"));
                _nRejected += 1;
            }
            else {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::GOOD);
                pexLogging::TTrace<4>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                      "Source kernel OK");
            }
        }
        else {
            kCandidate->setStatus(afwMath::SpatialCellCandidate::GOOD);
            pexLogging::TTrace<6>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                  "Sigma clipping not enabled");
        }
        
        /* Core resids */
        int candidateCoreRadius = _policy.getInt("candidateCoreRadius");
        _imstats.apply(diffim, candidateCoreRadius);
        pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate",
                              "Candidate %d core resids = %.2f +/- %.2f sigma (%d pix)",
                              kCandidate->getId(),
                              _imstats.getMean(),
                              _imstats.getRms(),
                              _imstats.getNpix());
        
    }

    typedef float PixelT;
    template class BuildSingleKernelVisitor<PixelT>;

}}}} // end of namespace lsst::ip::diffim::detail
