// -*- lsst-c++ -*-
/**
 * @file KernelSumVisitor.h
 *
 * @brief Implementation of KernelSumVisitor 
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */
#include <limits>

#include "lsst/afw/math.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/ip/diffim/KernelCandidate.h"
#include "lsst/ip/diffim/KernelSumVisitor.h"

namespace afwMath        = lsst::afw::math;
namespace pexLogging     = lsst::pex::logging; 
namespace pexPolicy      = lsst::pex::policy; 
namespace pexExcept      = lsst::pex::exceptions; 

namespace lsst { 
namespace ip { 
namespace diffim {
namespace detail {
    
    /**
     * @class KernelSumVisitor
     * @ingroup ip_diffim
     *
     * @brief A class to accumulate kernel sums across SpatialCells 
     *
     * @code
        Policy::Ptr policy(new Policy);
        policy->set("kernelSumClipping", false);
        policy->set("maxKsumSigma", 3.0);
     
        detail::KernelSumVisitor<PixelT> kernelSumVisitor(*policy);
        kernelSumVisitor.reset();
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::AGGREGATE);
        kernelCells.visitCandidates(&kernelSumVisitor, nStarPerCell);
        kernelSumVisitor.processKsumDistribution();
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::REJECT);
        kernelCells.visitCandidates(&kernelSumVisitor, nStarPerCell);
        int nRejected = kernelSumVisitor.getNRejected();
     * @endcode
     *
     *
     * @note The class has 2 processing modes; the first AGGREGATES kernel sums
     * across all candidates.  You must the process the distribution to set member
     * variables representing the mean and standard deviation of the kernel sums.
     * The second mode then REJECTs candidates with kernel sums outside the
     * acceptable range (set by the policy).  It does this by setting candidate
     * status to afwMath::SpatialCellCandidate::BAD.  In this mode it also
     * accumulates the number of candidates it sets as bad.
     *
     * @note The statistics call calculates sigma-clipped values (afwMath::MEANCLIP,
     * afwMath::STDEVCLIP)
     *
     */
    template<typename PixelT>
    KernelSumVisitor<PixelT>::KernelSumVisitor(
        pexPolicy::Policy const& policy ///< Policy file directing behavior
        ) :
        afwMath::CandidateVisitor(),
        _mode(AGGREGATE),
        _kSums(std::vector<double>()),
        _kSumMean(0.),
        _kSumStd(0.),
        _dkSumMax(0.),
        _kSumNpts(0),
        _nRejected(0),
        _policy(policy) 
    {};
    
    template<typename PixelT>
    void KernelSumVisitor<PixelT>::resetKernelSum() {
        _kSums.clear();
        _kSumMean =  0.;
        _kSumStd  =  0.;
        _dkSumMax =  0.;
        _kSumNpts =  0;
        _nRejected = 0;
    }

    template<typename PixelT>
    void KernelSumVisitor<PixelT>::processCandidate(lsst::afw::math::SpatialCellCandidate 
                                                    *candidate) {

        KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
        if (kCandidate == NULL) {
            throw LSST_EXCEPT(pexExcept::LogicErrorException,
                              "Failed to cast SpatialCellCandidate to KernelCandidate");
        }
        pexLogging::TTrace<6>("lsst.ip.diffim.KernelSumVisitor.processCandidate", 
                              "Processing candidate %d, mode %d", kCandidate->getId(), _mode);
        
        /* Grab all kernel sums and look for outliers */
        if (_mode == AGGREGATE) {
            _kSums.push_back(kCandidate->getKernelSolution(KernelCandidate<PixelT>::ORIG)->getKsum());
        }
        else if (_mode == REJECT) {
            if (_policy.getBool("kernelSumClipping")) {
                double kSum = 
                    kCandidate->getKernelSolution(KernelCandidate<PixelT>::ORIG)->getKsum();

                if (fabs(kSum - _kSumMean) > _dkSumMax) {
                    kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
                    pexLogging::TTrace<4>("lsst.ip.diffim.KernelSumVisitor.processCandidate", 
                                          "Rejecting candidate %d; bad source kernel sum : (%.2f)",
                                          kCandidate->getId(),
                                          kSum);
                    _nRejected += 1;
                }
            }
            else {
                pexLogging::TTrace<6>("lsst.ip.diffim.KernelSumVisitor.processCandidate", 
                                      "Sigma clipping not enabled");
            }
        }
    }
    
    template<typename PixelT>
    void KernelSumVisitor<PixelT>::processKsumDistribution() {
        if (_kSums.size() == 0) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              "Unable to determine kernel sum; 0 candidates");
        }
        else if (_kSums.size() == 1) {
            pexLogging::TTrace<2>("lsst.ip.diffim.KernelSumVisitor.processKsumDistribution", 
                                  "WARNING: only 1 kernel candidate");
            
            _kSumMean = _kSums[0];
            _kSumStd  = 0.0;
            _kSumNpts = 1;
        }
        else {
            try {
                afwMath::Statistics stats = afwMath::makeStatistics(_kSums, 
                                                                    afwMath::NPOINT | 
                                                                    afwMath::MEANCLIP | 
                                                                    afwMath::STDEVCLIP); 
                _kSumMean = stats.getValue(afwMath::MEANCLIP);
                _kSumStd  = stats.getValue(afwMath::STDEVCLIP);
                _kSumNpts = static_cast<int>(stats.getValue(afwMath::NPOINT));
            } catch (pexExcept::Exception &e) {
                LSST_EXCEPT_ADD(e, "Unable to calculate kernel sum statistics");
                throw e;
            }
            if (std::isnan(_kSumMean)) {
                throw LSST_EXCEPT(pexExcept::Exception, 
                                  str(boost::format("Mean kernel sum returns NaN (%d points)") 
                                      % _kSumNpts));
            }
            if (std::isnan(_kSumStd)) {
                throw LSST_EXCEPT(pexExcept::Exception, 
                                  str(boost::format("Kernel sum stdev returns NaN (%d points)") 
                                      % _kSumNpts));
            }
        }
        _dkSumMax = _policy.getDouble("maxKsumSigma") * _kSumStd;
        pexLogging::TTrace<2>("lsst.ip.diffim.KernelSumVisitor.processCandidate", 
                              "Kernel Sum Distribution : %.3f +/- %.3f (%d points)", 
                              _kSumMean, _kSumStd, _kSumNpts);
    }
    
    typedef float PixelT;

    template class KernelSumVisitor<PixelT>;

    template boost::shared_ptr<KernelSumVisitor<PixelT> > 
    makeKernelSumVisitor<PixelT>(lsst::pex::policy::Policy const&);

}}}} // end of namespace lsst::ip::diffim::detail
