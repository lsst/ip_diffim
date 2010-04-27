// -*- lsst-c++ -*-
/**
 * @file BuildSpatialKernelVisitor.h
 *
 * @brief Implementation of BuildSpatialKernelVisitor 
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include "boost/shared_ptr.hpp" 
#include "boost/timer.hpp" 

#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/LU"
#include "Eigen/QR"

#include "lsst/afw/math.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/ip/diffim/KernelCandidate.h"
#include "lsst/ip/diffim/KernelSolution.h"
#include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"

namespace afwMath        = lsst::afw::math;
namespace pexLogging     = lsst::pex::logging; 
namespace pexPolicy      = lsst::pex::policy; 
namespace pexExcept      = lsst::pex::exceptions; 

namespace lsst { 
namespace ip { 
namespace diffim {
namespace detail {
    /**
     * @class BuildSpatialKernelVisitor
     * @ingroup ip_diffim
     *
     * @brief Creates a spatial kernel and background from a list of candidates
     *
     * @code
        Policy::Ptr policy(new Policy);
        policy->set("spatialKernelOrder", spatialKernelOrder);
        policy->set("spatialBgOrder", spatialBgOrder);
        policy->set("kernelBasisSet", "delta-function");
        policy->set("usePcaForSpatialKernel", true);
    
        detail::BuildSpatialKernelVisitor<PixelT> spatialKernelFitter(*basisListToUse, 
                                                                      *policy);
        kernelCells.visitCandidates(&spatialKernelFitter, nStarPerCell);
        spatialKernelFitter.solveLinearEquation();
        std::pair<afwMath::LinearCombinationKernel::Ptr, 
            afwMath::Kernel::SpatialFunctionPtr> kb = spatialKernelFitter.getKernelSolution();
        spatialKernel     = kb.first;
        spatialBackground = kb.second; 
     * @endcode
     *
     * @note After visiting all candidates, solveLinearEquation() must be called to
     * trigger the matrix math.
     *
     * @note The user has the option to enfore conservation of the kernel sum across
     * the image through the policy.  In this case, all terms but the first are fit
     * for spatial variation.  This requires a little extra code to make sure the
     * matrices are the correct size, and that it is accessing the appropriate terms
     * in the matrices when creating the spatial models.
     * 
     */
    template<typename PixelT>
    BuildSpatialKernelVisitor<PixelT>::BuildSpatialKernelVisitor(
        lsst::afw::math::KernelList const& basisList, ///< Basis functions used in the fit
        lsst::pex::policy::Policy policy         ///< Policy file directing behavior
        ) :
        afwMath::CandidateVisitor(),
        _kernelSolution(new SpatialKernelSolution(basisList, policy)),
        _nCandidates(0) 
    {};
    
    
    template<typename PixelT>
    void BuildSpatialKernelVisitor<PixelT>::processCandidate(
        lsst::afw::math::SpatialCellCandidate *candidate
        ) {
        KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
        if (kCandidate == NULL) {
            throw LSST_EXCEPT(pexExcept::LogicErrorException,
                              "Failed to cast SpatialCellCandidate to KernelCandidate");
        }
        if (!(kCandidate->isInitialized())) {
            kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
            pexLogging::TTrace<3>("lsst.ip.diffim.BuildSpatialKernelVisitor.processCandidate", 
                                  "Cannot process candidate %d, continuing", kCandidate->getId());
            return;
        }
        
        pexLogging::TTrace<6>("lsst.ip.diffim.BuildSpatialKernelVisitor.processCandidate", 
                              "Processing candidate %d", kCandidate->getId());
        _nCandidates += 1;

        _kernelSolution->addConstraint(kCandidate->getXCenter(),
                                       kCandidate->getYCenter(),
                                       kCandidate->getKernelSolution(
                                           KernelCandidate<PixelT>::RECENT
                                           )->getM(),
                                       kCandidate->getKernelSolution(
                                           KernelCandidate<PixelT>::RECENT
                                           )->getB());
        
    }

    template<typename PixelT>
    void BuildSpatialKernelVisitor<PixelT>::solveLinearEquation() {
        _kernelSolution->solve();
    }

    template<typename PixelT>
    std::pair<afwMath::LinearCombinationKernel::Ptr, afwMath::Kernel::SpatialFunctionPtr>
    BuildSpatialKernelVisitor<PixelT>::getKernelSolution() {
        return _kernelSolution->getKernelSolution();
    }

    typedef float PixelT;
    template class BuildSpatialKernelVisitor<PixelT>;
  
}}}} // end of namespace lsst::ip::diffim::detail
    
