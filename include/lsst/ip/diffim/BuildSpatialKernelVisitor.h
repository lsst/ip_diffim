// -*- lsst-c++ -*-
/**
 * @file BuildSpatialKernelVisitor.h
 *
 * @brief Declaration of BuildSpatialKernelVisitor 
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_BUILDSPATIALKERNELVISITOR_H
#define LSST_IP_DIFFIM_BUILDSPATIALKERNELVISITOR_H

#include "Eigen/Core"
#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/pex/policy/Policy.h"

namespace lsst { 
namespace ip { 
namespace diffim { 
namespace detail {

    template<typename PixelT>
    class BuildSpatialKernelVisitor : public lsst::afw::math::CandidateVisitor {
    public:
        BuildSpatialKernelVisitor(
            lsst::afw::math::KernelList basisList,   ///< Basis functions used in the fit
            lsst::pex::policy::Policy policy         ///< Policy file directing behavior
            );

        int getNCandidates() { return _nCandidates; }

        void processCandidate(lsst::afw::math::SpatialCellCandidate *candidate);

        void solveLinearEquation();

        std::pair<lsst::afw::math::LinearCombinationKernel::Ptr, 
                  lsst::afw::math::Kernel::SpatialFunctionPtr> getSpatialModel();

    private:
        lsst::afw::math::KernelList _basisList; ///< List of kernel basis functions
        Eigen::MatrixXd _mMat;          ///< Least squares matrix
        Eigen::VectorXd _bVec;          ///< Least squares vector
        Eigen::VectorXd _sVec;          ///< Least squares solution
        int const _spatialKernelOrder;  ///< Spatial order of kernel variation
        int const _spatialBgOrder;      ///< Spatial order of background variation
        lsst::afw::math::Kernel::SpatialFunctionPtr _spatialKernelFunction; ///< Spatial kernel function
        lsst::afw::math::Kernel::SpatialFunctionPtr _spatialBgFunction;     ///< Spatial bg function
        unsigned int _nbases;      ///< Number of bases being fit for
        unsigned int _nkt;         ///< Number of kernel terms in spatial model
        unsigned int _nbt;         ///< Number of background terms in spatial model
        unsigned int _nt;          ///< Total number of terms in the solution; also dimensions of matrices
        lsst::pex::policy::Policy _policy; ///< Policy controlling behavior
        bool _constantFirstTerm;   ///< Is the first term spatially invariant?
        int _nCandidates;          ///< Number of candidates visited
    };

}}}} // end of namespace lsst::ip::diffim::detail

#endif
