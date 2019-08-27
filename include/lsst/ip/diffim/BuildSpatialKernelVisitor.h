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
#include "lsst/geom.h"
#include "lsst/ip/diffim.h"
#include "lsst/daf/base/PropertySet.h"

namespace lsst {
namespace ip {
namespace diffim {
namespace detail {

    template<typename PixelT>
    class BuildSpatialKernelVisitor : public lsst::afw::math::CandidateVisitor {
    public:
        typedef std::shared_ptr<BuildSpatialKernelVisitor<PixelT> > Ptr;

        BuildSpatialKernelVisitor(
            lsst::afw::math::KernelList const& basisList,  ///< Basis functions
            lsst::geom::Box2I const& regionBBox,  ///< Spatial region over which the function is fit
            lsst::daf::base::PropertySet ps           ///< PropertySet directing behavior
            );

        int getNCandidates() {return _nCandidates;}

        void processCandidate(lsst::afw::math::SpatialCellCandidate *candidate);

        void solveLinearEquation();

        inline std::shared_ptr<SpatialKernelSolution> getKernelSolution() {return _kernelSolution;}

        std::pair<std::shared_ptr<lsst::afw::math::LinearCombinationKernel>,
                  lsst::afw::math::Kernel::SpatialFunctionPtr> getSolutionPair();

    private:
        std::shared_ptr<SpatialKernelSolution> _kernelSolution;
        int _nCandidates;                  ///< Number of candidates visited
    };

    template<typename PixelT>
    std::shared_ptr<BuildSpatialKernelVisitor<PixelT> >
    makeBuildSpatialKernelVisitor(
        lsst::afw::math::KernelList const& basisList,
        lsst::geom::Box2I const& regionBBox,
        lsst::daf::base::PropertySet ps
        ) {

        return std::shared_ptr<BuildSpatialKernelVisitor<PixelT>>(
            new BuildSpatialKernelVisitor<PixelT>(basisList, regionBBox, ps)
            );
    }

}}}} // end of namespace lsst::ip::diffim::detail

#endif
