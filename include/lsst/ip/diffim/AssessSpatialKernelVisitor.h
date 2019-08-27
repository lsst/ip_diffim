// -*- lsst-c++ -*-
/**
 * @file AssessSpatialKernelVisitor.h
 *
 * @brief Declaration of AssessSpatialKernelVisitor
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_ASSESSSPATIALKERNELVISITOR_H
#define LSST_IP_DIFFIM_ASSESSSPATIALKERNELVISITOR_H

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/ip/diffim.h"
#include "lsst/daf/base/PropertySet.h"

namespace lsst {
namespace ip {
namespace diffim {
namespace detail {

    template<typename PixelT>
    class AssessSpatialKernelVisitor : public lsst::afw::math::CandidateVisitor {
        typedef lsst::afw::image::MaskedImage<PixelT> MaskedImageT;
    public:
        typedef std::shared_ptr<AssessSpatialKernelVisitor<PixelT> > Ptr;

        AssessSpatialKernelVisitor(
            std::shared_ptr<lsst::afw::math::LinearCombinationKernel> spatialKernel,   ///< Spatially varying kernel
            lsst::afw::math::Kernel::SpatialFunctionPtr spatialBackground, ///< Spatially varying background
            lsst::daf::base::PropertySet const& ps                        ///< config
            );
        virtual ~AssessSpatialKernelVisitor() {};

        void reset() {_nGood = 0; _nRejected = 0; _nProcessed = 0;}

        int getNGood() {return _nGood;}
        int getNRejected() {return _nRejected;}
        int getNProcessed() {return _nProcessed;}
        void processCandidate(lsst::afw::math::SpatialCellCandidate *candidate);

    private:
        std::shared_ptr<lsst::afw::math::LinearCombinationKernel> _spatialKernel;   ///< Spatial kernel function
        lsst::afw::math::Kernel::SpatialFunctionPtr _spatialBackground; ///< Spatial background function
        lsst::daf::base::PropertySet _ps;     ///< Configuration controlling behavior
        ImageStatistics<PixelT> _imstats;     ///< To calculate statistics of difference image
        int _nGood;                           ///< Number of good candidates remaining
        int _nRejected;                       ///< Number of candidates rejected during processCandidate()
        int _nProcessed;                      ///< Number of candidates processed during processCandidate()

        bool _useCoreStats;                   ///< Extracted from config
        int _coreRadius;                      ///< Extracted from config
    };

    template<typename PixelT>
    std::shared_ptr<AssessSpatialKernelVisitor<PixelT> >
    makeAssessSpatialKernelVisitor(
        std::shared_ptr<lsst::afw::math::LinearCombinationKernel> spatialKernel,
        lsst::afw::math::Kernel::SpatialFunctionPtr spatialBackground,
        lsst::daf::base::PropertySet const& ps
         ) {

        return std::shared_ptr<AssessSpatialKernelVisitor<PixelT>>(
            new AssessSpatialKernelVisitor<PixelT>(spatialKernel, spatialBackground, ps)
            );
    }

}}}} // end of namespace lsst::ip::diffim::detail

#endif
