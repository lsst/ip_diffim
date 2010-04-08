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
#include "lsst/pex/policy/Policy.h"

namespace lsst { 
namespace ip { 
namespace diffim { 
namespace detail {

    template<typename PixelT>
    class AssessSpatialKernelVisitor : public lsst::afw::math::CandidateVisitor {
        typedef lsst::afw::image::MaskedImage<PixelT> MaskedImageT;
    public:

        AssessSpatialKernelVisitor(
            lsst::afw::math::LinearCombinationKernel::Ptr spatialKernel,   ///< Spatially varying kernel 
            lsst::afw::math::Kernel::SpatialFunctionPtr spatialBackground, ///< Spatially varying background
            lsst::pex::policy::Policy const& policy                        ///< Policy file 
            );
        virtual ~AssessSpatialKernelVisitor() {};

        void reset() {_nGood = 0; _nRejected = 0;}

        int getNGood() {return _nGood;}
        int getNRejected() {return _nRejected;}
        void processCandidate(lsst::afw::math::SpatialCellCandidate *candidate);

    private:
        lsst::afw::math::LinearCombinationKernel::Ptr _spatialKernel;   ///< Spatial kernel function
        lsst::afw::math::Kernel::SpatialFunctionPtr _spatialBackground; ///< Spatial background function
        lsst::pex::policy::Policy _policy;            ///< Policy controlling behavior
        ImageStatistics<PixelT> _imstats;     ///< To calculate statistics of difference image
        int _nGood;                           ///< Number of good candidates remaining
        int _nRejected;                       ///< Number of candidates rejected during processCandidate()
    };

}}}} // end of namespace lsst::ip::diffim::detail

#endif
