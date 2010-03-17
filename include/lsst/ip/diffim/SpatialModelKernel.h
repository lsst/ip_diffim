// -*- lsst-c++ -*-
/**
 * @file SpatialModelKernel.h
 *
 * @brief Class used by SpatialModelCell for spatial Kernel fitting
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_SPATIALMODELKERNEL_H
#define LSST_IP_DIFFIM_SPATIALMODELKERNEL_H

#include "boost/shared_ptr.hpp"
#include "Eigen/Core"

#include "lsst/afw/math/SpatialCell.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/math/KernelFunctions.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/afw/detection/Footprint.h"

namespace lsst { 
namespace ip { 
namespace diffim {

    /**
     * @brief Fit for a spatial model of the kernel and background
     *
     * @param kernelCells  SpatialCellSet containing the candidate kernels
     * @param policy  Policy for configuration
     *
     * @ingroup ip_diffim
     */
    template<typename PixelT>
    std::pair<lsst::afw::math::LinearCombinationKernel::Ptr, lsst::afw::math::Kernel::SpatialFunctionPtr>
    fitSpatialKernelFromCandidates(
        lsst::afw::math::SpatialCellSet &kernelCells,
        lsst::pex::policy::Policy const& policy
        );


}}} // end of namespace lsst::ip::diffim

#endif

