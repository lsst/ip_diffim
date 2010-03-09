// -*- lsst-c++ -*-
/**
 * @file BasisSets.h
 *
 * @brief Subroutines associated with generating, normalising, and regularising Basis functions 
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_BASISSETS_H
#define LSST_IP_DIFFIM_BASISSETS_H

#include "boost/shared_ptr.hpp"

#include "Eigen/Core"

#include "lsst/afw/math/Kernel.h"

namespace lsst { 
namespace ip { 
namespace diffim {
   
    /**
     * @brief Build a set of Delta Function basis kernels
     * 
     * @note Total number of basis functions is width*height
     * 
     * @param width  Width of basis set (cols)
     * @param height Height of basis set (rows)
     *
     * @ingroup ip_diffim
     */    
    lsst::afw::math::KernelList generateDeltaFunctionBasisSet(
        unsigned int width,
        unsigned int height
        );

    /**
     * @brief Build a regularization matrix for Delta function kernels
     * 
     * @param width            Width of basis set you want to regularize
     * @param height           Height of basis set you want to regularize
     * @param order            Which derivative you expect to be smooth (derivative order+1 is penalized) 
     * @param boundary_style   0 = unwrapped, 1 = wrapped, 2 = order-tappered ('order' is highest used) 
     * @param difference_style 0 = forward, 1 = central
     * @param printB           debugging
     *
     * @ingroup ip_diffim
     */    
    boost::shared_ptr<Eigen::MatrixXd> generateFiniteDifferenceRegularization(
        unsigned int width,
        unsigned int height,
        unsigned int order,
	unsigned int boundary_style = 1, 
	unsigned int difference_style = 0,
	bool printB=false
        );

    boost::shared_ptr<Eigen::MatrixXd> foo(
        unsigned int width,
        unsigned int height,
        unsigned int order,
        int borderPenalty
        );

    /**
     * @brief Renormalize a list of basis kernels
     *
     * @note Renormalization means make Ksum_0 = 1.0, Ksum_i = 0.0, K_i.dot.K_i = 1.0
     * @note Output list of shared pointers to FixedKernels
     *
     * @param kernelListIn input list of basis kernels
     *
     * @note Images are checked for their current kernel sum.  If it is larger
     * than std::numeric_limits<double>::epsilon(), the kernel is first divided
     * by the kernel sum, giving it a kSum of 1.0, and then the first
     * (normalized) component is subtracted from it, giving it a kSum of 0.0.
     *
     * @ingroup ip_diffim
     */
    lsst::afw::math::KernelList renormalizeKernelList(
        lsst::afw::math::KernelList const &kernelListIn
        );

    /**
     * @brief Build a set of Alard/Lupton basis kernels
     *
     * @note Should consider implementing as SeparableKernels for additional speed,
     * but this will make the normalization a bit more complicated
     * 
     * @param halfWidth  size is 2*N + 1
     * @param nGauss     number of gaussians
     * @param sigGauss   Widths of the Gaussian Kernels
     * @param degGauss   Local spatial variation of bases
     *
     * @ingroup ip_diffim
     */    
    lsst::afw::math::KernelList generateAlardLuptonBasisSet(
        unsigned int halfWidth,                ///< size is 2*N + 1
        unsigned int nGauss,                   ///< number of gaussians
        std::vector<double> const& sigGauss,   ///< width of the gaussians
        std::vector<int>    const& degGauss    ///< local spatial variation of gaussians
        );

}}} // end of namespace lsst::ip::diffim

#endif
