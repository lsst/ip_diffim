// -*- lsst-c++ -*-
/**
 * @file BasisLists.h
 *
 * @brief Subroutines associated with generating, normalising, and regularising Basis functions
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_BASISSLISTS_H
#define LSST_IP_DIFFIM_BASISSLISTS_H

#include <memory>
#include <vector>

#include "Eigen/Core"

#include "lsst/pex/policy/Policy.h"
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
    lsst::afw::math::KernelList makeDeltaFunctionBasisList(
        int width,
        int height
        );

    /**
     * @brief Build a regularization matrix for Delta function kernels
     *
     * @param policy           Policy file dictating which type of matrix to make
     *
     * @ingroup ip_diffim
     *
     * @note Calls either makeForwardDifferenceMatrix or
     * makeCentralDifferenceMatrix based on the policy file.
     */
    Eigen::MatrixXd makeRegularizationMatrix(
        lsst::pex::policy::Policy policy
        );

    /**
     * @brief Build a forward difference regularization matrix for Delta function kernels
     *
     * @param width            Width of basis set you want to regularize
     * @param height           Height of basis set you want to regularize
     * @param orders           Which derivatives to penalize (1,2,3)
     * @param borderPenalty    Amount of penalty (if any) to apply to border pixels; > 0
     * @param fitForBackground Fit for differential background?
     *
     * @ingroup ip_diffim
     */
    Eigen::MatrixXd makeForwardDifferenceMatrix(
        int width,
        int height,
        std::vector<int> const & orders,
        float borderPenalty,
        bool fitForBackground
        );

    /**
     * @brief Build a central difference Laplacian regularization matrix for Delta function kernels
     *
     * @param width            Width of basis set you want to regularize
     * @param height           Height of basis set you want to regularize
     * @param stencil          Which type of Laplacian approximation to use
     * @param borderPenalty    Amount of penalty (if any) to apply to border pixels; > 0
     * @param fitForBackground Fit for differential background?
     *
     * @ingroup ip_diffim
     */
    Eigen::MatrixXd makeCentralDifferenceMatrix(
        int width,
        int height,
        int stencil,
        float borderPenalty,
        bool fitForBackground
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
        lsst::afw::math::KernelList const &kernelListIn);

    lsst::afw::math::KernelList renormalizeKernelListScaled(
        lsst::afw::math::KernelList const &kernelListIn, double scaleSig
        );

    int renormalizeImageToZero(lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel> & kernelImage,
            double normSig);

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
    lsst::afw::math::KernelList makeAlardLuptonBasisList(
        int halfWidth,
        int nGauss,
        std::vector<double> const& sigGauss,
        std::vector<int>    const& degGauss,
		double scaleSig
        );

}}} // end of namespace lsst::ip::diffim

#endif
