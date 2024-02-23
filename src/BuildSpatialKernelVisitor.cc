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

#include <memory>

#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/LU"
#include "Eigen/QR"

#include "lsst/afw/math.h"
#include "lsst/geom.h"
#include "lsst/log/Log.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/pex/exceptions/Runtime.h"

#include "lsst/ip/diffim/KernelCandidate.h"
#include "lsst/ip/diffim/KernelSolution.h"
#include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"

namespace afwMath        = lsst::afw::math;
namespace geom           = lsst::afw::geom;
namespace dafBase        = lsst::daf::base;
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
        std::shared_ptr<PropertySet> ps(new PropertySet);
        ps->set("spatialKernelOrder", spatialKernelOrder);
        ps->set("spatialBgOrder", spatialBgOrder);
        ps->set("kernelBasisSet", "delta-function");
        ps->set("usePcaForSpatialKernel", true);

        detail::BuildSpatialKernelVisitor<PixelT> spatialKernelFitter(*basisListToUse,
                                                                      *ps);
        kernelCells.visitCandidates(&spatialKernelFitter, nStarPerCell);
        spatialKernelFitter.solveLinearEquation();
        std::pair<std::shared_ptr<afwMath::LinearCombinationKernel>,
            afwMath::Kernel::SpatialFunctionPtr> kb = spatialKernelFitter.getKernelSolution();
        spatialKernel     = kb.first;
        spatialBackground = kb.second;
     * @endcode
     *
     * @note After visiting all candidates, solveLinearEquation() must be called to
     * trigger the matrix math.
     *
     * @note The user has the option to enfore conservation of the kernel sum across
     * the image through the property set.  In this case, all terms but the first are fit
     * for spatial variation.  This requires a little extra code to make sure the
     * matrices are the correct size, and that it is accessing the appropriate terms
     * in the matrices when creating the spatial models.
     *
     */
    template<typename PixelT>
    BuildSpatialKernelVisitor<PixelT>::BuildSpatialKernelVisitor(
        lsst::afw::math::KernelList const& basisList, ///< Basis functions used in the fit
        lsst::geom::Box2I const& regionBBox,     ///< Spatial region over which the function is fit
        lsst::daf::base::PropertySet const& ps        ///< ps file directing behavior
        ) :
        afwMath::CandidateVisitor(),
        _kernelSolution(),
        _nCandidates(0)
    {
        int spatialKernelOrder = ps.getAsInt("spatialKernelOrder");
        afwMath::Kernel::SpatialFunctionPtr spatialKernelFunction;

        int fitForBackground = ps.getAsBool("fitForBackground");
        int spatialBgOrder   = fitForBackground ? ps.getAsInt("spatialBgOrder") : 0;
        afwMath::Kernel::SpatialFunctionPtr background;

        std::string spatialModelType = ps.getAsString("spatialModelType");
        if (spatialModelType == "chebyshev1") {
            spatialKernelFunction = afwMath::Kernel::SpatialFunctionPtr(
                new afwMath::Chebyshev1Function2<double>(spatialKernelOrder, geom::Box2D(regionBBox))
                );
            background = afwMath::Kernel::SpatialFunctionPtr(
                new afwMath::Chebyshev1Function2<double>(spatialBgOrder, geom::Box2D(regionBBox))
                );

        }
        else if (spatialModelType == "polynomial") {
            spatialKernelFunction = afwMath::Kernel::SpatialFunctionPtr(
                new afwMath::PolynomialFunction2<double>(spatialKernelOrder)
                );
            background = afwMath::Kernel::SpatialFunctionPtr(
                new afwMath::PolynomialFunction2<double>(spatialBgOrder)
                );
        }
        else {
            throw LSST_EXCEPT(pexExcept::Exception,
                              str(boost::format("Invalid type (%s) for spatial models") %
                                  spatialModelType));
        }

        /* */

        _kernelSolution = std::shared_ptr<SpatialKernelSolution>(
            new SpatialKernelSolution(basisList, spatialKernelFunction, background, ps));
    };


    template<typename PixelT>
    void BuildSpatialKernelVisitor<PixelT>::processCandidate(
        lsst::afw::math::SpatialCellCandidate *candidate
        ) {
        KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
        if (kCandidate == NULL) {
            throw LSST_EXCEPT(pexExcept::LogicError,
                              "Failed to cast SpatialCellCandidate to KernelCandidate");
        }
        if (!(kCandidate->isInitialized())) {
            kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
            LOGL_DEBUG("TRACE2.ip.diffim.BuildSpatialKernelVisitor.processCandidate",
                       "Cannot process candidate %d, continuing", kCandidate->getId());
            return;
        }

        LOGL_DEBUG("TRACE5.ip.diffim.BuildSpatialKernelVisitor.processCandidate",
                   "Processing candidate %d", kCandidate->getId());
        _nCandidates += 1;

        /*
           Build the spatial kernel from the most recent fit, e.g. if its Pca
           you want to build a spatial model on the Pca basis, not original
           basis
        */
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
    std::pair<std::shared_ptr<afwMath::LinearCombinationKernel>, afwMath::Kernel::SpatialFunctionPtr>
    BuildSpatialKernelVisitor<PixelT>::getSolutionPair() {
        return _kernelSolution->getSolutionPair();
    }

    typedef float PixelT;
    template class BuildSpatialKernelVisitor<PixelT>;

}}}} // end of namespace lsst::ip::diffim::detail
