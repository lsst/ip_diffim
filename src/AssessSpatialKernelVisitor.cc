// -*- lsst-c++ -*-
/**
 * @file AssessSpatialKernelVisitor.cc
 *
 * @brief Implementation of AssessSpatialKernelVisitor
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/ip/diffim/ImageSubtract.h"
#include "lsst/ip/diffim/KernelCandidate.h"
#include "lsst/ip/diffim/AssessSpatialKernelVisitor.h"

namespace afwMath        = lsst::afw::math;
namespace afwImage       = lsst::afw::image;
namespace pexLogging     = lsst::pex::logging; 
namespace pexPolicy      = lsst::pex::policy; 
namespace pexExcept      = lsst::pex::exceptions; 

namespace lsst { 
namespace ip { 
namespace diffim {
namespace detail {
    /**
     * @class AssessSpatialKernelVisitor
     * @ingroup ip_diffim
     *
     * @brief Asseses the quality of a candidate given a spatial kernel and background model
     *
     * @code
        detail::AssessSpatialKernelVisitor<PixelT> spatialKernelAssessor(spatialKernel, 
                                                                         spatialBackground, 
                                                                         policy);
        spatialKernelAssessor.reset();
        kernelCells.visitCandidates(&spatialKernelAssessor, nStarPerCell);
        nRejected = spatialKernelAssessor.getNRejected();
     * @endcode
     *
     * @note Evaluates the spatial kernel and spatial background at the location of
     * each candidate, and computes the resulting difference image.  Sets candidate
     * as afwMath::SpatialCellCandidate::GOOD/BAD if requested by the Policy.
     * 
     */
    template<typename PixelT>
    AssessSpatialKernelVisitor<PixelT>::AssessSpatialKernelVisitor(
        lsst::afw::math::LinearCombinationKernel::Ptr spatialKernel,   ///< Spatially varying kernel model
        lsst::afw::math::Kernel::SpatialFunctionPtr spatialBackground, ///< Spatially varying backgound model
        lsst::pex::policy::Policy const& policy                        ///< Policy file directing behavior
        ) : 
        afwMath::CandidateVisitor(),
        _spatialKernel(spatialKernel),
        _spatialBackground(spatialBackground),
        _policy(policy),
        _imstats(ImageStatistics<PixelT>()),
        _nGood(0),
        _nRejected(0),
        _nProcessed(0)
    {};

    template<typename PixelT>
    void AssessSpatialKernelVisitor<PixelT>::processCandidate(
        lsst::afw::math::SpatialCellCandidate *candidate
        ) {
        
        KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
        if (kCandidate == NULL) {
            throw LSST_EXCEPT(pexExcept::LogicErrorException,
                              "Failed to cast SpatialCellCandidate to KernelCandidate");
        }
        if (!(kCandidate->isInitialized())) {
            kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
            pexLogging::TTrace<3>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                  "Cannot process candidate %d, continuing", kCandidate->getId());
            return;
        }
        
        pexLogging::TTrace<3>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                              "Processing candidate %d", kCandidate->getId());

        /* 
           Note - this is a hack until the Kernel API is upgraded by the
           Davis crew.  I need a "local" version of the spatially varying
           Kernel
        */
        afwImage::Image<double> kImage(_spatialKernel->getDimensions());
        double kSum = _spatialKernel->computeImage(kImage, false, 
                                                   afwImage::indexToPosition(
                                                       static_cast<int>(kCandidate->getXCenter())),
                                                   afwImage::indexToPosition(
                                                       static_cast<int>(kCandidate->getYCenter())));
        boost::shared_ptr<afwMath::Kernel>
            kernelPtr(new afwMath::FixedKernel(kImage));
        /* </hack> */
        
        double background = (*_spatialBackground)(afwImage::indexToPosition(
                                                      static_cast<int>(kCandidate->getXCenter())),
                                                  afwImage::indexToPosition(
                                                      static_cast<int>(kCandidate->getYCenter())));
        
        MaskedImageT diffim = kCandidate->getDifferenceImage(kernelPtr, background);
        
        /* Official resids */
        _imstats.apply(diffim);
        _nProcessed += 1;
        
        pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                              "Chi2 = %.3f", _imstats.getVariance());
        pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate",
                              "X = %.2f Y = %.2f",
                              kCandidate->getXCenter(), 
                              kCandidate->getYCenter());
        pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate",
                              "Kernel Sum = %.3f", kSum);
        pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate",
                              "Background = %.3f", background);
        pexLogging::TTrace<4>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate",
                              "Candidate %d resids = %.3f +/- %.3f sigma (%d pix)",
                              kCandidate->getId(),
                              _imstats.getMean(),
                              _imstats.getRms(),
                              _imstats.getNpix());
        
        bool meanIsNan = std::isnan(_imstats.getMean());
        bool rmsIsNan  = std::isnan(_imstats.getRms());
        if (meanIsNan || rmsIsNan) {
            kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
            pexLogging::TTrace<4>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                  "Rejecting candidate %d, encountered NaN",
                                  kCandidate->getId());
            _nRejected += 1;
            return;
        }
        
        if (_policy.getBool("spatialKernelClipping")) {            
            if (fabs(_imstats.getMean()) > _policy.getDouble("candidateResidualMeanMax")) {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
                pexLogging::TTrace<4>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                      "Rejecting candidate %d; bad mean residual : |%.3f| > %.3f",
                                      kCandidate->getId(),
                                      _imstats.getMean(),
                                      _policy.getDouble("candidateResidualMeanMax"));
                _nRejected += 1;
            }
            else if (_imstats.getRms() > _policy.getDouble("candidateResidualStdMax")) {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
                pexLogging::TTrace<4>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                      "Rejecting candidate %d; bad residual rms : %.3f > %.3f",
                                      kCandidate->getId(),
                                      _imstats.getRms(),
                                      _policy.getDouble("candidateResidualStdMax"));
                _nRejected += 1;
            }
            else {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::GOOD);
                pexLogging::TTrace<4>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                      "Spatial kernel OK");
                _nGood += 1;
            }
        }
        else {
            kCandidate->setStatus(afwMath::SpatialCellCandidate::GOOD);
            pexLogging::TTrace<6>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                  "Sigma clipping not enabled");
            _nGood += 1;
        }
        
        /* Core resids */
        int candidateCoreRadius = _policy.getInt("candidateCoreRadius");
        _imstats.apply(diffim, candidateCoreRadius);
        pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate",
                              "Candidate %d core resids = %.3f +/- %.3f sigma (%d pix)",
                              kCandidate->getId(),
                              _imstats.getMean(),
                              _imstats.getRms(),
                              _imstats.getNpix());
        
        
    }

    typedef float PixelT;
    template class AssessSpatialKernelVisitor<PixelT>;

}}}} // end of namespace lsst::ip::diffim::detail
