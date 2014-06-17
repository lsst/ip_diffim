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

#define DEBUG_IMAGES 0

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
        _imstats(ImageStatistics<PixelT>(_policy)),
        _nGood(0),
        _nRejected(0),
        _nProcessed(0),
        _useCoreStats(_policy.getBool("useCoreStats")),
        _coreRadius(_policy.getInt("candidateCoreRadius"))
    {};

    template<typename PixelT>
    void AssessSpatialKernelVisitor<PixelT>::processCandidate(
        lsst::afw::math::SpatialCellCandidate *candidate
        ) {
        
        KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
        if (kCandidate == NULL) {
            throw LSST_EXCEPT(pexExcept::LogicError,
                              "Failed to cast SpatialCellCandidate to KernelCandidate");
        }
        if (!(kCandidate->isInitialized())) {
            kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
            pexLogging::TTrace<3>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                  "Cannot process candidate %d, continuing", kCandidate->getId());
            return;
        }
        
        pexLogging::TTrace<2>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                              "Processing candidate %d", kCandidate->getId());

        /* 
           Note - this is a hack until the Kernel API is upgraded by the
           Davis crew.  I need a "local" version of the spatially varying
           Kernel
        */
        afwImage::Image<double> kImage(_spatialKernel->getDimensions());
        double kSum = _spatialKernel->computeImage(kImage, false, 
                                                   kCandidate->getXCenter(), kCandidate->getYCenter());
        boost::shared_ptr<afwMath::Kernel>
            kernelPtr(new afwMath::FixedKernel(kImage));
        /* </hack> */
        
        double background = (*_spatialBackground)(kCandidate->getXCenter(), kCandidate->getYCenter());
        
        MaskedImageT diffim = kCandidate->getDifferenceImage(kernelPtr, background);

        if (DEBUG_IMAGES) {
            kImage.writeFits(str(boost::format("askv_k%d.fits") % kCandidate->getId()));
            diffim.writeFits(str(boost::format("askv_d%d.fits") % kCandidate->getId()));
        }

        /* Official resids */
        try {
            if (_useCoreStats) 
                _imstats.apply(diffim, _coreRadius);
            else
                _imstats.apply(diffim);
        } catch (pexExcept::Exception& e) {
            pexLogging::TTrace<3>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                  "Unable to calculate imstats for Candidate %d", kCandidate->getId()); 
            kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
            return;
        }

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
        pexLogging::TTrace<3>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate",
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

        /* Core resids for debugging */
        if (!(_useCoreStats)) {
            try {
                _imstats.apply(diffim, _coreRadius);
            } catch (pexExcept::Exception& e) {
                pexLogging::TTrace<3>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                      "Unable to calculate core imstats for Candidate %d", 
                                      kCandidate->getId()); 
                kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
                return;
            }
            pexLogging::TTrace<4>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate",
                                  "Candidate %d core resids = %.3f +/- %.3f sigma (%d pix)",
                                  kCandidate->getId(),
                                  _imstats.getMean(),
                                  _imstats.getRms(),
                                  _imstats.getNpix());
        }
    }

    typedef float PixelT;
    template class AssessSpatialKernelVisitor<PixelT>;

}}}} // end of namespace lsst::ip::diffim::detail
