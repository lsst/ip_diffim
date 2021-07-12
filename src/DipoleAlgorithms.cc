/*
 * LSST Data Management System
 * Copyright 2008-2015 AURA/LSST
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

/**
 * @file
 */
#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <limits>       // std::numeric_limits
#include <cmath>        // std::sqrt

#if !defined(DOXYGEN)
#   include "Minuit2/FCNBase.h"
#   include "Minuit2/FunctionMinimum.h"
#   include "Minuit2/MnMigrad.h"
#   include "Minuit2/MnMinos.h"
#   include "Minuit2/MnPrint.h"
#endif

#include <memory>
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/table.h"
#include "lsst/afw/math.h"
#include "lsst/geom.h"
#include "lsst/ip/diffim/DipoleAlgorithms.h"
#include "ndarray/eigen.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace geom = lsst::geom;

namespace lsst { namespace ip { namespace diffim {

namespace {
meas::base::FlagDefinitionList dipoleFluxFlagDefinitions;
}

meas::base::FlagDefinition const DipoleFluxAlgorithm::FAILURE = dipoleFluxFlagDefinitions.addFailureFlag("general failure flag, set if anything went wrong");
meas::base::FlagDefinition const DipoleFluxAlgorithm::POS_FLAG = dipoleFluxFlagDefinitions.add("pos_flag", "failure flag for positive, set if anything went wrong");
meas::base::FlagDefinition const DipoleFluxAlgorithm::NEG_FLAG = dipoleFluxFlagDefinitions.add("neg_flag", "failure flag for negative, set if anything went wrong");

meas::base::FlagDefinitionList const & DipoleFluxAlgorithm::getFlagDefinitions() {
    return dipoleFluxFlagDefinitions;
}

namespace {
meas::base::FlagDefinitionList dipoleCentroidFlagDefinitions;
}

meas::base::FlagDefinition const DipoleCentroidAlgorithm::FAILURE = dipoleCentroidFlagDefinitions.addFailureFlag("general failure flag, set if anything went wrong");
meas::base::FlagDefinition const DipoleCentroidAlgorithm::POS_FLAG = dipoleCentroidFlagDefinitions.add("pos_flag", "failure flag for positive, set if anything went wrong");
meas::base::FlagDefinition const DipoleCentroidAlgorithm::NEG_FLAG = dipoleCentroidFlagDefinitions.add("neg_flag", "failure flag for negative, set if anything went wrong");

meas::base::FlagDefinitionList const & DipoleCentroidAlgorithm::getFlagDefinitions() {
    return dipoleCentroidFlagDefinitions;
}

    int const NEGCENTXPAR(0); // Parameter for the x-component of the negative lobe centroid
    int const NEGCENTYPAR(1); // Parameter for the y-component of the negative lobe centroid
    int const NEGFLUXPAR(2);  // Parameter for the flux of the negative lobe
    int const POSCENTXPAR(3); // Parameter for the x-component of the positive lobe centroid
    int const POSCENTYPAR(4); // Parameter for the y-component of the positive lobe centroid
    int const POSFLUXPAR(5);  // Parameter for the flux of the positive lobe


namespace {

void naiveCentroid(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const& exposure,
    geom::Point2I const & center,
    meas::base::CentroidResultKey const & keys
    )
{
    typedef afw::image::Image<float> ImageT;
    ImageT const& image = *exposure.getMaskedImage().getImage();
    // set to the input centroid, just in case all else fails
    source.set(keys.getX(), center.getX());
    source.set(keys.getY(), center.getY());

    int x = center.getX() - image.getX0();
    int y = center.getY() - image.getY0();

    if (x < 1 || x >= image.getWidth() - 1 || y < 1 || y >= image.getHeight() - 1) {
         throw LSST_EXCEPT(pex::exceptions::LengthError,
                           (boost::format("Object at (%d, %d) is too close to the edge")
                            % x % y).str());
    }

    ImageT::xy_locator im = image.xy_at(x, y);

    double const sum =
        (im(-1,  1) + im( 0,  1) + im( 1,  1) +
         im(-1,  0) + im( 0,  0) + im( 1,  0) +
         im(-1, -1) + im( 0, -1) + im( 1, -1));


    if (sum == 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeError,
                          (boost::format("Object at (%d, %d) has no counts") %
                           x % y).str());
    }

    double const sum_x =
        -im(-1,  1) + im( 1,  1) +
        -im(-1,  0) + im( 1,  0) +
        -im(-1, -1) + im( 1, -1);
    double const sum_y =
        (im(-1,  1) + im( 0,  1) + im( 1,  1)) -
        (im(-1, -1) + im( 0, -1) + im( 1, -1));

    float xx = afw::image::indexToPosition(x + image.getX0()) + sum_x / sum;
    float yy = afw::image::indexToPosition(y + image.getY0()) + sum_y / sum;
    source.set(keys.getX(), xx);
    source.set(keys.getY(), yy);
}

} // anonymous namespace


NaiveDipoleCentroid::NaiveDipoleCentroid(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : DipoleCentroidAlgorithm(ctrl, name, schema, "unweighted first moment centroid"),
    _ctrl(ctrl)
{ }

/**
 * Given an image and a pixel position, return a Centroid using a naive 3x3 weighted moment
 */
void NaiveDipoleCentroid::measure(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {
    afw::detection::PeakCatalog const& peaks = source.getFootprint()->getPeaks();

    int posInd = 0;
    double posValue = peaks[posInd].getPeakValue(), negValue = 0;
    if (posValue < 0.) { /* All peaks are negative so use the *most* negative value */
       posInd = peaks.size() - 1;
       posValue = peaks[posInd].getPeakValue();
    }
    naiveCentroid(source, exposure, peaks[posInd].getI(),
                 (posValue >= 0 ? getPositiveKeys() : getNegativeKeys()));

    if (posValue > 0. && posInd == 0 && peaks.size() > 1) { /* See if there's also a negative peak */
       int negInd = peaks.size() - 1;
       negValue = peaks[negInd].getPeakValue();
       if (posValue > 0. && negValue < 0.) {
          naiveCentroid(source, exposure, peaks[negInd].getI(),
                        (negValue >= 0 ? getPositiveKeys() : getNegativeKeys()));
       }
    }

    mergeCentroids(source, posValue, negValue);

}

void NaiveDipoleCentroid::mergeCentroids(afw::table::SourceRecord & source,
                                         double posValue, double negValue) const {

   double pos_x, pos_y, pos_f;
   double neg_x, neg_y, neg_f;

    pos_x = source.get(getPositiveKeys().getX());
    pos_y = source.get(getPositiveKeys().getY());
    pos_f = posValue;

    neg_x = source.get(getNegativeKeys().getX());
    neg_y = source.get(getNegativeKeys().getY());
    neg_f = -negValue;

    if(std::isfinite(pos_x) && std::isfinite(pos_y) &&
       std::isfinite(neg_x) && std::isfinite(neg_y)) {
       source.set(getCenterKeys().getX(), (pos_x * pos_f + neg_x * neg_f) / (pos_f + neg_f));
       source.set(getCenterKeys().getY(), (pos_y * pos_f + neg_y * neg_f) / (pos_f + neg_f));
    } else if (std::isfinite(pos_x) && std::isfinite(pos_y)) {
        source.set(getCenterKeys().getX(), pos_x);
        source.set(getCenterKeys().getY(), pos_y);
    } else {
        source.set(getCenterKeys().getX(), neg_x);
        source.set(getCenterKeys().getY(), neg_y);
    }
}

void NaiveDipoleCentroid::fail(afw::table::SourceRecord & measRecord,
                               meas::base::MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}


namespace {

class NaiveDipoleFootprinter {
public:
    explicit NaiveDipoleFootprinter() : _sumPositive(0.0), _sumNegative(0.0), _numPositive(0),
                                        _numNegative(0) {}

    /// Reset everything for a new Footprint
    void reset() {
        _sumPositive = _sumNegative = 0.0;
        _numPositive = _numNegative = 0;
    }
    void reset(afwDet::Footprint const&) {}

    /// method called for each pixel by applyFunctor
    void operator() (geom::Point2I const & pos,
                     afw::image::MaskedImage<float>::Image::Pixel const & ival,
                     afw::image::MaskedImage<float>::Image::Pixel const & vval) {
        if (ival >= 0.0) {
            _sumPositive += ival;
            _varPositive += vval;
            ++_numPositive;
        } else {
            _sumNegative += ival;
            _varPositive += vval;
            ++_numNegative;
        }
    }

    double getSumPositive() const { return _sumPositive; }
    double getSumNegative() const { return _sumNegative; }
    double getVarPositive() const { return _sumPositive; }
    double getVarNegative() const { return _sumNegative; }
    int getNumPositive() const { return _numPositive; }
    int getNumNegative() const { return _numNegative; }

private:
    double _sumPositive;
    double _sumNegative;
    double _varPositive;
    double _varNegative;
    int _numPositive;
    int _numNegative;
};

} // anonymous namespace


/**
 * Given an image and a pixel position, return a Centroid using a naive 3x3 weighted moment
 */
void NaiveDipoleFlux::measure(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {
    typedef afw::image::Exposure<float>::MaskedImageT MaskedImageT;

    NaiveDipoleFootprinter functor;
    source.getFootprint()->getSpans()->applyFunctor(functor, *(exposure.getMaskedImage().getImage()),
                                                     *(exposure.getMaskedImage().getVariance()));

    source.set(getPositiveKeys().getInstFlux(), functor.getSumPositive());
    source.set(getPositiveKeys().getInstFluxErr(), ::sqrt(functor.getVarPositive()));
    source.set(_numPositiveKey, functor.getNumPositive());

    source.set(getNegativeKeys().getInstFlux(), functor.getSumNegative());
    source.set(getNegativeKeys().getInstFluxErr(), ::sqrt(functor.getVarNegative()));
    source.set(_numNegativeKey, functor.getNumNegative());
    functor.reset();
}

void NaiveDipoleFlux::fail(afw::table::SourceRecord & measRecord, meas::base::MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}


/**
 * Class to minimize PsfDipoleFlux; this is the object that Minuit minimizes
 */
class MinimizeDipoleChi2 : public ROOT::Minuit2::FCNBase {
public:
    explicit MinimizeDipoleChi2(PsfDipoleFlux const& psfDipoleFlux,
                                afw::table::SourceRecord & source,
                                afw::image::Exposure<float> const& exposure
                                ) : _errorDef(1.0),
                                    _nPar(6),
                                    _maxPix(1e4),
                                    _bigChi2(1e10),
                                    _psfDipoleFlux(psfDipoleFlux),
                                    _source(source),
                                    _exposure(exposure)
    {}
    double Up() const { return _errorDef; }
    void setErrorDef(double def) { _errorDef = def; }
    int getNpar() const { return _nPar; }
    int getMaxPix() const { return _maxPix; }
    void setMaxPix(int maxPix) { _maxPix = maxPix; }

    // Evaluate our cost function (in this case chi^2)
    virtual double operator()(std::vector<double> const & params) const {
        double negCenterX = params[NEGCENTXPAR];
        double negCenterY = params[NEGCENTYPAR];
        double negFlux    = params[NEGFLUXPAR];
        double posCenterX = params[POSCENTXPAR];
        double posCenterY = params[POSCENTYPAR];
        double posFlux    = params[POSFLUXPAR];

        /* Restrict negative dipole to be negative; positive to be positive */
        if ((negFlux > 0.0) || (posFlux < 0.0)) {
            return _bigChi2;
        }

        std::pair<double,int> fit = _psfDipoleFlux.chi2(_source, _exposure, negCenterX, negCenterY, negFlux,
                                                        posCenterX, posCenterY, posFlux);
        double chi2 = fit.first;
        int nPix = fit.second;
        if (nPix > _maxPix) {
            return _bigChi2;
        }

        return chi2;
    }

private:
    double _errorDef;       // how much cost function has changed at the +- 1 error points
    int _nPar;              // number of parameters in the fit; hard coded for MinimizeDipoleChi2
    int _maxPix;            // maximum number of pixels that shoud be in the footprint;
                            // prevents too much centroid wander
    double _bigChi2;        // large value to tell fitter when it has gone into bad region of parameter space

    PsfDipoleFlux const& _psfDipoleFlux;
    afw::table::SourceRecord & _source;
    afw::image::Exposure<float> const& _exposure;
};

std::pair<double,int> PsfDipoleFlux::chi2(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const& exposure,
    double negCenterX, double negCenterY, double negFlux,
    double posCenterX, double posCenterY, double posFlux
) const {

    geom::Point2D negCenter(negCenterX, negCenterY);
    geom::Point2D posCenter(posCenterX, posCenterY);

    std::shared_ptr<afw::detection::Footprint const> footprint = source.getFootprint();

    /*
     * Fit for the superposition of Psfs at the two centroids.
     */
    std::shared_ptr<afwDet::Psf const> psf = exposure.getPsf();
    std::shared_ptr<afwImage::Image<afwMath::Kernel::Pixel>> negPsf = psf->computeImage(negCenter);
    std::shared_ptr<afwImage::Image<afwMath::Kernel::Pixel>> posPsf = psf->computeImage(posCenter);

    afwImage::Image<double> negModel(footprint->getBBox());
    afwImage::Image<double> posModel(footprint->getBBox());
    afwImage::Image<float> data(*(exposure.getMaskedImage().getImage()),footprint->getBBox());
    afwImage::Image<afwImage::VariancePixel> var(*(exposure.getMaskedImage().getVariance()),
                                                 footprint->getBBox());

    geom::Box2I negPsfBBox = negPsf->getBBox();
    geom::Box2I posPsfBBox = posPsf->getBBox();
    geom::Box2I negModelBBox = negModel.getBBox();
    geom::Box2I posModelBBox = posModel.getBBox();

    // Portion of the negative Psf that overlaps the model
    int negXmin = std::max(negPsfBBox.getMinX(), negModelBBox.getMinX());
    int negYmin = std::max(negPsfBBox.getMinY(), negModelBBox.getMinY());
    int negXmax = std::min(negPsfBBox.getMaxX(), negModelBBox.getMaxX());
    int negYmax = std::min(negPsfBBox.getMaxY(), negModelBBox.getMaxY());
    geom::Box2I negBBox = geom::Box2I(geom::Point2I(negXmin, negYmin),
                                      geom::Point2I(negXmax, negYmax));
    afwImage::Image<afwMath::Kernel::Pixel> negSubim(*negPsf, negBBox);
    afwImage::Image<double> negModelSubim(negModel, negBBox);
    negModelSubim += negSubim;

    // Portion of the positive Psf that overlaps the model
    int posXmin = std::max(posPsfBBox.getMinX(), posModelBBox.getMinX());
    int posYmin = std::max(posPsfBBox.getMinY(), posModelBBox.getMinY());
    int posXmax = std::min(posPsfBBox.getMaxX(), posModelBBox.getMaxX());
    int posYmax = std::min(posPsfBBox.getMaxY(), posModelBBox.getMaxY());
    geom::Box2I posBBox = geom::Box2I(geom::Point2I(posXmin, posYmin),
                                      geom::Point2I(posXmax, posYmax));
    afwImage::Image<afwMath::Kernel::Pixel> posSubim(*posPsf, posBBox);
    afwImage::Image<double> posModelSubim(posModel, posBBox);
    posModelSubim += posSubim;

    negModel  *= negFlux;   // scale negative model to image
    posModel  *= posFlux;   // scale positive model to image
    afwImage::Image<double> residuals(negModel, true); // full model contains negative lobe...
    residuals += posModel;  // plus positive lobe...
    residuals -= data;      // minus the data...
    residuals *= residuals; // squared...
    residuals /= var;       // divided by the variance : [(model-data)/sigma]**2
    afwMath::Statistics stats = afwMath::makeStatistics(residuals, afwMath::SUM | afwMath::NPOINT);
    double chi2 = stats.getValue(afwMath::SUM);
    int nPix = stats.getValue(afwMath::NPOINT);
    return std::pair<double,int>(chi2, nPix);
}

void PsfDipoleFlux::measure(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {

    typedef afw::image::Exposure<float>::MaskedImageT MaskedImageT;

    std::shared_ptr<afw::detection::Footprint const> footprint = source.getFootprint();
    if (!footprint) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError,
                          (boost::format("No footprint for source %d") % source.getId()).str());
    }

    afw::detection::PeakCatalog peakCatalog = afw::detection::PeakCatalog(footprint->getPeaks());

    if (peakCatalog.size() == 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError,
                          (boost::format("No peak for source %d") % source.getId()).str());
    }
    else if (peakCatalog.size() == 1) {
        // No deblending to do
        return;
    }

    // For N>=2, just measure the brightest-positive and brightest-negative
    // peaks.  peakCatalog is automatically ordered by peak flux, with the most
    // positive one (brightest) being first
    afw::detection::PeakRecord const& positivePeak = peakCatalog.front();
    afw::detection::PeakRecord const& negativePeak = peakCatalog.back();

    // Set up fit parameters and param names
    ROOT::Minuit2::MnUserParameters fitPar;

    fitPar.Add((boost::format("P%d")%NEGCENTXPAR).str(), negativePeak.getFx(), _ctrl.stepSizeCoord);
    fitPar.Add((boost::format("P%d")%NEGCENTYPAR).str(), negativePeak.getFy(), _ctrl.stepSizeCoord);
    fitPar.Add((boost::format("P%d")%NEGFLUXPAR).str(), negativePeak.getPeakValue(), _ctrl.stepSizeFlux);
    fitPar.Add((boost::format("P%d")%POSCENTXPAR).str(), positivePeak.getFx(), _ctrl.stepSizeCoord);
    fitPar.Add((boost::format("P%d")%POSCENTYPAR).str(), positivePeak.getFy(), _ctrl.stepSizeCoord);
    fitPar.Add((boost::format("P%d")%POSFLUXPAR).str(), positivePeak.getPeakValue(), _ctrl.stepSizeFlux);

    // Create the minuit object that knows how to minimise our functor
    //
    MinimizeDipoleChi2 minimizerFunc(*this, source, exposure);
    minimizerFunc.setErrorDef(_ctrl.errorDef);

    //
    // tell minuit about it
    //
    ROOT::Minuit2::MnMigrad migrad(minimizerFunc, fitPar);

    //
    // And let it loose
    //
    ROOT::Minuit2::FunctionMinimum min = migrad(_ctrl.maxFnCalls);

    float minChi2 = min.Fval();
    bool const isValid = min.IsValid() && std::isfinite(minChi2);

    if (true || isValid) {              // calculate coeffs even in minuit is unhappy

        /* I need to call chi2 one more time to grab nPix to calculate chi2/dof.
           Turns out that the Minuit operator method has to be const, and the
           measurement _apply method has to be const, so I can't store nPix as a
           private member variable anywhere.  Consted into a corner.
        */
        std::pair<double,int> fit = chi2(source, exposure,
                                         min.UserState().Value(NEGCENTXPAR),
                                         min.UserState().Value(NEGCENTYPAR),
                                         min.UserState().Value(NEGFLUXPAR),
                                         min.UserState().Value(POSCENTXPAR),
                                         min.UserState().Value(POSCENTYPAR),
                                         min.UserState().Value(POSFLUXPAR));
        double evalChi2 = fit.first;
        int nPix = fit.second;

        std::shared_ptr<geom::Point2D> minNegCentroid(new geom::Point2D(min.UserState().Value(NEGCENTXPAR),
                                                            min.UserState().Value(NEGCENTYPAR)));
        source.set(getNegativeKeys().getInstFlux(), min.UserState().Value(NEGFLUXPAR));
        source.set(getNegativeKeys().getInstFluxErr(), min.UserState().Error(NEGFLUXPAR));

        std::shared_ptr<geom::Point2D> minPosCentroid(new geom::Point2D(min.UserState().Value(POSCENTXPAR),
                                                            min.UserState().Value(POSCENTYPAR)));
        source.set(getPositiveKeys().getInstFlux(), min.UserState().Value(POSFLUXPAR));
        source.set(getPositiveKeys().getInstFluxErr(), min.UserState().Error(POSFLUXPAR));

        source.set(_chi2dofKey, evalChi2 / (nPix - minimizerFunc.getNpar()));
        source.set(_negCentroid.getX(), minNegCentroid->getX());
        source.set(_negCentroid.getY(), minNegCentroid->getY());
        source.set(_posCentroid.getX(), minPosCentroid->getX());
        source.set(_posCentroid.getY(), minPosCentroid->getY());
        source.set(_avgCentroid.getX(), 0.5*(minNegCentroid->getX() + minPosCentroid->getX()));
        source.set(_avgCentroid.getY(), 0.5*(minNegCentroid->getY() + minPosCentroid->getY()));

    }
}

void PsfDipoleFlux::fail(afw::table::SourceRecord & measRecord, meas::base::MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}
}}}  // namespace lsst::ip::diffim
