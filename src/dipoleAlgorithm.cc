// -*- LSST-C++ -*-

/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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
#include <functional>   // std::binary_function
#include <limits>       // std::numeric_limits
#include <cmath>        // std::sqrt

#if !defined(DOXYGEN)
#   include "Minuit2/FCNBase.h"
#   include "Minuit2/FunctionMinimum.h"
#   include "Minuit2/MnMigrad.h"
#   include "Minuit2/MnMinos.h"
#   include "Minuit2/MnPrint.h"
#endif

#include "boost/shared_ptr.hpp"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/afw/table.h"
#include "lsst/afw/math.h"
#include "lsst/afw/geom.h"
#include "lsst/meas/algorithms.h"
#include "lsst/ip/diffim/DipoleAlgorithms.h"
#include "ndarray/eigen.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwTable = lsst::afw::table;
namespace afwMath = lsst::afw::math;
namespace afwGeom = lsst::afw::geom;
namespace measAlgorithms = lsst::meas::algorithms;

namespace lsst {
namespace ip {
namespace diffim {

    int const NEGCENTXPAR(0); // Parameter for the x-component of the negative lobe centroid
    int const NEGCENTYPAR(1); // Parameter for the y-component of the negative lobe centroid 
    int const NEGFLUXPAR(2);  // Parameter for the flux of the negative lobe
    int const POSCENTXPAR(3); // Parameter for the x-component of the positive lobe centroid
    int const POSCENTYPAR(4); // Parameter for the y-component of the positive lobe centroid
    int const POSFLUXPAR(5);  // Parameter for the flux of the positive lobe

/**
 * A class that knows how to calculate centroids as a simple unweighted first
 * moment of the 3x3 region around the peaks
 */
class NaiveDipoleCentroid : public DipoleCentroidAlgorithm {
public:

    NaiveDipoleCentroid(NaiveDipoleCentroidControl const & ctrl, afw::table::Schema & schema) :
        DipoleCentroidAlgorithm(ctrl, schema, "unweighted 3x3 first moment centroid")
    {}

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(NaiveDipoleCentroid);
};

/**
 * A class that knows how to calculate centroids as a simple unweighted first
 * moment of the 3x3 region around the peaks
 */
class NaiveDipoleFlux : public DipoleFluxAlgorithm {
public:

    NaiveDipoleFlux(NaiveDipoleFluxControl const & ctrl, afw::table::Schema & schema) :
        DipoleFluxAlgorithm(ctrl, schema, "raw flux counts"),
        _numPositiveKey(schema.addField<int>(ctrl.name+".npos", "number of positive pixels", "dn")),
        _numNegativeKey(schema.addField<int>(ctrl.name+".nneg", "number of negative pixels", "dn"))
    {}

private:
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(NaiveDipoleFlux);

    afw::table::Key<int> _numPositiveKey;
    afw::table::Key<int> _numNegativeKey;
};



namespace {

template<typename PixelT>
void naiveCentroid(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2I const & center,
    afw::table::KeyTuple<afw::table::Centroid> keys
    )
{
    source.set(keys.meas, afw::geom::Point2D(center));
    typedef afw::image::Image<PixelT> ImageT;
    ImageT const& image = *exposure.getMaskedImage().getImage();

    int x = center.getX() - image.getX0();
    int y = center.getY() - image.getY0();

    if (x < 1 || x >= image.getWidth() - 1 || y < 1 || y >= image.getHeight() - 1) {
         throw LSST_EXCEPT(lsst::pex::exceptions::LengthError,
                           (boost::format("Object at (%d, %d) is too close to the edge") 
                            % x % y).str());
    }

    typename ImageT::xy_locator im = image.xy_at(x, y);

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

    source.set(keys.flag, false);
    source.set(
        keys.meas, 
        afw::geom::Point2D(
            lsst::afw::image::indexToPosition(x + image.getX0()) + sum_x / sum,
            lsst::afw::image::indexToPosition(y + image.getY0()) + sum_y / sum
        )
    );
}

} // anonymous namespace


/**
 * Given an image and a pixel position, return a Centroid using a naive 3x3 weighted moment
 */
template<typename PixelT>
void NaiveDipoleCentroid::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getPositiveKeys().flag, true); // say we've failed so that's the result if we throw
    source.set(getNegativeKeys().flag, true); // say we've failed so that's the result if we throw

    afw::detection::Footprint::PeakList const& peaks = source.getFootprint()->getPeaks();

    naiveCentroid(source, exposure, peaks[0]->getI(), (peaks[0]->getPeakValue() >= 0 ? 
                                                       getPositiveKeys() : 
                                                       getNegativeKeys()));
    if (peaks.size() > 1) {
        naiveCentroid(source, exposure, peaks[1]->getI(), (peaks[1]->getPeakValue() >= 0 ? 
                                                           getPositiveKeys() :
                                                           getNegativeKeys()));
    }
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(NaiveDipoleCentroid);

PTR(meas::algorithms::AlgorithmControl) NaiveDipoleCentroidControl::_clone() const {
    return boost::make_shared<NaiveDipoleCentroidControl>(*this);
}

PTR(meas::algorithms::Algorithm) NaiveDipoleCentroidControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<NaiveDipoleCentroid>(*this, boost::ref(schema));
}



namespace {

template <typename MaskedImageT>
class NaiveDipoleFootprinter : public afw::detection::FootprintFunctor<MaskedImageT> {
public:
    explicit NaiveDipoleFootprinter(MaskedImageT const& mimage ///< The image the source lives in
        ) : afw::detection::FootprintFunctor<MaskedImageT>(mimage), 
            _sumPositive(0.0), _sumNegative(0.0), _numPositive(0), _numNegative(0) {}

    /// Reset everything for a new Footprint
    void reset() {
        _sumPositive = _sumNegative = 0.0;
        _numPositive = _numNegative = 0;
    }
    void reset(afwDet::Footprint const&) {}

    /// method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int,                                   ///< column-position of pixel
                    int                                    ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel ival = loc.image(0, 0);
        typename MaskedImageT::Image::Pixel vval = loc.variance(0, 0);
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
template<typename PixelT>
void NaiveDipoleFlux::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getPositiveKeys().flag, true); // say we've failed so that's the result if we throw
    source.set(getNegativeKeys().flag, true); // say we've failed so that's the result if we throw

    typedef typename afw::image::Exposure<PixelT>::MaskedImageT MaskedImageT;

    NaiveDipoleFootprinter<MaskedImageT> functor(exposure.getMaskedImage());
    functor.apply(*source.getFootprint());

    source.set(getPositiveKeys().meas, functor.getSumPositive());
    source.set(getPositiveKeys().err, ::sqrt(functor.getVarPositive()));
    source.set(_numPositiveKey, functor.getNumPositive());
    source.set(getPositiveKeys().flag, false);

    source.set(getNegativeKeys().meas, functor.getSumNegative());
    source.set(getNegativeKeys().err, ::sqrt(functor.getVarNegative()));
    source.set(_numNegativeKey, functor.getNumNegative());
    source.set(getNegativeKeys().flag, false);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(NaiveDipoleFlux);

PTR(meas::algorithms::AlgorithmControl) NaiveDipoleFluxControl::_clone() const {
    return boost::make_shared<NaiveDipoleFluxControl>(*this);
}

PTR(meas::algorithms::Algorithm) NaiveDipoleFluxControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<NaiveDipoleFlux>(*this, boost::ref(schema));
}



/**
 * Implementation of Psf dipole flux
 */
class PsfDipoleFlux : public DipoleFluxAlgorithm {
public:

    PsfDipoleFlux(PsfDipoleFluxControl const & ctrl, afw::table::Schema & schema) :
        DipoleFluxAlgorithm(ctrl, schema, "jointly fitted psf flux counts"),
        _chi2dofKey(schema.addField<float>(ctrl.name+".chi2dof", 
                                           "chi2 per degree of freedom of fit")),
        _avgCentroid(
            addCentroidFields(schema, ctrl.name+".centroid", 
                              "average of the postive and negative lobe positions")),
        _negCentroid(
            addCentroidFields(schema, ctrl.name+".neg.centroid", 
                              "psf fitted center of negative lobe")),
        _posCentroid(
            addCentroidFields(schema, ctrl.name+".pos.centroid", 
                              "psf fitted center of positive lobe")),
        _flagMaxPixelsKey(schema.addField<afw::table::Flag>(ctrl.name+".flags.maxpix",
                                                            "set if too large a footprint was sent to the algorithm"))
    {}
    template <typename PixelT>
    std::pair<double,int> chi2(afw::table::SourceRecord & source,
                afw::image::Exposure<PixelT> const & exposure,
                double negCenterX, double negCenterY, double negFlux,
                double posCenterX, double poCenterY, double posFlux
                ) const;


private:
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(PsfDipoleFlux);

    afw::table::Key<float> _chi2dofKey;
    afw::table::KeyTuple<afw::table::Centroid> _avgCentroid;
    afw::table::KeyTuple<afw::table::Centroid> _negCentroid;
    afw::table::KeyTuple<afw::table::Centroid> _posCentroid;

    afw::table::Key< afw::table::Flag > _flagMaxPixelsKey;
};

/**
 * Class to minimize PsfDipoleFlux; this is the object that Minuit minimizes
 */
template<typename PixelT>
class MinimizeDipoleChi2 : public ROOT::Minuit2::FCNBase {
public:
    explicit MinimizeDipoleChi2(PsfDipoleFlux const& psfDipoleFlux,
                                afw::table::SourceRecord & source,
                                afw::image::Exposure<PixelT> const& exposure
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

        std::pair<double,int> fit = _psfDipoleFlux.chi2(_source, _exposure, negCenterX, negCenterY, negFlux, posCenterX, posCenterY, posFlux);
        double chi2 = fit.first;
        int nPix = fit.second;
        if (nPix > _maxPix) {
            return _bigChi2;
        }

        return chi2;
    }
    
private:
    double _errorDef;               // how much cost function has changed at the +- 1 error points
    int _nPar;                      // number of parameters in the fit; hard coded for MinimizeDipoleChi2
    int _maxPix;                    // maximum number of pixels that shoud be in the footprint; prevents too much centroid wander
    double _bigChi2;                // large value to tell fitter when it has gone into bad region of parameter space
    
    PsfDipoleFlux const& _psfDipoleFlux;
    afw::table::SourceRecord & _source;
    afw::image::Exposure<PixelT> const& _exposure;
};

template<typename PixelT>
std::pair<double,int> PsfDipoleFlux::chi2(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    double negCenterX, double negCenterY, double negFlux,
    double posCenterX, double posCenterY, double posFlux
) const { 
    
    afw::geom::Point2D negCenter(negCenterX, negCenterY);
    afw::geom::Point2D posCenter(posCenterX, posCenterY);

    CONST_PTR(afw::detection::Footprint) footprint = source.getFootprint();

    /* Fit for the superposition of Psfs at the two centroids.  
     *
     */
    CONST_PTR(afwDet::Psf) psf = exposure.getPsf();
    PTR(afwImage::Image<afwMath::Kernel::Pixel>) negPsf = psf->computeImage(negCenter);
    PTR(afwImage::Image<afwMath::Kernel::Pixel>) posPsf = psf->computeImage(posCenter);
    
    afwImage::Image<double> negModel(footprint->getBBox());
    afwImage::Image<double> posModel(footprint->getBBox());
    afwImage::Image<PixelT> data(*(exposure.getMaskedImage().getImage()), 
                                 footprint->getBBox());
    afwImage::Image<afwImage::VariancePixel> var(*(exposure.getMaskedImage().getVariance()), 
                                                 footprint->getBBox());
    
    afwGeom::Box2I negPsfBBox = negPsf->getBBox();
    afwGeom::Box2I posPsfBBox = posPsf->getBBox();
    afwGeom::Box2I negModelBBox = negModel.getBBox();
    afwGeom::Box2I posModelBBox = posModel.getBBox();
    
    // Portion of the negative Psf that overlaps the model
    int negXmin = std::max(negPsfBBox.getMinX(), negModelBBox.getMinX());
    int negYmin = std::max(negPsfBBox.getMinY(), negModelBBox.getMinY());
    int negXmax = std::min(negPsfBBox.getMaxX(), negModelBBox.getMaxX());
    int negYmax = std::min(negPsfBBox.getMaxY(), negModelBBox.getMaxY());
    afwGeom::Box2I negBBox = afwGeom::Box2I(afwGeom::Point2I(negXmin, negYmin), 
                                            afwGeom::Point2I(negXmax, negYmax));
    afwImage::Image<afwMath::Kernel::Pixel> negSubim(*negPsf, negBBox);
    afwImage::Image<double> negModelSubim(negModel, negBBox);
    negModelSubim += negSubim;
    
    // Portion of the positive Psf that overlaps the model
    int posXmin = std::max(posPsfBBox.getMinX(), posModelBBox.getMinX());
    int posYmin = std::max(posPsfBBox.getMinY(), posModelBBox.getMinY());
    int posXmax = std::min(posPsfBBox.getMaxX(), posModelBBox.getMaxX());
    int posYmax = std::min(posPsfBBox.getMaxY(), posModelBBox.getMaxY());
    afwGeom::Box2I posBBox = afwGeom::Box2I(afwGeom::Point2I(posXmin, posYmin), 
                                            afwGeom::Point2I(posXmax, posYmax));
    afwImage::Image<afwMath::Kernel::Pixel> posSubim(*posPsf, posBBox);
    afwImage::Image<double> posModelSubim(posModel, posBBox);
    posModelSubim += posSubim;
    
    negModel  *= negFlux;  // scale negative model to image
    posModel  *= posFlux;  // scale positive model to image
    afwImage::Image<double> residuals(negModel, true); // full model contains negative lobe...
    residuals += posModel; // plus positive lobe...
    residuals -= data;     // minus the data...
    residuals *= residuals;// squared...
    residuals /= var;      // divided by the variance : [(model-data)/sigma]**2
    afwMath::Statistics stats = afwMath::makeStatistics(residuals, afwMath::SUM | afwMath::NPOINT);
    double chi2 = stats.getValue(afwMath::SUM);
    int nPix = stats.getValue(afwMath::NPOINT);
    return std::pair<double,int>(chi2, nPix);
}    
 
template<typename PixelT>
void PsfDipoleFlux::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center // Not used; source required to have footprint&peaks
) const {

    source.set(getPositiveKeys().flag, true); // say we've failed so that's the result if we throw
    source.set(getNegativeKeys().flag, true); // say we've failed so that's the result if we throw
    source.set(_flagMaxPixelsKey, true);

    typedef typename afw::image::Exposure<PixelT>::MaskedImageT MaskedImageT;

    CONST_PTR(afw::detection::Footprint) footprint = source.getFootprint();
    if (!footprint) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError,
                          (boost::format("No footprint for source %d") % source.getId()).str());
    }

    PsfDipoleFluxControl const & ctrl = static_cast<PsfDipoleFluxControl const &>(getControl());
    if (footprint->getArea() > ctrl.maxPixels) {
        // Too big
        return;
    }
    source.set(_flagMaxPixelsKey, false);

    afw::detection::Footprint::PeakList peakList = 
        afw::detection::Footprint::PeakList(footprint->getPeaks());

    if (peakList.size() == 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError,
                          (boost::format("No peak for source %d") % source.getId()).str());
    }
    else if (peakList.size() == 1) {
        // No deblending to do 
        return;
    }

    // For N>=2, just measure the brightest-positive and brightest-negative
    // peaks.  peakList is automatically ordered by peak flux, with the most
    // positive one (brightest) being first
    PTR(afw::detection::Peak) positivePeak = peakList.front();
    PTR(afw::detection::Peak) negativePeak = peakList.back();

    // Set up fit parameters and param names
    ROOT::Minuit2::MnUserParameters fitPar;

    fitPar.Add((boost::format("P%d")%NEGCENTXPAR).str(), negativePeak->getFx(), ctrl.stepSizeCoord);
    fitPar.Add((boost::format("P%d")%NEGCENTYPAR).str(), negativePeak->getFy(), ctrl.stepSizeCoord);
    fitPar.Add((boost::format("P%d")%NEGFLUXPAR).str(), negativePeak->getPeakValue(), ctrl.stepSizeFlux);
    fitPar.Add((boost::format("P%d")%POSCENTXPAR).str(), positivePeak->getFx(), ctrl.stepSizeCoord);
    fitPar.Add((boost::format("P%d")%POSCENTYPAR).str(), positivePeak->getFy(), ctrl.stepSizeCoord);
    fitPar.Add((boost::format("P%d")%POSFLUXPAR).str(), positivePeak->getPeakValue(), ctrl.stepSizeFlux);

    // Create the minuit object that knows how to minimise our functor
    //
    MinimizeDipoleChi2<PixelT> minimizerFunc(*this, source, exposure);
    minimizerFunc.setErrorDef(ctrl.errorDef);

    //
    // tell minuit about it
    //    
    ROOT::Minuit2::MnMigrad migrad(minimizerFunc, fitPar);

    //
    // And let it loose
    //
    ROOT::Minuit2::FunctionMinimum min = migrad(ctrl.maxFnCalls);

    float minChi2 = min.Fval();
    bool const isValid = min.IsValid() && std::isfinite(minChi2);

    if (true || isValid) {              // calculate coeffs even in minuit is unhappy

        /* I need to call chi2 one more time to grab nPix to calculate chi2/dof.
           Turns out that the Minuit operator method has to be const, and the
           measurement _apply method has to be const, so I can't store nPix as a
           private member variable anywhere.  Consted into a corner.
        */
        std::pair<double,int> fit = chi2(source, exposure, 
                                         min.UserState().Value(NEGCENTXPAR), min.UserState().Value(NEGCENTYPAR), 
                                         min.UserState().Value(NEGFLUXPAR), min.UserState().Value(POSCENTXPAR), 
                                         min.UserState().Value(POSCENTYPAR), min.UserState().Value(POSFLUXPAR));
        double evalChi2 = fit.first;
        int nPix = fit.second;
        
        PTR(afw::geom::Point2D) minNegCentroid(new afw::geom::Point2D(min.UserState().Value(NEGCENTXPAR), 
                                                                      min.UserState().Value(NEGCENTYPAR)));
        source.set(getNegativeKeys().meas, min.UserState().Value(NEGFLUXPAR));
        source.set(getNegativeKeys().err, min.UserState().Error(NEGFLUXPAR));
        source.set(getNegativeKeys().flag, false);
        
        PTR(afw::geom::Point2D) minPosCentroid(new afw::geom::Point2D(min.UserState().Value(POSCENTXPAR), 
                                                                      min.UserState().Value(POSCENTYPAR))); 
        source.set(getPositiveKeys().meas, min.UserState().Value(POSFLUXPAR));
        source.set(getPositiveKeys().err, min.UserState().Error(POSFLUXPAR));
        source.set(getPositiveKeys().flag, false);

        source.set(_chi2dofKey, evalChi2 / (nPix - minimizerFunc.getNpar()));
        source.set(_negCentroid.meas, *minNegCentroid);
        source.set(_posCentroid.meas, *minPosCentroid);
        source.set(_avgCentroid.meas, 
                   afw::geom::Point2D(0.5*(minNegCentroid->getX() + minPosCentroid->getX()),
                                      0.5*(minNegCentroid->getY() + minPosCentroid->getY())));

    }

}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(PsfDipoleFlux);

PTR(meas::algorithms::AlgorithmControl) PsfDipoleFluxControl::_clone() const {
    return boost::make_shared<PsfDipoleFluxControl>(*this);
}

PTR(meas::algorithms::Algorithm) PsfDipoleFluxControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<PsfDipoleFlux>(*this, boost::ref(schema));
}



}}}  // namespace lsst::ip::diffim
