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

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/meas/algorithms/CentroidControl.h"

#include "lsst/ip/diffim/DipoleAlgorithms.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;

namespace lsst {
namespace ip {
namespace diffim {


/**
 * @brief A class that knows how to calculate centroids as a simple unweighted first moment
 * of the 3x3 region around the peaks
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
 * @brief A class that knows how to calculate centroids as a simple unweighted first moment
 * of the 3x3 region around the peaks
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
    afw::table::KeyTuple<afw::table::Centroid> keys,
    float background
    )
{
    source.set(keys.meas, afw::geom::Point2D(center));
    typedef afw::image::Image<PixelT> ImageT;
    ImageT const& image = *exposure.getMaskedImage().getImage();

    int x = center.getX() - image.getX0();
    int y = center.getY() - image.getY0();

    if (x < 1 || x >= image.getWidth() - 1 || y < 1 || y >= image.getHeight() - 1) {
         throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                           (boost::format("Object at (%d, %d) is too close to the edge") % x % y).str());
    }

    typename ImageT::xy_locator im = image.xy_at(x, y);

    double const sum =
        (im(-1,  1) + im( 0,  1) + im( 1,  1) +
         im(-1,  0) + im( 0,  0) + im( 1,  0) +
         im(-1, -1) + im( 0, -1) + im( 1, -1))
        - 9 * background;

    if (sum == 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
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
 * @brief Given an image and a pixel position, return a Centroid using a naive 3x3 weighted moment
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
    float background = static_cast<NaiveDipoleCentroidControl const &>(getControl()).background;


    naiveCentroid(source, exposure, peaks[0]->getI(), (peaks[0]->getPeakValue() >= 0  ? getPositiveKeys() :
                                                       getNegativeKeys()), background);
    if (peaks.size() > 1) {
        naiveCentroid(source, exposure, peaks[1]->getI(), (peaks[1]->getPeakValue() >= 0  ? getPositiveKeys() :
                                                           getNegativeKeys()), background);
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
    explicit NaiveDipoleFootprinter(MaskedImageT const& mimage, ///< The image the source lives in
                                    float background            ///< Background tweak value
        ) : afw::detection::FootprintFunctor<MaskedImageT>(mimage), _background(background),
            _sumPositive(0.0), _sumNegative(0.0), _numPositive(0), _numNegative(0) {}

    /// @brief Reset everything for a new Footprint
    void reset() {
        _sumPositive = _sumNegative = 0.0;
        _numPositive = _numNegative = 0;
    }
    void reset(afwDet::Footprint const&) {}

    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int,                                   ///< column-position of pixel
                    int                                    ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel ival = loc.image(0, 0) - _background;
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
    float _background;
    double _sumPositive;
    double _sumNegative;
    double _varPositive;
    double _varNegative;
    int _numPositive;
    int _numNegative;
};

} // anonymous namespace



/**
 * @brief Given an image and a pixel position, return a Centroid using a naive 3x3 weighted moment
 */
template<typename PixelT>
void NaiveDipoleFlux::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getPositiveKeys().flag, true); // say we've failed so that's the result if we throw
    source.set(getNegativeKeys().flag, true); // say we've failed so that's the result if we throw

    float background = static_cast<NaiveDipoleFluxControl const &>(getControl()).background;
    typedef typename afw::image::Exposure<PixelT>::MaskedImageT MaskedImageT;

    NaiveDipoleFootprinter<MaskedImageT> functor(exposure.getMaskedImage(), background);
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
 * @brief Implementation of Psf dipole flux
 */
class PsfDipoleFlux : public DipoleFluxAlgorithm {
public:

    PsfDipoleFlux(PsfDipoleFluxControl const & ctrl, afw::table::Schema & schema) :
        DipoleFluxAlgorithm(ctrl, schema, "raw psf flux counts")
    {}

private:
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(PsfDipoleFlux);

};

template<typename PixelT>
void PsfDipoleFlux::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center // Not used; source required to have footprint&peaks
) const {
    source.set(getPositiveKeys().flag, true); // say we've failed so that's the result if we throw
    source.set(getNegativeKeys().flag, true); // say we've failed so that's the result if we throw

    float background = static_cast<PsfDipoleFluxControl const &>(getControl()).background;
    typedef typename afw::image::Exposure<PixelT>::MaskedImageT MaskedImageT;

    CONST_PTR(afw::detection::Footprint) foot = source.getFootprint();
    if (!foot) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException,
                          (boost::format("No footprint for source %d") % source.getId()).str());
    }
    afw::detection::Footprint::PeakList const& peakList = foot->getPeaks();
    if (peakList.size() == 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException,
                          (boost::format("No peak for source %d") % source.getId()).str());
    }

    PTR(afw::detection::Peak) peak1 = peakList[0];
    afw::geom::Point2D center1(peak1->getFx(), peak1->getFy());
    PTR(afw::detection::Peak) peak2 = peakList[1];
    afw::geom::Point2D center2(peak2->getFx(), peak2->getFy());
    
    /*
    afwTable::Schema schema = afwTable::SourceTable::makeMinimalSchema();
    measAlgorithms::PsfFluxControl fluxControlPos(_ctrl.name+".pos");
    measAlgorithms::MeasureSources msPos =
        measAlgorithms::MeasureSourcesBuilder()
        .addAlgorithm(fluxControlPos)
        .build(schema);

    afwTable::Schema schema = afwTable::SourceTable::makeMinimalSchema();
    measAlgorithms::PsfFluxControl fluxControlNeg(_ctrl.name+".neg");
    measAlgorithms::MeasureSources msNeg =
        measAlgorithms::MeasureSourcesBuilder()
        .addAlgorithm(fluxControlNeg)
        .build(schema);

    msPos.applyWithPixel(sourcePos, exposure);
    msNeg.applyWithPixel(sourceNeg, exposure);
    */

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
