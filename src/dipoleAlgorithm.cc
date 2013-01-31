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

namespace {
    struct CmpPeakFlux : public std::binary_function<afw::detection::Peak::Ptr,
                                                     afw::detection::Peak::Ptr,
                                                     bool> { 
        bool operator()(afw::detection::Peak::Ptr p1, afw::detection::Peak::Ptr p2) {
            return p1->getPeakValue() < p2->getPeakValue();
        }
    };
} // end anonymous
 
template<typename PixelT>
void PsfDipoleFlux::_apply(
    afw::table::SourceRecord & source, 
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center // Not used; source required to have footprint&peaks
) const {
    source.set(getPositiveKeys().flag, true); // say we've failed so that's the result if we throw
    source.set(getNegativeKeys().flag, true); // say we've failed so that's the result if we throw

    //float background = static_cast<PsfDipoleFluxControl const &>(getControl()).background;
    typedef typename afw::image::Exposure<PixelT>::MaskedImageT MaskedImageT;

    CONST_PTR(afw::detection::Footprint) footprint = source.getFootprint();
    if (!footprint) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException,
                          (boost::format("No footprint for source %d") % source.getId()).str());
    }
    afw::detection::Footprint::PeakList peakList = afw::detection::Footprint::PeakList(footprint->getPeaks());

    if (peakList.size() == 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException,
                          (boost::format("No peak for source %d") % source.getId()).str());
    }
    else if (peakList.size() == 1) {
        std::cout << "TESTING: Only 1 peak! " << std::endl;
        // No deblending to do 
        return;
    }
    // For N>=2, just measure the brightest and faintest peaks
    
    // Order by peak flux, most negative one is first
    std::sort(peakList.begin(), peakList.end(), CmpPeakFlux());
    PTR(afw::detection::Peak) negativePeak = peakList.front();
    PTR(afw::detection::Peak) positivePeak = peakList.back();

    afw::geom::Point2D negCenter(negativePeak->getFx(), negativePeak->getFy());
    afw::geom::Point2D posCenter(positivePeak->getFx(), positivePeak->getFy());

    /* 1) Naive implementation: fit each peak separately using the above
     * centroid
     *
     */
    
    /* Minimal source to receive the Psf flux */
    afwTable::Schema schema = afwTable::SourceTable::makeMinimalSchema();
    measAlgorithms::PsfFluxControl fluxControl;
    measAlgorithms::MeasureSources ms =
        measAlgorithms::MeasureSourcesBuilder()
        .addAlgorithm(fluxControl)
        .build(schema);
    PTR(afwTable::SourceTable) table = afwTable::SourceTable::make(schema);

    PTR(afwTable::SourceRecord) negativeSource = table->makeRecord();
    ms.apply(*negativeSource, exposure, negCenter, false);

    PTR(afwTable::SourceRecord) positiveSource = table->makeRecord();
    ms.apply(*positiveSource, exposure, posCenter, false);

    afwTable::Flux::MeasKey fluxKey = table->getSchema()[fluxControl.name];
    std::cout << "TESTINGA: Negative flux " << negativeSource->get(fluxKey) << std::endl;
    std::cout << "TESTINGA: Positive flux " << positiveSource->get(fluxKey) << std::endl;

    /* 2) Slightly better implementation: fit for the superposition of Psfs at
     * the two centroids.  Still does not take into account possible centroid
     * offsets.
     *
     */
    try {
        CONST_PTR(afwDet::Psf) psf = exposure.getPsf();
        PTR(afwImage::Image<afwMath::Kernel::Pixel>) negPsf = psf->computeImage(negCenter, true);
        PTR(afwImage::Image<afwMath::Kernel::Pixel>) posPsf = psf->computeImage(posCenter, true);
        double negSum = negPsf->getArray().asEigen().sum();
        double posSum = posPsf->getArray().asEigen().sum();
        
        afwImage::Image<double> negModel(footprint->getBBox());
        afwImage::Image<double> posModel(footprint->getBBox());
        afwImage::Image<PixelT> data(*(exposure.getMaskedImage().getImage()), footprint->getBBox(), afwImage::PARENT);
        
        afwGeom::Box2I negPsfBBox = negPsf->getBBox(afwImage::PARENT);
        afwGeom::Box2I posPsfBBox = posPsf->getBBox(afwImage::PARENT);
        afwGeom::Box2I negModelBBox = negModel.getBBox(afwImage::PARENT);
        afwGeom::Box2I posModelBBox = posModel.getBBox(afwImage::PARENT);
        
        // Portion of the negative Psf that overlaps the model
        int negXmin = std::max(negPsfBBox.getMinX(), negModelBBox.getMinX());
        int negYmin = std::max(negPsfBBox.getMinY(), negModelBBox.getMinY());
        int negXmax = std::min(negPsfBBox.getMaxX(), negModelBBox.getMaxX());
        int negYmax = std::min(negPsfBBox.getMaxY(), negModelBBox.getMaxY());
        afwGeom::Box2I negBBox = afwGeom::Box2I(afwGeom::Point2I(negXmin, negYmin), 
                                                afwGeom::Point2I(negXmax, negYmax));
        afwImage::Image<afwMath::Kernel::Pixel> negSubim(*negPsf, negBBox, afwImage::PARENT);
        afwImage::Image<double> negModelSubim(negModel, negBBox, afwImage::PARENT);
        negModelSubim += negSubim;
        
        
        // Portion of the positive Psf that overlaps the model
        int posXmin = std::max(posPsfBBox.getMinX(), posModelBBox.getMinX());
        int posYmin = std::max(posPsfBBox.getMinY(), posModelBBox.getMinY());
        int posXmax = std::min(posPsfBBox.getMaxX(), posModelBBox.getMaxX());
        int posYmax = std::min(posPsfBBox.getMaxY(), posModelBBox.getMaxY());
        afwGeom::Box2I posBBox = afwGeom::Box2I(afwGeom::Point2I(posXmin, posYmin), 
                                                afwGeom::Point2I(posXmax, posYmax));
        afwImage::Image<afwMath::Kernel::Pixel> posSubim(*posPsf, posBBox, afwImage::PARENT);
        afwImage::Image<double> posModelSubim(posModel, posBBox, afwImage::PARENT);
        posModelSubim += posSubim;
        
        // Set up a linear least squares fit with no centroid shift
        ndarray::Array<double, 2, 2> Mt = ndarray::allocate(2, footprint->getArea());
        ndarray::Array<double, 1, 1> b = ndarray::allocate(footprint->getArea());
        
        afwDet::flattenArray(*footprint, negModel.getArray(), Mt[0].shallow(), negModel.getXY0());
        afwDet::flattenArray(*footprint, posModel.getArray(), Mt[1].shallow(), posModel.getXY0());
        afwDet::flattenArray(*footprint, data.getArray(), b, data.getXY0());
        afw::math::LeastSquares lstsq = afwMath::LeastSquares::fromDesignMatrix(Mt.transpose().shallow(), b);
        double fluxNeg = lstsq.getSolution()[0];
        double fluxPos = lstsq.getSolution()[1];
        std::cout << "TESTINGB: Negative flux " << fluxNeg << " " << negSum << " " << fluxNeg*negSum << std::endl;
        std::cout << "TESTINGB: Positive flux " << fluxPos << " " << posSum << " " << fluxPos*posSum << std::endl;
        
        negModel.writeFits("/tmp/neg.fits");
        posModel.writeFits("/tmp/pos.fits");
        data.writeFits("/tmp/data.fits");
        negModel *= fluxNeg;
        posModel *= fluxPos;
        negModel += posModel;
        negModel.writeFits("/tmp/model.fits");
        negModel -= data;
        negModel.writeFits("/tmp/resids.fits");
                
    } catch (pex::exceptions::Exception const& e) {
        // Swallow all exceptions, because one bad measurement shouldn't affect all others
        std::cout << str(boost::format("Measuring on source %d at (%f,%f): %s") %
                         source.getId() % center.getX() % center.getY() % e.what()) << std::endl;
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
