// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 
#ifndef LSST_IP_DIFFIM_DIPOLEALGORITHMS_H
#define LSST_IP_DIFFIM_DIPOLEALGORITHMS_H
//!
// Control/algorithm hierarchy for dipole measurement.
//

#include "lsst/base.h"
#include "lsst/pex/config.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst {
namespace ip {
namespace diffim {

class DipoleCentroidControl;
class DipoleFluxControl;


/**
 *  @brief Intermediate base class for algorithms that compute a centroid.
 */
class DipoleCentroidAlgorithm : public meas::algorithms::Algorithm {
public:

    /**
     *  @brief Tuple type that holds the keys that define a standard centroid algorithm.
     *
     *  Algorithms are encouraged to add additional flags as appropriate, but these are required.
     */
    typedef afw::table::KeyTuple<afw::table::Centroid> KeyTuple;

    /// @copydoc meas::algorithms::Algorithm::getControl
    DipoleCentroidControl const & getControl() const;

    /// @brief Return the standard centroid keys registered by this algorithm.
    KeyTuple const & getPositiveKeys() const { return _positiveKeys; }
    KeyTuple const & getNegativeKeys() const { return _negativeKeys; }

protected:

    /// @brief Initialize with a manually-constructed key tuple.
    DipoleCentroidAlgorithm(DipoleCentroidControl const & ctrl, KeyTuple const & positiveKeys,
        KeyTuple const & negativeKeys);

    /// @brief Initialize using afw::table::addCentroid field to fill out repetitive descriptions.
    DipoleCentroidAlgorithm(DipoleCentroidControl const & ctrl, 
                            afw::table::Schema & schema, char const * doc);

private:
    KeyTuple _positiveKeys;
    KeyTuple _negativeKeys;
};

class DipoleCentroidControl : public meas::algorithms::AlgorithmControl {
public:

    PTR(DipoleCentroidControl) clone() const {
        return boost::static_pointer_cast<DipoleCentroidControl>(_clone());
    }

    PTR(DipoleCentroidAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)(),
        meas::algorithms::AlgorithmMap const & others = meas::algorithms::AlgorithmMap()
    ) const {
        return boost::static_pointer_cast<DipoleCentroidAlgorithm>(
            _makeAlgorithm(schema, metadata, others));
    }

protected:
    explicit DipoleCentroidControl(std::string const & name_, double priority=0.0) :
        AlgorithmControl(name_, priority) {}
};


/**
 *  @brief Intermediate base class for algorithms that compute a flux.
 */
class DipoleFluxAlgorithm : public meas::algorithms::Algorithm {
public:

    /**
     *  @brief Tuple type that holds the keys that define a standard flux algorithm.
     *
     *  Algorithms are encouraged to add additional flags as appropriate, but these are required.
     */
    typedef afw::table::KeyTuple<afw::table::Flux> KeyTuple;

    /// @copydoc meas::algorithms::Algorithm::getControl
    DipoleFluxControl const & getControl() const;

    /// @brief Return the standard flux keys registered by this algorithm.
    KeyTuple const & getPositiveKeys() const { return _positiveKeys; }
    KeyTuple const & getNegativeKeys() const { return _negativeKeys; }

protected:

    /// @brief Initialize with a manually-constructed key tuple.
    DipoleFluxAlgorithm(DipoleFluxControl const & ctrl, KeyTuple const & positiveKeys,
                        KeyTuple const & negativeKeys);

    /// @brief Initialize using afw::table::addFlux field to fill out repetitive descriptions.
    DipoleFluxAlgorithm(DipoleFluxControl const & ctrl, 
                        afw::table::Schema & schema, char const * doc);

private:
    KeyTuple _positiveKeys;
    KeyTuple _negativeKeys;
};

class DipoleFluxControl : public meas::algorithms::AlgorithmControl {
public:

    PTR(DipoleFluxControl) clone() const { 
        return boost::static_pointer_cast<DipoleFluxControl>(_clone()); }

    PTR(DipoleFluxAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)(),
        meas::algorithms::AlgorithmMap const & others = meas::algorithms::AlgorithmMap()
    ) const {
        return boost::static_pointer_cast<DipoleFluxAlgorithm>(
            _makeAlgorithm(schema, metadata, others));
    }

protected:
    explicit DipoleFluxControl(std::string const & name_, double priority=2.0) :
        AlgorithmControl(name_, priority) {}
};

inline DipoleCentroidAlgorithm::DipoleCentroidAlgorithm(
    DipoleCentroidControl const & ctrl, KeyTuple const & positiveKeys, KeyTuple const & negativeKeys
    ) :
    Algorithm(ctrl), _positiveKeys(positiveKeys), _negativeKeys(negativeKeys)
{}

inline DipoleCentroidAlgorithm::DipoleCentroidAlgorithm(
    DipoleCentroidControl const & ctrl, afw::table::Schema & schema, char const * doc
    ) :
    Algorithm(ctrl), _positiveKeys(afw::table::addCentroidFields(schema, ctrl.name + ".pos", doc)),
    _negativeKeys(afw::table::addCentroidFields(schema, ctrl.name + ".neg", doc))
{}

inline DipoleCentroidControl const & DipoleCentroidAlgorithm::getControl() const {
    return static_cast<DipoleCentroidControl const &>(Algorithm::getControl());
}

inline DipoleFluxAlgorithm::DipoleFluxAlgorithm(
    DipoleFluxControl const & ctrl, KeyTuple const & positiveKeys, KeyTuple const & negativeKeys
    ) :
    Algorithm(ctrl), _positiveKeys(positiveKeys), _negativeKeys(negativeKeys)
{}

inline DipoleFluxAlgorithm::DipoleFluxAlgorithm(
    DipoleFluxControl const & ctrl, afw::table::Schema & schema, char const * doc
    ) :
    Algorithm(ctrl), _positiveKeys(afw::table::addFluxFields(schema, ctrl.name + ".pos", doc)),
    _negativeKeys(afw::table::addFluxFields(schema, ctrl.name + ".neg", doc))
{}

inline DipoleFluxControl const & DipoleFluxAlgorithm::getControl() const {
    return static_cast<DipoleFluxControl const &>(Algorithm::getControl());
}




/**
 *  @brief C++ control object for naive dipole centroid.
 */
class NaiveDipoleCentroidControl : public DipoleCentroidControl {
public:
    NaiveDipoleCentroidControl() : DipoleCentroidControl("centroid.dipole.naive") {}

private:
    virtual PTR(meas::algorithms::AlgorithmControl) _clone() const;
    virtual PTR(meas::algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

/**
 *  @brief C++ control object for naive dipole fluxes.
 */
class NaiveDipoleFluxControl : public DipoleFluxControl {
public:
    NaiveDipoleFluxControl() : DipoleFluxControl("flux.dipole.naive") {}

private:
    virtual PTR(meas::algorithms::AlgorithmControl) _clone() const;
    virtual PTR(meas::algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

/**
 *  @brief C++ control object for PSF dipole fluxes.
 */
class PsfDipoleFluxControl : public DipoleFluxControl {
public:
    LSST_CONTROL_FIELD(maxPixels, int, "Maximum number of pixels to apply the measurement to");

    PsfDipoleFluxControl() : DipoleFluxControl("flux.dipole.psf"), maxPixels(500) {}

private:
    virtual PTR(meas::algorithms::AlgorithmControl) _clone() const;
    virtual PTR(meas::algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};


}}}// namespace lsst::ip::diffim

#endif // !LSST_IP_DIFFIM_DIPOLEALGORITHMS_H
