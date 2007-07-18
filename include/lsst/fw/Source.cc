// -*- LSST-C++ -*-
/**
 * \file
 *
 * Definition of inline member functions declared in Source.h
 *
 * This file is meant to be included by lsst/fw/Source.h
 *
 * \author Andrew Becker
 *
 * \ingroup fw
 */

// Inline Member Functions
inline double
lsst::fw::Source::getId() const {
    return _id;
}
inline double
lsst::fw::Source::getColc() const {
    return _colc;
}
inline double
lsst::fw::Source::getRowc() const {
    return _rowc;
}
inline double
lsst::fw::Source::getDcol() const {
    return _dcol;
}
inline double
lsst::fw::Source::getDrow() const {
    return _drow;
}
lsst::fw::Source::Source()
:
    lsst::fw::LsstBase(typeid(this)),
    _id(0.),
    _colc(0.),
    _rowc(0.),
    _dcol(0.),
    _drow(0.)
{ }
lsst::fw::Source::Source(
    int id, ///< unique identifier
    double colc, ///< column center
    double rowc, ///< row center
    double dcol, ///< column extent
    double drow) ///< row extent
:
    lsst::fw::LsstBase(typeid(this)),
    _id(id),
    _colc(colc),
    _rowc(rowc),
    _dcol(dcol),
    _drow(drow)
{ }
