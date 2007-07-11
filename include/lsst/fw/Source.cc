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
inline unsigned
lsst::fw::Source::getColc() const {
    return _colc;
}
inline unsigned
lsst::fw::Source::getRowc() const {
    return _rowc;
}
inline unsigned
lsst::fw::Source::getDcol() const {
    return _dcol;
}
inline unsigned
lsst::fw::Source::getDrow() const {
    return _drow;
}
lsst::fw::Source::Source()
:
    lsst::fw::LsstBase(typeid(this)),
    _rowc(0),
    _colc(0),
    _drow(0),
    _dcol(0)
{ }
lsst::fw::Source::Source(
    unsigned rowc, ///< row center
    unsigned colc, ///< column center
    unsigned drow, ///< row extent
    unsigned dcol) ///< column extent
:
    lsst::fw::LsstBase(typeid(this)),
    _rowc(rowc),
    _colc(colc),
    _drow(drow),
    _dcol(dcol)
{ }
