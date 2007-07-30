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
    _drow(0.),
    // Below is LSST DB Schema
    _sourceId(0),
    _ampExposureId(0),
    _filterId(0),
    _objectId(0),
    _movingObjectId(0),
    _procHistoryId(0),
    _ra(0),
    _decl(0),
    ___zoneId_placeholder(0),
    _raErr4wcs(0),
    _decErr4wcs(0),
    _raErr4detection(0),
    _decErr4detection(0),
    _rowcErr(0),
    _colcErr(0),
    _cx(0),
    _cy(0),
    _cz(0),
    _taiMidPoint(0),
    _taiRange(0),
    _fwhmA(0),
    _fwhmB(0),
    _fwhmTheta(0),
    _flux(0),
    _fluxErr(0),
    _psfMag(0),
    _psfMagErr(0),
    _apMag(0),
    _apMagErr(0),
    _apDia(0),
    _petroMag(0),
    _petroMagErr(0),
    _snr(0),
    _chi2(0),
    _sky(0),
    _skyErr(0),
    _moment0(0),
    _moment1_x(0),
    _moment1_y(0),
    _moment2_xx(0),
    _moment2_xy(0),
    _moment2_yy(0),
    _moment3_xxx(0),
    _moment3_xxy(0),
    _moment3_xyy(0),
    _moment3_yyy(0),
    _moment4_xxxx(0),
    _moment4_xxxy(0),
    _moment4_xxyy(0),
    _moment4_xyyy(0),
    _moment4_yyyy(0),
    _flag4association(0),
    _flag4detection(0),
    _flag4wcs(0)
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
