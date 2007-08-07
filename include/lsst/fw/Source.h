// -*- LSST-C++ -*-
#ifndef LSST_FW_Source_H
#define LSST_FW_Source_H
/**
 * \file
 *
 * A VERY beta version of Source to enable image subtraction testing
 *
 * \author Andrew Becker
 *
 * \ingroup fw
 */

namespace lsst {
namespace fw {

    class Source : private lsst::mwi::LsstBase {

    public:
        explicit Source();
        virtual ~Source() {};
        explicit Source(
            long int sourceId,
            double colc,
            double rowc,
            double dcol,
            double drow
            );
       inline double getSourceId() const;
       inline double getColc() const;
       inline double getRowc() const;
       inline double getDcol() const;
       inline double getDrow() const;
    private:
       long int _sourceId;
       double _colc;
       double _rowc;
       double _dcol;
       double _drow;

       // Below is LSST DB Schema
       long int _ampExposureId;
       short int _filterId;
       long int _objectId;
       long int _movingObjectId;
       int _procHistoryId;
       double _ra;
       double _decl;
       int ___zoneId_placeholder;
       double _raErr4wcs;
       double _decErr4wcs;
       double _raErr4detection;
       double _decErr4detection;
       // _xpos DECIMAL(7,2) NOT NULL; // _rowc above
       // _ypos DECIMAL(7,2) NOT NULL; // _colc above
       double _rowcErr;                // WAS : _xposErr DECIMAL(4,2) NOT NULL;
       double _colcErr;                // WAS : _yposErr DECIMAL(4,2) NOT NULL;
       double _cx;
       double _cy;
       double _cz;
       double _taiMidPoint;
       double _taiRange;
       float _fwhmA;
       float _fwhmB;
       float _fwhmTheta;
       double _flux;
       double _fluxErr;
       float _psfMag;
       float _psfMagErr;
       float _apMag;
       float _apMagErr;
       float _apDia;
       float _petroMag;
       float _petroMagErr;
       float _snr;
       float _chi2;
       float _sky;
       float _skyErr;
       float _moment0;
       float _moment1_x;
       float _moment1_y;
       float _moment2_xx;
       float _moment2_xy;
       float _moment2_yy;
       float _moment3_xxx;
       float _moment3_xxy;
       float _moment3_xyy;
       float _moment3_yyy;
       float _moment4_xxxx;
       float _moment4_xxxy;
       float _moment4_xxyy;
       float _moment4_xyyy;
       float _moment4_yyyy;
       short int _flag4association;
       short int _flag4detection;
       short int _flag4wcs;
    };
    
}   // namespace fw
}   // namespace lsst

#include <lsst/fw/Source.cc>

#endif // !defined(LSST_FW_Source_H)
