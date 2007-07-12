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

    class Source : private lsst::fw::LsstBase {

    public:
        explicit Source();
        virtual ~Source() {};
        explicit Source(
            double rowc,
            double colc,
            double drow,
            double dcol
            );
       inline double getColc() const;
       inline double getRowc() const;
       inline double getDcol() const;
       inline double getDrow() const;
    private:
       double _rowc;
       double _colc;
       double _drow;
       double _dcol;
    };
    
}   // namespace fw
}   // namespace lsst

#include <lsst/fw/Source.cc>

#endif // !defined(LSST_FW_Source_H)
