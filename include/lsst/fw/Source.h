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

    class Source : private lsst::fw::LsstBase
    {

    public:
        explicit Source();
        virtual ~Source() {};
        explicit Source(
            unsigned rowc,
            unsigned colc,
            unsigned drow,
            unsigned dcol
            );
       inline unsigned getColc() const;
       inline unsigned getRowc() const;
       inline unsigned getDcol() const;
       inline unsigned getDrow() const;
    private:
       unsigned _rowc;
       unsigned _colc;
       unsigned _drow;
       unsigned _dcol;
    };
    
}   // namespace fw
}   // namespace lsst

#endif // !defined(LSST_FW_Source_H)
