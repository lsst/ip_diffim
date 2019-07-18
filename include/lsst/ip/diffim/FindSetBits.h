// -*- lsst-c++ -*-

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

/**
 * @file FindSetBits.h
 *
 * @brief Image Subtraction helper functions
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_FINDSETBITS_H
#define LSST_IP_DIFFIM_FINDSETBITS_H

#include "lsst/afw/image.h"

namespace lsst {
namespace ip {
namespace diffim {


    /**
     * @brief Class to accumulate Mask bits
     *
     * @note Search through a Mask for any set bits.
     *
     * @ingroup ip_diffim
     */
    template <typename MaskT>
    class FindSetBits {
    public:
        typedef typename MaskT::x_iterator x_iterator;

        FindSetBits() :
            _bits(0) {;}
        virtual ~FindSetBits() {} ;

        // Clear the accumulators
        void reset() { _bits = 0;}

        // Return the bits set
        typename MaskT::Pixel getBits() const { return _bits; }

        // Work your magic
        void apply(MaskT const& mask) {
            reset();
            for (int y = 0; y != mask.getHeight(); ++y) {
                for (x_iterator ptr = mask.row_begin(y), end = mask.row_end(y); ptr != end; ++ptr) {
                    _bits |= (*ptr);
                }
            }
        }

    private:
        typename MaskT::Pixel _bits;
    };

}}} // end of namespace lsst::ip::diffim


#endif
