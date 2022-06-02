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
 * @file ImageStatistics.h
 *
 * @brief Image Subtraction helper functions
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_IMAGESTATISTICS_H
#define LSST_IP_DIFFIM_IMAGESTATISTICS_H

#include <cmath>
#include <limits>
#include <memory>

#include "lsst/afw/image.h"
#include "lsst/log/Log.h"
#include "lsst/daf/base/PropertySet.h"

namespace lsst {
namespace ip {
namespace diffim {

    /**
     * @brief Class to calculate difference image statistics
     *
     * @note Find mean and unbiased variance of pixel residuals in units of
     * sqrt(variance)
     *
     * @ingroup ip_diffim
     */
    template <typename PixelT>
    class ImageStatistics {
    public:
        typedef std::shared_ptr<ImageStatistics> Ptr;
        typedef typename lsst::afw::image::MaskedImage<PixelT>::x_iterator x_iterator;

        ImageStatistics(lsst::daf::base::PropertySet const& ps) :
        _xsum(0.), _x2sum(0.), _npix(0), _bpMask(0) {

            std::vector<std::string> detBadMaskPlanes = ps.getArray<std::string>("badMaskPlanes");
            for (std::vector<std::string>::iterator mi = detBadMaskPlanes.begin();
                 mi != detBadMaskPlanes.end(); ++mi){

                try {
                    _bpMask |= lsst::afw::image::Mask<lsst::afw::image::MaskPixel>::getPlaneBitMask(*mi);
                } catch (pexExcept::Exception& e) {
                    LOGL_DEBUG("TRACE4.ip.diffim.ImageStatistics",
                               "Cannot update bad bit mask with %s", (*mi).c_str());
                    LOGL_DEBUG("TRACE5.ip.diffim.ImageStatistics",
                               e.what());
                }
            }
        } ;
        virtual ~ImageStatistics() {} ;

        // Clear the accumulators
        void reset() { _xsum = _x2sum = 0.; _npix = 0;}

        // Work your magic
        void apply(lsst::afw::image::MaskedImage<PixelT> const& image) {
            apply(image, -1);
        }

        void apply(lsst::afw::image::MaskedImage<PixelT> const& image, int core) {
            reset();
            int y0, y1, x0, x1;
            if (core == -1) {
                y0 = 0;
                y1 = image.getHeight();
                x0 = 0;
                x1 = image.getWidth();
            }
            else {
                y0 = std::max(0, image.getHeight()/2 - core);
                y1 = std::min(image.getHeight(), image.getHeight()/2 + core + 1);
                x0 = std::max(0, image.getWidth()/2 - core);
                x1 = std::min(image.getWidth(), image.getWidth()/2 + core + 1);
            }

            for (int y = y0; y != y1; ++y) {
                for (x_iterator ptr = image.x_at(x0, y), end = image.x_at(x1, y);
                     ptr != end; ++ptr) {
                    if (!((*ptr).mask() & _bpMask)) {
                        double const ivar = 1. / (double) (*ptr).variance();
                        if (std::isfinite(ivar)) {
                            _xsum  += (double) (*ptr).image() * sqrt(ivar);
                            _x2sum += (double) (*ptr).image() * (double) (*ptr).image() * ivar;
                            _npix  += 1;
                        }
                    }
                }
            }
            if ((!std::isfinite(_xsum)) || (!std::isfinite(_x2sum))) {
                throw LSST_EXCEPT(pexExcept::Exception,
                                  "Nan/Inf in ImageStatistics.apply");
            }
        }

        void setBpMask(lsst::afw::image::MaskPixel bpMask) {_bpMask = bpMask;}
        lsst::afw::image::MaskPixel getBpMask() {return _bpMask;}

        // Mean of distribution
        double getMean() const {
            return (_npix > 0) ? _xsum/_npix : std::numeric_limits<double>::quiet_NaN();
        }
        // Variance of distribution
        double getVariance() const {
            return (_npix > 1) ? (_x2sum/_npix - _xsum/_npix * _xsum/_npix) * _npix/(_npix-1.) :
                std::numeric_limits<double>::quiet_NaN();
        }
        // RMS
        double getRms() const {
            return sqrt(getVariance());
        }
        // Return the number of good pixels
        int getNpix() const { return _npix; }

        // Return Sdqa rating
        bool evaluateQuality(lsst::daf::base::PropertySet const& ps) {
            if ( fabs(getMean())     > ps.getAsDouble("maximumFootprintResidualMean") ) return false;
            if ( getRms()            > ps.getAsDouble("maximumFootprintResidualStd")  ) return false;
            return true;
        }

    private:
        double _xsum;
        double _x2sum;
        int    _npix;
        lsst::afw::image::MaskPixel _bpMask;
    };


}}} // end of namespace lsst::ip::diffim


#endif
