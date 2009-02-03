// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Implementation of SpatialModelBase class
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/Image.h>

#include <lsst/ip/diffim/SpatialModelBase.h>

namespace lsst {
namespace ip {
namespace diffim {
    
template <typename ImageT, typename MaskT>
SpatialModelBase<ImageT, MaskT>::SpatialModelBase() :
    _id(-1),
    _isBuilt(false),
    _isGood(false),
    _colcNorm(0.0),
    _rowcNorm(0.0)
{;}

// Explicit instantiations
template class SpatialModelBase<float, lsst::afw::image::MaskPixel>;
template class SpatialModelBase<double, lsst::afw::image::MaskPixel>;

}}} // end of namespace lsst::ip::diffim
