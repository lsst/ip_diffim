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
    
template <typename ImageT>
SpatialModelBase<ImageT>::SpatialModelBase() :
    _id(-1),
    _isBuilt(false),
    _isGood(false),
    _colc(0.0),
    _rowc(0.0)
{;}

// Explicit instantiations
template class SpatialModelBase<float>;
template class SpatialModelBase<double>;

template class cmpSpatialModels<float>;
template class cmpSpatialModels<double>;

}}} // end of namespace lsst::ip::diffim
