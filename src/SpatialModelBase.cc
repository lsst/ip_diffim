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

#include <lsst/ip/diffim/SpatialModelBase.h>

namespace lsst {
namespace ip {
namespace diffim {
    
SpatialModelBase::SpatialModelBase() :
    _isGood(false),
    _isBuilt(false)
{;}

}}} // end of namespace lsst::ip::diffim
