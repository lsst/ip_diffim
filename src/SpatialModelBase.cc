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
#include <lsst/ip/diffim/SpatialModelBase.h>

namespace lsst {
namespace ip {
namespace diffim {
    
/** Build the model
 *
 * @note Should be overridden by all derived classes
 */
template <typename ImageT, typename MaskT>
bool SpatialModelBase<ImageT, MaskT>::buildModel() {
    return false;
}

/** Return quality of the model
 *
 * @note Not yet implemented 
 *
 * @note Currently each derived class overrides this
 */
template <typename ImageT, typename MaskT>
double SpatialModelBase<ImageT, MaskT>::returnSdqaRating() {
    return 0.0;
}

// Explicit instantiations
template class SpatialModelBase<float, lsst::afw::image::maskPixelType>;
template class SpatialModelBase<double, lsst::afw::image::maskPixelType>;

}}} // end of namespace lsst::ip::diffim
