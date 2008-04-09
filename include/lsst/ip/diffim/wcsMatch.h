// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file 
  *
  * \ingroup diffim
  *
  * \brief Implementation of the templated utility function, wcsMatch, for
  * Astrometric Image Remapping for the LSST.
  *
  * \author Nicole M. Silvestri, University of Washington
  *
  * Contact: nms@astro.washington.edu 
  *
  * \version
  *
  * LSST Legalese here...
  */

#ifndef LSST_IP_DIFFIM_WCSMATCH_H
#define LSST_IP_DIFFIM_WCSMATCH_H

#include <string>

#include <boost/shared_ptr.hpp>

#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/image/Wcs.h>

namespace lsst {
namespace ip {
namespace diffim {
       
    typedef boost::uint16_t maskPixelType;

    template<typename ImageT, typename MaskT> 
    int wcsMatch(
        lsst::afw::image::Exposure<ImageT, MaskT> &remapExposure,
        lsst::afw::image::Exposure<ImageT, MaskT> const &origExposure,
        std::string const kernelType, 
        int kernelCols, 
        int kernelRows
        );

    template<typename ImageT, typename MaskT> 
    lsst::afw::image::Exposure<ImageT, MaskT> wcsMatch(
        int &numEdgePixels,
        lsst::afw::image::Wcs const &remapWcs,
        int remapCols, 
        int remapRows, 
        lsst::afw::image::Exposure<ImageT, MaskT> const &origExposure,
        std::string const kernelType, 
        int kernelCols, 
        int kernelRows
        );
       
}}} // namespace diffim   ip   lsst

#endif // !defined(LSST_IP_DIFFIM_WCSMATCH_H)
