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

#ifndef LSST_FW_wcsMatch_H
#define LSST_FW_wcsMatch_H

#include <string>

#include <boost/shared_ptr.hpp>

#include <lsst/afw/Exposure.h>
#include <lsst/afw/WCS.h>

namespace lsst {
namespace ip {
namespace diffim {
       
    typedef boost::uint16_t maskPixelType;

    template<typename ImageT, typename MaskT> 
    int wcsMatch(
        lsst::afw::Exposure<ImageT, MaskT> &remapExposure,
        lsst::afw::Exposure<ImageT, MaskT> const &origExposure,
        std::string const kernelType, 
        int kernelCols, 
        int kernelRows
        );

    template<typename ImageT, typename MaskT> 
    lsst::afw::Exposure<ImageT, MaskT> wcsMatch(
        int &numEdgePixels,
        lsst::afw::WCS const &remapWcs,
        int remapCols, 
        int remapRows, 
        lsst::afw::Exposure<ImageT, MaskT> const &origExposure,
        std::string const kernelType, 
        int kernelCols, 
        int kernelRows
        );
       
}}} // namespace diffim   ip   lsst

#endif // !defined(LSST_FW_wcsMatch_H)
