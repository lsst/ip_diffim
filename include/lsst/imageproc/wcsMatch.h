// -*- LSST-C++ -*- // fixed format comment for emacs
/**
  * \file 
  *
  * \ingroup imageproc
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

#include <lsst/fw/Exposure.h>
#include <lsst/fw/WCS.h>

namespace lsst {
namespace imageproc {
       
    typedef boost::uint16_t maskPixelType;

    template<typename ImageT, typename MaskT> 
    int wcsMatch(
        lsst::fw::Exposure<ImageT, MaskT> &remapExposure,
        lsst::fw::Exposure<ImageT, MaskT> const &origExposure,
        std::string const kernelType, 
        int kernelCols, 
        int kernelRows
        );

    template<typename ImageT, typename MaskT> 
    lsst::fw::Exposure<ImageT, MaskT> wcsMatch(
        int &numEdgePixels,
        lsst::fw::WCS const &remapWcs,
        int remapCols, 
        int remapRows, 
        lsst::fw::Exposure<ImageT, MaskT> const &origExposure,
        std::string const kernelType, 
        int kernelCols, 
        int kernelRows
        );
       
}} // namespace fw::lsst

#ifndef SWIG // don't bother SWIG with .cc files
#include <lsst/imageproc/wcsMatch.cc>
#endif

#endif // !defined(LSST_FW_wcsMatch_H)
