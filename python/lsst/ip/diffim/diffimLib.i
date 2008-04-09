// -*- lsst-c++ -*-
%define diffimLib_DOCSTRING
"
Python bindings for lsst::ip::diffim code
"
%enddef

%feature("autodoc", "1");
%module(docstring=diffimLib_DOCSTRING) diffimLib

// Everything we will need in the _wrap.cc file
%{
#include <lsst/ip/diffim/ImageSubtract.h>
#include <lsst/ip/diffim/wcsMatch.h>
%}

%init %{
%}

namespace boost {
    class bad_any_cast; // remove warning: Nothing known about 'boost::bad_any_cast'
}

// Everything whose bindings we will have to know about
%include "lsst/p_lsstSwig.i"    // this needs to go first otherwise i do not know about e.g. boost
%include "lsst/afw/image/lsstImageTypes.i"  // vw and Image/Mask types and typedefs
%include "lsst/detection/detectionLib.i"    // need for FootprintContainerT

// handle C++ arguments that should be outputs in python
%apply int& OUTPUT { int& };
%apply float& OUTPUT { float& };
%apply double& OUTPUT { double& };

%pythoncode %{
def version(HeadURL = r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/ip/diffim/tickets/7/python/lsst/ip/diffim/diffimLib.i $"):
    """Return a version given a HeadURL string; default: ip_diffim's version"""
    return guessSvnVersion(HeadURL)

%}

// Here you give the names of the functions you want to Swig

/*******************/
/* ImageSubtract.h */

%include "lsst/ip/diffim/ImageSubtract.h"
%template(DiffImContainerD)                lsst::ip::diffim::DiffImContainer<float, lsst::afw::image::maskPixelType>;
%template(vectorDiffImContainerD)          std::vector<lsst::ip::diffim::DiffImContainer<float, lsst::afw::image::maskPixelType> >;

%template(computeDiffImStats)
    lsst::ip::diffim::computeDiffImStats<float, lsst::afw::image::maskPixelType>;
%template(computeDiffImStats)
    lsst::ip::diffim::computeDiffImStats<double, lsst::afw::image::maskPixelType>;

%template(computePsfMatchingKernelForPostageStamp)
    lsst::ip::diffim::computePsfMatchingKernelForPostageStamp<float, lsst::afw::image::maskPixelType>;
%template(computePsfMatchingKernelForPostageStamp)
    lsst::ip::diffim::computePsfMatchingKernelForPostageStamp<double, lsst::afw::image::maskPixelType>;

%template(getCollectionOfFootprintsForPsfMatching)
    lsst::ip::diffim::getCollectionOfFootprintsForPsfMatching<float, lsst::afw::image::maskPixelType>;
%template(getCollectionOfFootprintsForPsfMatching)
    lsst::ip::diffim::getCollectionOfFootprintsForPsfMatching<double, lsst::afw::image::maskPixelType>;

%template(computePcaKernelBasis)            lsst::ip::diffim::computePcaKernelBasis<float, lsst::afw::image::maskPixelType>;
%template(computePcaKernelBasis)            lsst::ip::diffim::computePcaKernelBasis<double, lsst::afw::image::maskPixelType>;

%template(computeSpatiallyVaryingPsfMatchingKernel)
    lsst::ip::diffim::computeSpatiallyVaryingPsfMatchingKernel<float, lsst::afw::image::maskPixelType>;
%template(computeSpatiallyVaryingPsfMatchingKernel)
    lsst::ip::diffim::computeSpatiallyVaryingPsfMatchingKernel<double, lsst::afw::image::maskPixelType>;

%template(maskOk)                           lsst::ip::diffim::maskOk<lsst::afw::image::maskPixelType>;

%template(calculateMaskedImageResiduals)
    lsst::ip::diffim::calculateMaskedImageResiduals<float, lsst::afw::image::maskPixelType>;
%template(calculateMaskedImageResiduals)
    lsst::ip::diffim::calculateMaskedImageResiduals<double, lsst::afw::image::maskPixelType>;

%template(calculateImageResiduals)          lsst::ip::diffim::calculateImageResiduals<float>;
%template(calculateImageResiduals)          lsst::ip::diffim::calculateImageResiduals<double>;
%template(addFunction)                      lsst::ip::diffim::addFunction<double, double>;
%template(addFunction)                      lsst::ip::diffim::addFunction<float, double>;
%template(addFunction)                      lsst::ip::diffim::addFunction<double, float>;
%template(addFunction)                      lsst::ip::diffim::addFunction<float, float>;

/* ImageSubtract.h */
/*******************/
/* PCA.h */

%include "lsst/ip/diffim/PCA.h"
%template(computePca)                    lsst::ip::diffim::computePca<vw::math::Matrix<float>, vw::math::Vector<float> >;
%template(computePca)                    lsst::ip::diffim::computePca<vw::math::Matrix<double>, vw::math::Vector<double> >;
%template(decomposeMatrixUsingBasis)     lsst::ip::diffim::decomposeMatrixUsingBasis<vw::math::Matrix<float> >;
%template(decomposeMatrixUsingBasis)     lsst::ip::diffim::decomposeMatrixUsingBasis<vw::math::Matrix<double> >;
%template(approximateMatrixUsingBasis)   lsst::ip::diffim::approximateMatrixUsingBasis<vw::math::Matrix<float> >;
%template(approximateMatrixUsingBasis)   lsst::ip::diffim::approximateMatrixUsingBasis<vw::math::Matrix<double> >;

/* PCA.h */

/*******************/
/* WCSMatch.h */
%include "lsst/ip/diffim/wcsMatch.h"
%template(wcsMatch)    lsst::ip::diffim::wcsMatch<boost::uint16_t, lsst::afw::image::maskPixelType>;
%template(wcsMatch)    lsst::ip::diffim::wcsMatch<float, lsst::afw::image::maskPixelType>;
%template(wcsMatch)    lsst::ip::diffim::wcsMatch<double, lsst::afw::image::maskPixelType>;

/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
