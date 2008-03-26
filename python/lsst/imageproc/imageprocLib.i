// -*- lsst-c++ -*-
%define imageprocLib_DOCSTRING
"
Python bindings for imageproc module
"
%enddef

%feature("autodoc", "1");
%module(docstring=imageprocLib_DOCSTRING) imageprocLib

// Everything we will need in the _wrap.cc file
%{
#include <lsst/imageproc/ImageSubtract.h>
#include <lsst/imageproc/wcsMatch.h>
%}

%init %{
%}

// Everything whose bindings we will have to know about
%include "lsst/mwi/p_lsstSwig.i"             // this needs to go first otherwise i don't know about e.g. boost
%include "lsst/fw/Core/lsstImageTypes.i"     // vw and Image/Mask types and typedefs
%include "lsst/detection/detectionLib.i"     // need otherwise FootprintContainerT not known about

// handle C++ arguments that should be outputs in python
%apply int& OUTPUT { int& };
%apply float& OUTPUT { float& };
%apply double& OUTPUT { double& };

%pythoncode %{
def version(HeadURL = r"$HeadURL: svn+ssh://svn.lsstcorp.org/DC2/imageproc/tickets/7/python/lsst/imageproc/imageprocLib.i $"):
    """Return a version given a HeadURL string; default: imageproc's version"""
    return guessSvnVersion(HeadURL)

%}

// Here you give the names of the functions you want to Swig

/*******************/
/* ImageSubtract.h */

%include "lsst/imageproc/ImageSubtract.h"
%template(DiffImContainerD)                lsst::imageproc::DiffImContainer<float, lsst::fw::maskPixelType>;
%template(vectorDiffImContainerD)          std::vector<lsst::imageproc::DiffImContainer<float, lsst::fw::maskPixelType> >;

%template(computeDiffImStats)
    lsst::imageproc::computeDiffImStats<float, lsst::fw::maskPixelType>;
%template(computeDiffImStats)
    lsst::imageproc::computeDiffImStats<double, lsst::fw::maskPixelType>;

%template(computePsfMatchingKernelForPostageStamp)
    lsst::imageproc::computePsfMatchingKernelForPostageStamp<float, lsst::fw::maskPixelType>;
%template(computePsfMatchingKernelForPostageStamp)
    lsst::imageproc::computePsfMatchingKernelForPostageStamp<double, lsst::fw::maskPixelType>;

%template(getCollectionOfFootprintsForPsfMatching)
    lsst::imageproc::getCollectionOfFootprintsForPsfMatching<float, lsst::fw::maskPixelType>;
%template(getCollectionOfFootprintsForPsfMatching)
    lsst::imageproc::getCollectionOfFootprintsForPsfMatching<double, lsst::fw::maskPixelType>;

%template(computePcaKernelBasis)            lsst::imageproc::computePcaKernelBasis<float, lsst::fw::maskPixelType>;
%template(computePcaKernelBasis)            lsst::imageproc::computePcaKernelBasis<double, lsst::fw::maskPixelType>;

%template(computeSpatiallyVaryingPsfMatchingKernel)
    lsst::imageproc::computeSpatiallyVaryingPsfMatchingKernel<float, lsst::fw::maskPixelType>;
%template(computeSpatiallyVaryingPsfMatchingKernel)
    lsst::imageproc::computeSpatiallyVaryingPsfMatchingKernel<double, lsst::fw::maskPixelType>;

%template(maskOk)                           lsst::imageproc::maskOk<lsst::fw::maskPixelType>;

%template(calculateMaskedImageResiduals)
    lsst::imageproc::calculateMaskedImageResiduals<float, lsst::fw::maskPixelType>;
%template(calculateMaskedImageResiduals)
    lsst::imageproc::calculateMaskedImageResiduals<double, lsst::fw::maskPixelType>;

%template(calculateImageResiduals)          lsst::imageproc::calculateImageResiduals<float>;
%template(calculateImageResiduals)          lsst::imageproc::calculateImageResiduals<double>;
%template(addFunction)                      lsst::imageproc::addFunction<double, double>;
%template(addFunction)                      lsst::imageproc::addFunction<float, double>;
%template(addFunction)                      lsst::imageproc::addFunction<double, float>;
%template(addFunction)                      lsst::imageproc::addFunction<float, float>;

/* ImageSubtract.h */
/*******************/
/* PCA.h */

%include "lsst/imageproc/PCA.h"
%template(computePca)                    lsst::imageproc::computePca<vw::math::Matrix<float>, vw::math::Vector<float> >;
%template(computePca)                    lsst::imageproc::computePca<vw::math::Matrix<double>, vw::math::Vector<double> >;
%template(decomposeMatrixUsingBasis)     lsst::imageproc::decomposeMatrixUsingBasis<vw::math::Matrix<float> >;
%template(decomposeMatrixUsingBasis)     lsst::imageproc::decomposeMatrixUsingBasis<vw::math::Matrix<double> >;
%template(approximateMatrixUsingBasis)   lsst::imageproc::approximateMatrixUsingBasis<vw::math::Matrix<float> >;
%template(approximateMatrixUsingBasis)   lsst::imageproc::approximateMatrixUsingBasis<vw::math::Matrix<double> >;

/* PCA.h */

/*******************/
/* WCSMatch.h */
%include "lsst/imageproc/wcsMatch.h"
%template(wcsMatch)    lsst::imageproc::wcsMatch<boost::uint16_t, lsst::fw::maskPixelType>;
%template(wcsMatch)    lsst::imageproc::wcsMatch<float, lsst::fw::maskPixelType>;
%template(wcsMatch)    lsst::imageproc::wcsMatch<double, lsst::fw::maskPixelType>;

/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
