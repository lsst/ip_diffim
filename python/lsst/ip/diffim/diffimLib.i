// -*- lsst-c++ -*-
%define diffimLib_DOCSTRING
"
Python bindings for lsst::ip::diffim code
"
%enddef

%feature("autodoc", "1");
%module(docstring=diffimLib_DOCSTRING, package="lsst.ip.diffim") diffimLib

/* Everything needed in the _wrap.cc file */
%{
#include <lsst/ip/diffim/ImageSubtract.h>
#include <lsst/ip/diffim/wcsMatch.h>
#include <lsst/ip/diffim/Pca.h>
%}

/* Remove warnings */
/* Nothing known about 'boost::bad_any_cast' */
namespace boost {
    class bad_any_cast; 
}

%init %{
%}

/* Everything whose bindings are needed in the wrapper */
%include "lsst/p_lsstSwig.i"                // this needs to go first otherwise i do not know about e.g. boost
%include "lsst/afw/image/lsstImageTypes.i"  // vw and Image/Mask types and typedefs
%include "lsst/detection/detectionLib.i"    // need for FootprintContainerT

/* Everything whose bindings we need to know about but not wrap, e.g. MaskedImageF */
%import "lsst/afw/image/imageLib.i" 

/* handle C++ arguments that should be outputs in python */
%apply int& OUTPUT { int& };
%apply float& OUTPUT { float& };
%apply double& OUTPUT { double& };

%pythoncode %{
def version(HeadURL = r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/ip/diffim/tickets/7/python/lsst/ip/diffim/diffimLib.i $"):
    """Return a version given a HeadURL string; default: ip_diffim's version"""
    return guessSvnVersion(HeadURL)

%}

/* Here you give the names of the functions you want to wrap */
/* "include" your header files here */

/*******************/
/* ImageSubtract.h */

%include "lsst/ip/diffim/ImageSubtract.h"

/* classes */
%template(DifferenceImageStatisticsF)
    lsst::ip::diffim::DifferenceImageStatistics<float, lsst::afw::image::maskPixelType>;

%template(DifferenceImageFootprintInformationF)  
    lsst::ip::diffim::DifferenceImageFootprintInformation<float, lsst::afw::image::maskPixelType>;
%boost_shared_ptr(DifiPtrF, 
                  lsst::ip::diffim::DifferenceImageFootprintInformation<float, lsst::afw::image::maskPixelType>);
%template(DifiPtrListF)
    std::vector<lsst::ip::diffim::DifferenceImageFootprintInformation<float, lsst::afw::image::maskPixelType>::Ptr>;

/* subroutines */
%template(getGoodFootprints)
    lsst::ip::diffim::getGoodFootprints<float, lsst::afw::image::maskPixelType>;
%template(getGoodFootprints)
    lsst::ip::diffim::getGoodFootprints<double, lsst::afw::image::maskPixelType>;

%template(convolveAndSubtract)
    lsst::ip::diffim::convolveAndSubtract<float, lsst::afw::image::maskPixelType>;
%template(convolveAndSubtract)
    lsst::ip::diffim::convolveAndSubtract<double, lsst::afw::image::maskPixelType>;

%template(computePsfMatchingKernelForFootprint)
    lsst::ip::diffim::computePsfMatchingKernelForFootprint<float, lsst::afw::image::maskPixelType>;
%template(computePsfMatchingKernelForFootprint)
    lsst::ip::diffim::computePsfMatchingKernelForFootprint<double, lsst::afw::image::maskPixelType>;

%template(computePsfMatchingKernelForFootprint2)
    lsst::ip::diffim::computePsfMatchingKernelForFootprint2<float, lsst::afw::image::maskPixelType>;
%template(computePsfMatchingKernelForFootprint2)
    lsst::ip::diffim::computePsfMatchingKernelForFootprint2<double, lsst::afw::image::maskPixelType>;

%template(getCollectionOfFootprintsForPsfMatching)
    lsst::ip::diffim::getCollectionOfFootprintsForPsfMatching<float, lsst::afw::image::maskPixelType>;
%template(getCollectionOfFootprintsForPsfMatching)
    lsst::ip::diffim::getCollectionOfFootprintsForPsfMatching<double, lsst::afw::image::maskPixelType>;

%template(maskOk)                           lsst::ip::diffim::maskOk<lsst::afw::image::maskPixelType>;

%template(calculateMaskedImageStatistics)
    lsst::ip::diffim::calculateMaskedImageStatistics<float, lsst::afw::image::maskPixelType>;
%template(calculateMaskedStatistics)
    lsst::ip::diffim::calculateMaskedImageStatistics<double, lsst::afw::image::maskPixelType>;

%template(calculateImageStatistics)         lsst::ip::diffim::calculateImageStatistics<float>;
%template(calculateImageStatistics)         lsst::ip::diffim::calculateImageStatistics<double>;
%template(addFunctionToImage)               lsst::ip::diffim::addFunctionToImage<double, double>;
%template(addFunctionToImage)               lsst::ip::diffim::addFunctionToImage<float, double>;
%template(addFunctionToImage)               lsst::ip::diffim::addFunctionToImage<double, float>;
%template(addFunctionToImage)               lsst::ip::diffim::addFunctionToImage<float, float>;

/* ImageSubtract.h */
/*******************/
/* Pca.h */

%include "lsst/ip/diffim/Pca.h"
%template(computePca)                    lsst::ip::diffim::computePca<vw::math::Matrix<float>, vw::math::Vector<float> >;
%template(computePca)                    lsst::ip::diffim::computePca<vw::math::Matrix<double>, vw::math::Vector<double> >;
%template(decomposeMatrixUsingBasis)     lsst::ip::diffim::decomposeMatrixUsingBasis<vw::math::Matrix<float> >;
%template(decomposeMatrixUsingBasis)     lsst::ip::diffim::decomposeMatrixUsingBasis<vw::math::Matrix<double> >;
%template(approximateMatrixUsingBasis)   lsst::ip::diffim::approximateMatrixUsingBasis<vw::math::Matrix<float> >;
%template(approximateMatrixUsingBasis)   lsst::ip::diffim::approximateMatrixUsingBasis<vw::math::Matrix<double> >;

/* Pca.h */

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
