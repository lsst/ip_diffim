// -*- lsst-c++ -*-
%define diffimLib_DOCSTRING
"
Python bindings for lsst::ip::diffim code
"
%enddef

%feature("autodoc", "1");
%module(docstring=diffimLib_DOCSTRING, package="lsst.ip.diffim") diffimLib

// Reference for this file is at http://dev.lsstcorp.org/trac/wiki/SwigFAQ 
// Nice practical example is at http://dev.lsstcorp.org/trac/browser/DMS/afw/trunk/python/lsst/afw/image/imageLib.i 

// Suppress swig complaints
#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
#pragma SWIG nowarn=362                 // operator=  ignored

// Remove warnings
// Nothing known about 'boost::bad_any_cast'
namespace boost {
    class bad_any_cast; 
    //typedef unsigned short boost::uint16_t;
}

/* handle C++ arguments that should be outputs in python */
%apply int& OUTPUT { int& };
%apply float& OUTPUT { float& };
%apply double& OUTPUT { double& };

%{
#include <boost/shared_ptr.hpp>

#include <lsst/afw/image.h>
#include <lsst/afw/math.h>
#include <lsst/afw/detection.h>
#include <lsst/afw/math/ConvolveImage.h>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>

#include <vw/Math/Functions.h> 
#include <vw/Math/Vector.h> 
#include <vw/Math/Matrix.h> 
#include <vw/Math/LinearAlgebra.h> 

#include <lsst/pex/policy/Policy.h>
%}

/******************************************************************************/

%include "lsst/p_lsstSwig.i"
%include "lsst/daf/base/persistenceMacros.i"
%import  "lsst/afw/image/imageLib.i" 
%import  "lsst/afw/detection/detectionLib.i"
%import  "lsst/afw/math/kernelLib.i"

%lsst_exceptions();

%pythoncode %{
import lsst.utils

def version(HeadURL = r"$HeadURL$"):
    """Return a version given a HeadURL string. If a different version is setup, return that too"""

    version_svn = lsst.utils.guessSvnVersion(HeadURL)

    try:
        import eups
    except ImportError:
        return version_svn
    else:
        try:
            version_eups = eups.setup("ip_diffim")
        except AttributeError:
            return version_svn

    if version_eups == version_svn:
        return version_svn
    else:
        return "%s (setup: %s)" % (version_svn, version_eups)

%}

%{
#include "lsst/ip/diffim/ImageSubtract.h"
%}

%apply boost::shared_ptr<lsst::afw::math::Kernel>& OUTPUT { boost::shared_ptr<lsst::afw::math::Kernel>& };
%include "lsst/ip/diffim/ImageSubtract.h"

/*
%template(FindSetBitsU)
%    lsst::ip::diffim::FindSetBits<lsst::afw::image::Mask>;

%template(FindCountsF)
    lsst::ip::diffim::FindCounts<lsst::afw::image::MaskedImage<float> >;
%template(FindCountsD)
    lsst::ip::diffim::FindCounts<lsst::afw::image::MaskedImage<double> >;

%template(ImageStatisticsF)
    lsst::ip::diffim::ImageStatistics<lsst::afw::image::MaskedImage<float> >;
%template(ImageStatisticsD)
    lsst::ip::diffim::ImageStatistics<lsst::afw::image::MaskedImage<double> >;
*/

%template(convolveAndSubtract)
    lsst::ip::diffim::convolveAndSubtract<float>;
%template(convolveAndSubtract)
    lsst::ip::diffim::convolveAndSubtract<double>;

%template(convolveAndSubtract)
    lsst::ip::diffim::convolveAndSubtract<float, double>;
%template(convolveAndSubtract)
    lsst::ip::diffim::convolveAndSubtract<double, double>;

%template(computePsfMatchingKernelForFootprint)
    lsst::ip::diffim::computePsfMatchingKernelForFootprint<float, float>;
%template(computePsfMatchingKernelForFootprint)
    lsst::ip::diffim::computePsfMatchingKernelForFootprint<double, float>;

%template(computePsfMatchingKernelForFootprintEigen)
    lsst::ip::diffim::computePsfMatchingKernelForFootprintEigen<float, float>;
%template(computePsfMatchingKernelForFootprintEigen)
    lsst::ip::diffim::computePsfMatchingKernelForFootprintEigen<double, float>;

%template(computePsfMatchingKernelForFootprintVW)
    lsst::ip::diffim::computePsfMatchingKernelForFootprintVW<float, float>;
%template(computePsfMatchingKernelForFootprintVW)
    lsst::ip::diffim::computePsfMatchingKernelForFootprintVW<double, float>;
%clear boost::shared_ptr<lsst::afw::math::Kernel>& OUTPUT;

%template(getCollectionOfFootprintsForPsfMatching)
    lsst::ip::diffim::getCollectionOfFootprintsForPsfMatching<float>;
%template(getCollectionOfFootprintsForPsfMatching)
    lsst::ip::diffim::getCollectionOfFootprintsForPsfMatching<double>;

%template(addFunctionToImage)               lsst::ip::diffim::addFunctionToImage<double, double>;
%template(addFunctionToImage)               lsst::ip::diffim::addFunctionToImage<float, double>;
%template(addFunctionToImage)               lsst::ip::diffim::addFunctionToImage<double, float>;
%template(addFunctionToImage)               lsst::ip::diffim::addFunctionToImage<float, float>;

/******************************************************************************/

%{
#include "lsst/ip/diffim/SpatialModelBase.h"
%}

SWIG_SHARED_PTR(SpatialModelBaseF, lsst::ip::diffim::SpatialModelBase<float>);
SWIG_SHARED_PTR(SpatialModelBaseD, lsst::ip::diffim::SpatialModelBase<double>);

%include "lsst/ip/diffim/SpatialModelBase.h"

%template(SpatialModelBaseF) lsst::ip::diffim::SpatialModelBase<float>;
%template(SpatialModelBaseD) lsst::ip::diffim::SpatialModelBase<double>;

%template(VectorSpatialModelBaseF) std::vector<lsst::ip::diffim::SpatialModelBase<float>::Ptr >;
%template(VectorSpatialModelBaseD) std::vector<lsst::ip::diffim::SpatialModelBase<double>::Ptr >;

/******************************************************************************/

%{
#include "lsst/ip/diffim/SpatialModelKernel.h"
%}

SWIG_SHARED_PTR(SpatialModelKernelF, lsst::ip::diffim::SpatialModelKernel<float>);
SWIG_SHARED_PTR(SpatialModelKernelD, lsst::ip::diffim::SpatialModelKernel<double>);

%include "lsst/ip/diffim/SpatialModelKernel.h"

%template(SpatialModelKernelF) lsst::ip::diffim::SpatialModelKernel<float>;
%template(SpatialModelKernelD) lsst::ip::diffim::SpatialModelKernel<double>;

%template(VectorSpatialModelKernelF) std::vector<lsst::ip::diffim::SpatialModelKernel<float>::Ptr >;
%template(VectorSpatialModelKernelD) std::vector<lsst::ip::diffim::SpatialModelKernel<double>::Ptr >;


/******************************************************************************/

%{
#include "lsst/ip/diffim/SpatialModelCell.h"
%}

SWIG_SHARED_PTR(SpatialModelCellF, lsst::ip::diffim::SpatialModelCell<float>);
SWIG_SHARED_PTR(SpatialModelCellD, lsst::ip::diffim::SpatialModelCell<double>);

%include "lsst/ip/diffim/SpatialModelCell.h"

%template(SpatialModelCellF) lsst::ip::diffim::SpatialModelCell<float>;
%template(SpatialModelCellD) lsst::ip::diffim::SpatialModelCell<double>;

%template(VectorSpatialModelCellF) std::vector<lsst::ip::diffim::SpatialModelCell<float>::Ptr >;
%template(VectorSpatialModelCellD) std::vector<lsst::ip::diffim::SpatialModelCell<double>::Ptr >;

/******************************************************************************/

/*
%{
#include "lsst/ip/diffim/Pca.h"
%}

%include "lsst/ip/diffim/Pca.h"

%template(computePca)                    lsst::ip::diffim::computePca<vw::math::Matrix<float>, vw::math::Vector<float> >;
%template(computePca)                    lsst::ip::diffim::computePca<vw::math::Matrix<double>, vw::math::Vector<double> >;
%template(decomposeMatrixUsingBasis)     lsst::ip::diffim::decomposeMatrixUsingBasis<vw::math::Matrix<float> >;
%template(decomposeMatrixUsingBasis)     lsst::ip::diffim::decomposeMatrixUsingBasis<vw::math::Matrix<double> >;
%template(approximateMatrixUsingBasis)   lsst::ip::diffim::approximateMatrixUsingBasis<vw::math::Matrix<float> >;
%template(approximateMatrixUsingBasis)   lsst::ip::diffim::approximateMatrixUsingBasis<vw::math::Matrix<double> >;
*/
