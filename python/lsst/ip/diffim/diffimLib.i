// -*- lsst-c++ -*-
%define diffimLib_DOCSTRING
"
Python bindings for lsst::ip::diffim code
"
%enddef

%feature("autodoc", "1");
%module(docstring=diffimLib_DOCSTRING, package="lsst.ip.diffim") diffimLib

// Avoid Swig bug exposed in Ticket 640
// http://dev.lsstcorp.org/trac/ticket/640
%ignore lsst::ip::diffim::FindSetBits::operator();
%ignore lsst::ip::diffim::ImageStatistics::operator();

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

#include <lsst/afw.h>

#include <lsst/ip/diffim/ImageSubtract.h>
#include <lsst/ip/diffim/BasisSets.h>
#include <lsst/ip/diffim/PsfMatchingFunctor.h>
#include <lsst/ip/diffim/SpatialModelKernel.h>

#include <lsst/pex/policy/Policy.h>
%}

/******************************************************************************/


%include "lsst/p_lsstSwig.i"
%include "lsst/daf/base/persistenceMacros.i"
%import  "lsst/afw/image/imageLib.i"
%import  "lsst/afw/math/mathLib.i"
%import  "lsst/afw/detection/detectionLib.i"

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

/******************************************************************************/

%{
#include "lsst/ip/diffim/PsfMatchingFunctor.h"
%}

SWIG_SHARED_PTR(PsfMatchingFunctorF, lsst::ip::diffim::PsfMatchingFunctor<float>);
SWIG_SHARED_PTR(PsfMatchingFunctorD, lsst::ip::diffim::PsfMatchingFunctor<double>);

%include "lsst/ip/diffim/PsfMatchingFunctor.h"

%template(PsfMatchingFunctorF)
    lsst::ip::diffim::PsfMatchingFunctor<float>;
%template(PsfMatchingFunctorD)
    lsst::ip::diffim::PsfMatchingFunctor<double>;
%template(makePsfMatchingFunctorF) lsst::ip::diffim::makePsfMatchingFunctor<float>;
%template(makePsfMatchingFunctorD) lsst::ip::diffim::makePsfMatchingFunctor<double>;

%template(pair_Kernel_double)   std::pair<lsst::afw::math::Kernel::Ptr, double>;

/******************************************************************************/

%{
#include "lsst/ip/diffim/ImageSubtract.h"
%}

%include "lsst/ip/diffim/ImageSubtract.h"

%template(FindSetBitsU)
    lsst::ip::diffim::FindSetBits<lsst::afw::image::Mask<lsst::afw::image::MaskPixel> >;

%template(ImageStatisticsI)
    lsst::ip::diffim::ImageStatistics<int>;
%template(ImageStatisticsF)
    lsst::ip::diffim::ImageStatistics<float>;
%template(ImageStatisticsD)
    lsst::ip::diffim::ImageStatistics<double>;

%define %convolveAndSubtract(PIXEL_T)
   %template(convolveAndSubtract)
       lsst::ip::diffim::convolveAndSubtract<PIXEL_T, double>;
   %template(convolveAndSubtract)
       lsst::ip::diffim::convolveAndSubtract<PIXEL_T, lsst::afw::math::Function2<double> const&>;
%enddef

%convolveAndSubtract(float);
#if 0
%convolveAndSubtract(double);           // image subtraction on double images??!?  Not instantiated in .cc either
#endif

%template(getCollectionOfFootprintsForPsfMatching)
    lsst::ip::diffim::getCollectionOfFootprintsForPsfMatching<float>;
%template(getCollectionOfFootprintsForPsfMatching)
    lsst::ip::diffim::getCollectionOfFootprintsForPsfMatching<double>;

%template(addSomethingToImage)   lsst::ip::diffim::addSomethingToImage<float,  lsst::afw::math::PolynomialFunction2<double> >;
%template(addSomethingToImage)   lsst::ip::diffim::addSomethingToImage<double, lsst::afw::math::PolynomialFunction2<double> >;
%template(addSomethingToImage)   lsst::ip::diffim::addSomethingToImage<float>;
%template(addSomethingToImage)   lsst::ip::diffim::addSomethingToImage<double>;

/******************************************************************************/

%{
#include "lsst/ip/diffim/SpatialModelKernel.h"
%}

%define %IMAGE(PIXTYPE)
lsst::afw::image::Image<PIXTYPE>
%enddef

%define %KernelCandidatePtr(NAME, TYPE)
SWIG_SHARED_PTR_DERIVED(KernelCandidate##NAME,
                        lsst::afw::math::SpatialCellImageCandidate<%IMAGE(lsst::afw::math::Kernel::Pixel)>,
                        lsst::ip::diffim::KernelCandidate<TYPE>);
%enddef

%define %KernelCandidate(NAME, TYPE)
%template(KernelCandidate##NAME) lsst::ip::diffim::KernelCandidate<TYPE>;
%template(makeKernelCandidate) lsst::ip::diffim::makeKernelCandidate<TYPE>;
%inline %{
    lsst::ip::diffim::KernelCandidate<TYPE> *
        cast_KernelCandidate##NAME(lsst::afw::math::SpatialCellCandidate * candidate) {
        return dynamic_cast<lsst::ip::diffim::KernelCandidate<TYPE> *>(candidate);
    }
%}
%enddef

%KernelCandidatePtr(F, float);

%include "lsst/ip/diffim/SpatialModelKernel.h"

%KernelCandidate(F, float);
%template(fitSpatialKernelFromCandidates) lsst::ip::diffim::fitSpatialKernelFromCandidates<float>;

%template(pair_Kernel_Function) std::pair<lsst::afw::math::LinearCombinationKernel::Ptr, lsst::afw::math::Kernel::SpatialFunctionPtr>;

/******************************************************************************/

//%template(eigenMatrix)          boost::shared_ptr<Eigen::MatrixXd>;

