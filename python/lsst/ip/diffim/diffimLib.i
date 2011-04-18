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
// Nice practical example is at 
//     http://dev.lsstcorp.org/trac/browser/DMS/afw/trunk/python/lsst/afw/image/imageLib.i 

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
#include "boost/shared_ptr.hpp"

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"

#include "lsst/pex/policy/Policy.h"
%}

/******************************************************************************/

/* Eigen / numpy.  Info comes from $AFW_DIR/include/lsst/afw/numpyTypemaps.h */
%{
#define PY_ARRAY_UNIQUE_SYMBOL LSST_IP_DIFFIM_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "lsst/ndarray/python.h"
#include "lsst/ndarray/python/eigen.h"
%}

%init %{
    import_array();
%}

%include "lsst/p_lsstSwig.i"

%declareEigenMatrix(Eigen::MatrixXd);
%declareEigenMatrix(Eigen::VectorXd);
/* Eigen / numpy.  Info comes from $AFW_DIR/include/lsst/afw/numpyTypemaps.h */

%include "lsst/daf/base/persistenceMacros.i"

%import  "lsst/afw/image/imageLib.i"
%import  "lsst/afw/math/mathLib.i"
%import  "lsst/afw/math/objectVectors.i"
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

%template(pair_Kernel_double)   std::pair<lsst::afw::math::Kernel::Ptr, double>;

/******************************************************************************/

%{
#include "lsst/ip/diffim/ImageSubtract.h"
%}

%include "lsst/ip/diffim/ImageSubtract.h"

%template(fitSpatialKernelFromCandidates) lsst::ip::diffim::fitSpatialKernelFromCandidates<float>;

%template(pair_Kernel_Function) std::pair<lsst::afw::math::LinearCombinationKernel::Ptr, 
                                          lsst::afw::math::Kernel::SpatialFunctionPtr>;

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
%convolveAndSubtract(double);  // image subtraction on double images??!?  Not instantiated in .cc either
#endif

/******************************************************************************/

%{
#include "lsst/ip/diffim/KernelSolution.h"
%}


%define %KernelSolutionPtrs(NAME, TYPE)
SWIG_SHARED_PTR_DERIVED(StaticKernelSolution##NAME,
                        lsst::ip::diffim::KernelSolution,
                        lsst::ip::diffim::StaticKernelSolution<TYPE>);
SWIG_SHARED_PTR_DERIVED(MaskedKernelSolution##NAME,
                        lsst::ip::diffim::StaticKernelSolution<TYPE>,
                        lsst::ip::diffim::MaskedKernelSolution<TYPE>);
SWIG_SHARED_PTR_DERIVED(RegularizedKernelSolution##NAME,
                        lsst::ip::diffim::StaticKernelSolution<TYPE>,
                        lsst::ip::diffim::RegularizedKernelSolution<TYPE>);
%enddef

%define %KernelCandidates(NAME, TYPE)
%template(StaticKernelSolution##NAME) lsst::ip::diffim::StaticKernelSolution<TYPE>;
%template(MaskedKernelSolution##NAME) lsst::ip::diffim::MaskedKernelSolution<TYPE>;
%template(RegularizedKernelSolution##NAME) lsst::ip::diffim::RegularizedKernelSolution<TYPE>;
%enddef

SWIG_SHARED_PTR(KernelSolution, lsst::ip::diffim::KernelSolution);
SWIG_SHARED_PTR(StaticKernelSolution, lsst::ip::diffim::StaticKernelSolution);
SWIG_SHARED_PTR(SpatialKernelSolution, lsst::ip::diffim::SpatialKernelSolution);
%KernelSolutionPtrs(F, float);

%include "lsst/ip/diffim/KernelSolution.h"

%KernelCandidates(F, float);

/******************************************************************************/

%{
#include "lsst/ip/diffim/KernelCandidateDetection.h"
%}

%define %KernelCandidateDetectionPtr(NAME, TYPE)
SWIG_SHARED_PTR(KernelCandidateDetection##NAME,
                lsst::ip::diffim::KernelCandidateDetection<TYPE>);
%enddef

%define %KernelCandidateDetection(NAME, TYPE)
%template(KernelCandidateDetection##NAME) lsst::ip::diffim::KernelCandidateDetection<TYPE>;
%enddef

%KernelCandidateDetectionPtr(F, float);

%include "lsst/ip/diffim/KernelCandidateDetection.h"

%KernelCandidateDetection(F, float);
/******************************************************************************/

%{
#include "lsst/ip/diffim/KernelCandidate.h"
%}

%define %IMAGE(PIXTYPE)
lsst::afw::image::Image<PIXTYPE>
%enddef

%define %KernelCandidatePtr(NAME, TYPE)
SWIG_SHARED_PTR_DERIVED(KernelCandidate##NAME,
                        lsst::afw::math::SpatialCellImageCandidate<%IMAGE(lsst::afw::math::Kernel::Pixel)>,
                        lsst::ip::diffim::KernelCandidate<TYPE>);

/* Same problem as with meas algorithms makePsfCandidate */
%inline %{
namespace lsst { namespace ip { namespace diffim { namespace lsstSwig {
template <typename PixelT>
typename KernelCandidate<PixelT>::Ptr
makeKernelCandidateForSwig(float const xCenter,
                           float const yCenter, 
                           boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > const& 
                               miToConvolvePtr,
                           boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > const& 
                               miToNotConvolvePtr,
                           lsst::pex::policy::Policy const& policy) {
    
    return typename KernelCandidate<PixelT>::Ptr(new KernelCandidate<PixelT>(xCenter, yCenter,
                                                                             miToConvolvePtr,
                                                                             miToNotConvolvePtr,
                                                                             policy));
}
}}}}
%}

%ignore makeKernelCandidate;
%enddef

%define %KernelCandidate(NAME, TYPE)
%template(KernelCandidate##NAME) lsst::ip::diffim::KernelCandidate<TYPE>;
%template(makeKernelCandidate) lsst::ip::diffim::lsstSwig::makeKernelCandidateForSwig<TYPE>;
%inline %{
    lsst::ip::diffim::KernelCandidate<TYPE>::Ptr
        cast_KernelCandidate##NAME(lsst::afw::math::SpatialCellCandidate::Ptr candidate) {
        return boost::shared_dynamic_cast<lsst::ip::diffim::KernelCandidate<TYPE> >(candidate);
    }
%}
%enddef

%KernelCandidatePtr(F, float);

%include "lsst/ip/diffim/KernelCandidate.h"

%KernelCandidate(F, float);

/******************************************************************************/

%{
#include "lsst/ip/diffim/BasisLists.h"
%}

%include "lsst/ip/diffim/BasisLists.h"

/******************************************************************************/
/* I shouldn't have to do this but it exists noplace else, so... */

//%{
//#include "Eigen/Core"
//%}
//
//%template(eigenMatrixPtr) boost::shared_ptr<Eigen::MatrixXd>;
//%template(eigenVectorPtr) boost::shared_ptr<Eigen::VectorXd>;

/******************************************************************************/

%include "lsst/ip/diffim/detail.i"

