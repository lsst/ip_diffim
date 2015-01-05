// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

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
}

/* handle C++ arguments that should be outputs in python */
%apply int& OUTPUT { int& };
%apply float& OUTPUT { float& };
%apply double& OUTPUT { double& };

%{
#include "boost/shared_ptr.hpp"


#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/math.h"
#include "lsst/afw/geom.h" 
#include "lsst/afw/image.h"
#include "lsst/afw/table.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/meas/base.h"

#include "lsst/ip/diffim.h"
%}

/******************************************************************************/

/* Eigen / numpy.  Info comes from $AFW_DIR/include/lsst/afw/numpyTypemaps.h */
%{
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%init %{
    import_array();
%}

%include "lsst/p_lsstSwig.i"
%include "ndarray.i"

%template(PtrEigenMatrixXd) boost::shared_ptr<Eigen::MatrixXd>;
%template(PtrEigenVectorXd) boost::shared_ptr<Eigen::VectorXd>;
%declareNumPyConverters(Eigen::MatrixXd);
%declareNumPyConverters(Eigen::VectorXd);
%shared_ptr(Eigen::MatrixXd);
%shared_ptr(Eigen::VectorXd);

%include "lsst/daf/base/persistenceMacros.i"
%include "lsst/pex/config.h"            // LSST_CONTROL_FIELD.
%include "lsst/meas/base/constants.h"
%include "lsst/meas/base/exceptions.i"
%include "lsst/meas/base/utilities.i"
%include "lsst/meas/base/Algorithm.h"

%import  "lsst/afw/image/imageLib.i"
%import  "lsst/afw/math/mathLib.i"
%import  "lsst/afw/math/objectVectors.i"
%import  "lsst/afw/detection/detectionLib.i"
%include "lsst/afw/image/Exposure.h"
%include "lsst/afw/table/Source.h"
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
%template(pair_Kernel_Function) std::pair<lsst::afw::math::LinearCombinationKernel::Ptr,
                                          lsst::afw::math::Kernel::SpatialFunctionPtr>;

/******************************************************************************/

%{
#include "lsst/ip/diffim/FindSetBits.h"
%}

%include "lsst/ip/diffim/FindSetBits.h"

%template(FindSetBitsU)
    lsst::ip::diffim::FindSetBits<lsst::afw::image::Mask<lsst::afw::image::MaskPixel> >;

/******************************************************************************/

%{
#include "lsst/ip/diffim/ImageStatistics.h"
%}

%include "lsst/ip/diffim/ImageStatistics.h"

%template(ImageStatisticsI)
    lsst::ip::diffim::ImageStatistics<int>;
%template(ImageStatisticsF)
    lsst::ip::diffim::ImageStatistics<float>;
%template(ImageStatisticsD)
    lsst::ip::diffim::ImageStatistics<double>;

/******************************************************************************/

%{
#include "lsst/ip/diffim/ImageSubtract.h"
%}

%include "lsst/ip/diffim/ImageSubtract.h"


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
%shared_ptr(lsst::ip::diffim::StaticKernelSolution<TYPE>);
%shared_ptr(lsst::ip::diffim::MaskedKernelSolution<TYPE>);
%shared_ptr(lsst::ip::diffim::RegularizedKernelSolution<TYPE>);
%enddef

%define %KernelSolutions(NAME, TYPE)
%template(StaticKernelSolution##NAME) lsst::ip::diffim::StaticKernelSolution<TYPE>;
%template(MaskedKernelSolution##NAME) lsst::ip::diffim::MaskedKernelSolution<TYPE>;
%template(RegularizedKernelSolution##NAME) lsst::ip::diffim::RegularizedKernelSolution<TYPE>;
%enddef

%shared_ptr(lsst::ip::diffim::KernelSolution);
%shared_ptr(lsst::ip::diffim::SpatialKernelSolution);

%KernelSolutionPtrs(F, float);

%include "lsst/ip/diffim/KernelSolution.h"

%KernelSolutions(F, float);

/******************************************************************************/

%{
#include "lsst/ip/diffim/KernelCandidateDetection.h"
%}

%define %KernelCandidateDetectionPtr(NAME, TYPE)
%shared_ptr(lsst::ip::diffim::KernelCandidateDetection<TYPE>);
%enddef

%define %KernelCandidateDetection(NAME, TYPE)
%template(KernelCandidateDetection##NAME) lsst::ip::diffim::KernelCandidateDetection<TYPE>;
%enddef

%KernelCandidateDetectionPtr(F, float);

%include "lsst/ip/diffim/KernelCandidateDetection.h"

%KernelCandidateDetection(F, float);
/******************************************************************************/

%define %IMAGE(PIXTYPE)
lsst::afw::image::Image<PIXTYPE>
%enddef

%define %KernelCandidatePtr(NAME, TYPE)
%shared_ptr(lsst::ip::diffim::KernelCandidate<TYPE>);

%enddef

%define %KernelCandidate(NAME, TYPE)
%template(KernelCandidate##NAME) lsst::ip::diffim::KernelCandidate<TYPE>;
%template(makeKernelCandidate) lsst::ip::diffim::makeKernelCandidate<TYPE>;
%inline %{
    lsst::ip::diffim::KernelCandidate<TYPE>::Ptr
        cast_KernelCandidate##NAME(lsst::afw::math::SpatialCellCandidate::Ptr candidate) {
        return boost::dynamic_pointer_cast<lsst::ip::diffim::KernelCandidate<TYPE> >(candidate);
    }
%}
%enddef

%include "lsst/ip/diffim/KernelCandidate.h"
%KernelCandidatePtr(F, float);

%include "lsst/ip/diffim/KernelCandidate.h"

%KernelCandidate(F, float);

/******************************************************************************/

%{
#include "lsst/ip/diffim/BasisLists.h"
%}

%include "lsst/ip/diffim/BasisLists.h"


%include  "lsst/ip/diffim/dipole.i"

/******************************************************************************/

%include "lsst/ip/diffim/detail.i"

/******************************************************************************/
