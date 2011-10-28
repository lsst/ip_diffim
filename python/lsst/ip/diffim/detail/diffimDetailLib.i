// -*- lsst-c++ -*-
%define diffimDetailLib_DOCSTRING
"
Python bindings for lsst::ip::diffim::detail code
"
%enddef

%feature("autodoc", "1");
%module(docstring=diffimDetailLib_DOCSTRING, package="lsst.ip.diffim.detail") diffimDetailLib

//%{
//#include <boost/shared_ptr.hpp>
//#include <lsst/pex/policy/Policy.h>
//#include <lsst/afw.h>
//#include <lsst/ip/diffim.h>
//%}

%include "lsst/p_lsstSwig.i"
%import  "lsst/afw/image/imageLib.i"
%import  "lsst/afw/math/mathLib.i"
//%lsst_exceptions();

%{
#include "lsst/ip/diffim/AssessSpatialKernelVisitor.h"
#include "lsst/ip/diffim/BuildSingleKernelVisitor.h"
#include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"
#include "lsst/ip/diffim/KernelPcaVisitor.h"
#include "lsst/ip/diffim/KernelSumVisitor.h"
%}

SWIG_SHARED_PTR(AssessSpatialKernelVisitor, lsst::ip::diffim::detail::AssessSpatialKernelVisitor);
SWIG_SHARED_PTR(BuildSingleKernelVisitor, lsst::ip::diffim::detail::BuildSingleKernelVisitor);
SWIG_SHARED_PTR(BuildSpatialKernelVisitor, lsst::ip::diffim::detail::BuildSpatialKernelVisitor);
SWIG_SHARED_PTR(KernelPcaVisitor, lsst::ip::diffim::detail::KernelPcaVisitor);
SWIG_SHARED_PTR(KernelSumVisitorF, lsst::ip::diffim::detail::KernelSumVisitor<float>);

%include "lsst/ip/diffim/AssessSpatialKernelVisitor.h"
%include "lsst/ip/diffim/BuildSingleKernelVisitor.h"
%include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"
%include "lsst/ip/diffim/KernelPcaVisitor.h"
%include "lsst/ip/diffim/KernelSumVisitor.h"

%template(KernelSumVisitorF) lsst::ip::diffim::detail::KernelSumVisitor<float>;
