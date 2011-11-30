// -*- lsst-c++ -*-
%define diffimDetailLib_DOCSTRING
"
Python bindings for lsst::ip::diffim::detail code
"
%enddef

%feature("autodoc", "1");
%module(docstring=diffimDetailLib_DOCSTRING, package="lsst.ip.diffim.detail") diffimDetailLib

%{
#include "lsst/pex/logging.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/ip/diffim/AssessSpatialKernelVisitor.h"
#include "lsst/ip/diffim/BuildSingleKernelVisitor.h"
#include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"
#include "lsst/ip/diffim/KernelPcaVisitor.h"
#include "lsst/ip/diffim/KernelSumVisitor.h"
%}

%include "lsst/p_lsstSwig.i"
%import  "lsst/afw/image/imageLib.i"
%import  "lsst/afw/math/mathLib.i"
%lsst_exceptions();

%shared_ptr(lsst::ip::diffim::detail::AssessSpatialKernelVisitor);
%shared_ptr(lsst::ip::diffim::detail::BuildSingleKernelVisitor);
%shared_ptr(lsst::ip::diffim::detail::BuildSpatialKernelVisitor);
%shared_ptr(lsst::ip::diffim::detail::KernelPcaVisitor);
%shared_ptr(lsst::ip::diffim::detail::KernelSumVisitor<float>);

%include "lsst/ip/diffim/AssessSpatialKernelVisitor.h"
%include "lsst/ip/diffim/BuildSingleKernelVisitor.h"
%include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"
%include "lsst/ip/diffim/KernelPcaVisitor.h"
%include "lsst/ip/diffim/KernelSumVisitor.h"

%template(KernelSumVisitorF) lsst::ip::diffim::detail::KernelSumVisitor<float>;
