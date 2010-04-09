// -*- lsst-c++ -*-

/******************************************************************************/

%{
#include "lsst/ip/diffim/KernelSumVisitor.h"
%}

%define %KernelSumVisitorPtr(NAME, TYPE)
SWIG_SHARED_PTR_DERIVED(KernelSumVisitor##NAME, 
                        lsst::afw::math::CandidateVisitor, 
                        lsst::ip::diffim::detail::KernelSumVisitor<TYPE>);
%enddef

%define %KernelSumVisitor(NAME, TYPE)
%template(KernelSumVisitor##NAME) lsst::ip::diffim::detail::KernelSumVisitor<TYPE>;
%template(makeKernelSumVisitor) lsst::ip::diffim::detail::makeKernelSumVisitor<TYPE>;
%enddef

%KernelSumVisitorPtr(F, float)

%include "lsst/ip/diffim/KernelSumVisitor.h"

%KernelSumVisitor(F, float)

/******************************************************************************/

%{
#include "lsst/ip/diffim/KernelPcaVisitor.h"
%}

%define %KernelPcaVisitorPtr(NAME, TYPE)
SWIG_SHARED_PTR_DERIVED(KernelPcaVisitor##NAME, 
                        lsst::afw::math::CandidateVisitor, 
                        lsst::ip::diffim::detail::KernelPcaVisitor<TYPE>);
%enddef

%define %KernelPcaVisitor(NAME, TYPE)
%template(KernelPcaVisitor##NAME) lsst::ip::diffim::detail::KernelPcaVisitor<TYPE>;
%template(makeKernelPcaVisitor) lsst::ip::diffim::detail::makeKernelPcaVisitor<TYPE>;
%enddef

%KernelPcaVisitorPtr(F, float)

%include "lsst/ip/diffim/KernelPcaVisitor.h"

%KernelPcaVisitor(F, float)

/******************************************************************************/

%{
#include "lsst/ip/diffim/BuildSingleKernelVisitor.h"
%}

%define %BuildSingleKernelVisitorPtr(NAME, TYPE)
SWIG_SHARED_PTR_DERIVED(BuildSingleKernelVisitor##NAME, 
                        lsst::afw::math::CandidateVisitor, 
                        lsst::ip::diffim::detail::BuildSingleKernelVisitor<TYPE>);
%enddef

%define %BuildSingleKernelVisitor(NAME, TYPE)
%template(BuildSingleKernelVisitor##NAME) lsst::ip::diffim::detail::BuildSingleKernelVisitor<TYPE>;
%template(makeBuildSingleKernelVisitor) lsst::ip::diffim::detail::makeBuildSingleKernelVisitor<TYPE>;
%enddef

%BuildSingleKernelVisitorPtr(F, float)

%include "lsst/ip/diffim/BuildSingleKernelVisitor.h"

%BuildSingleKernelVisitor(F, float)

/******************************************************************************/

