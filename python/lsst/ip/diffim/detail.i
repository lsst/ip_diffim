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

%{
#include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"
%}

%define %BuildSpatialKernelVisitorPtr(NAME, TYPE)
SWIG_SHARED_PTR_DERIVED(BuildSpatialKernelVisitor##NAME, 
                        lsst::afw::math::CandidateVisitor, 
                        lsst::ip::diffim::detail::BuildSpatialKernelVisitor<TYPE>);
%enddef

%define %BuildSpatialKernelVisitor(NAME, TYPE)
%template(BuildSpatialKernelVisitor##NAME) lsst::ip::diffim::detail::BuildSpatialKernelVisitor<TYPE>;
%template(makeBuildSpatialKernelVisitor) lsst::ip::diffim::detail::makeBuildSpatialKernelVisitor<TYPE>;
%enddef

%BuildSpatialKernelVisitorPtr(F, float)

%include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"

%BuildSpatialKernelVisitor(F, float)


/******************************************************************************/

%{
#include "lsst/ip/diffim/AssessSpatialKernelVisitor.h"
%}

%define %AssessSpatialKernelVisitorPtr(NAME, TYPE)
SWIG_SHARED_PTR_DERIVED(AssessSpatialKernelVisitor##NAME, 
                        lsst::afw::math::CandidateVisitor, 
                        lsst::ip::diffim::detail::AssessSpatialKernelVisitor<TYPE>);
%enddef

%define %AssessSpatialKernelVisitor(NAME, TYPE)
%template(AssessSpatialKernelVisitor##NAME) lsst::ip::diffim::detail::AssessSpatialKernelVisitor<TYPE>;
%template(makeAssessSpatialKernelVisitor) lsst::ip::diffim::detail::makeAssessSpatialKernelVisitor<TYPE>;
%enddef

%AssessSpatialKernelVisitorPtr(F, float)

%include "lsst/ip/diffim/AssessSpatialKernelVisitor.h"

%AssessSpatialKernelVisitor(F, float)
