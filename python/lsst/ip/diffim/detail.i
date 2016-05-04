// -*- lsst-c++ -*-

/******************************************************************************/

%{
#include "lsst/ip/diffim/KernelSumVisitor.h"
%}

%define %KernelSumVisitorPtr(TYPE)
    %shared_ptr(lsst::ip::diffim::detail::KernelSumVisitor<TYPE>);
%enddef

%define %KernelSumVisitor(NAME, TYPE)
    %template(KernelSumVisitor##NAME) lsst::ip::diffim::detail::KernelSumVisitor<TYPE>;
    %template(makeKernelSumVisitor) lsst::ip::diffim::detail::makeKernelSumVisitor<TYPE>;
%enddef

%KernelSumVisitorPtr(float)

%include "lsst/ip/diffim/KernelSumVisitor.h"

%KernelSumVisitor(F, float)

/******************************************************************************/

%{
#include "lsst/ip/diffim/KernelPca.h"
#include <memory>
%}

%define %KernelPcaVisitorPtr(TYPE)
    %shared_ptr(lsst::ip::diffim::detail::KernelPcaVisitor<TYPE>);
%enddef

%define %KernelPcaPtr(TYPE)
    %shared_ptr(lsst::ip::diffim::detail::KernelPca<lsst::afw::image::Image<TYPE> >);
%enddef

%define %KernelPcaVisitor(NAME, TYPE)
    %template(KernelPcaVisitor##NAME) lsst::ip::diffim::detail::KernelPcaVisitor<TYPE>;
    %template(makeKernelPcaVisitor) lsst::ip::diffim::detail::makeKernelPcaVisitor<TYPE>;
%enddef

%define %KernelPca(NAME, TYPE)
    %template(KernelPca##NAME) lsst::ip::diffim::detail::KernelPca<lsst::afw::image::Image<TYPE> >;
    %inline %{
        lsst::ip::diffim::detail::KernelPca<lsst::afw::image::Image<TYPE> >::Ptr
            cast_KernelPca##NAME(lsst::afw::image::ImagePca<lsst::afw::image::Image<TYPE> >::Ptr imagePca) {
            return std::dynamic_pointer_cast<lsst::ip::diffim::detail::KernelPca<lsst::afw::image::Image<TYPE> > >(imagePca);
        }
    %}
%enddef

%KernelPcaPtr(lsst::afw::math::Kernel::Pixel)
%KernelPcaVisitorPtr(float)

%include "lsst/ip/diffim/KernelPca.h"

%KernelPca(D, lsst::afw::math::Kernel::Pixel)
%KernelPcaVisitor(F, float)


/******************************************************************************/

%{
#include "lsst/ip/diffim/BuildSingleKernelVisitor.h"
%}

%define %BuildSingleKernelVisitorPtr(TYPE)
    %shared_ptr(lsst::ip::diffim::detail::BuildSingleKernelVisitor<TYPE>);
%enddef

%define %BuildSingleKernelVisitor(NAME, TYPE)
    %template(BuildSingleKernelVisitor##NAME) lsst::ip::diffim::detail::BuildSingleKernelVisitor<TYPE>;
    %template(makeBuildSingleKernelVisitor) lsst::ip::diffim::detail::makeBuildSingleKernelVisitor<TYPE>;
%enddef

%BuildSingleKernelVisitorPtr(float)

%include "lsst/ip/diffim/BuildSingleKernelVisitor.h"

%BuildSingleKernelVisitor(F, float)

/******************************************************************************/

%{
#include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"
%}

%define %BuildSpatialKernelVisitorPtr(TYPE)
    %shared_ptr(lsst::ip::diffim::detail::BuildSpatialKernelVisitor<TYPE>);
%enddef

%define %BuildSpatialKernelVisitor(NAME, TYPE)
    %template(BuildSpatialKernelVisitor##NAME) lsst::ip::diffim::detail::BuildSpatialKernelVisitor<TYPE>;
    %template(makeBuildSpatialKernelVisitor) lsst::ip::diffim::detail::makeBuildSpatialKernelVisitor<TYPE>;
%enddef

%BuildSpatialKernelVisitorPtr(float)

%include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"

%BuildSpatialKernelVisitor(F, float)


/******************************************************************************/

%{
#include "lsst/ip/diffim/AssessSpatialKernelVisitor.h"
%}

%define %AssessSpatialKernelVisitorPtr(TYPE)
    %shared_ptr(lsst::ip::diffim::detail::AssessSpatialKernelVisitor<TYPE>);
%enddef

%define %AssessSpatialKernelVisitor(NAME, TYPE)
    %template(AssessSpatialKernelVisitor##NAME) lsst::ip::diffim::detail::AssessSpatialKernelVisitor<TYPE>;
    %template(makeAssessSpatialKernelVisitor) lsst::ip::diffim::detail::makeAssessSpatialKernelVisitor<TYPE>;
%enddef

%AssessSpatialKernelVisitorPtr(float)

%include "lsst/ip/diffim/AssessSpatialKernelVisitor.h"

%AssessSpatialKernelVisitor(F, float)
