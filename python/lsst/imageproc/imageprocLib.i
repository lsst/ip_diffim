// -*- lsst-c++ -*-
%define imageprocLib_DOCSTRING
"
Python bindings for imageproc module
"
%enddef

%feature("autodoc", "1");
%module(docstring=imageprocLib_DOCSTRING) imageprocLib

// Everything we will need in the _wrap.cc file
%{
#include <lsst/imageproc/ImageSubtract.h>
#include <lsst/detection/Footprint.h>
%}

%inline %{
namespace lsst { namespace fw { } }
namespace lsst { namespace imageproc { } }

using namespace lsst;
using namespace lsst::imageproc;
%}

%init %{
%}

// Everything whose bindings we will have to know about
%include "lsst/mwi/p_lsstSwig.i"             // this needs to go first otherwise i don't know about e.g. boost
%include "lsst/fw/Core/lsstImageTypes.i"     // vw and Image/Mask types and typedefs
%include "lsst/detection/detectionLib.i"     // need otherwise FootprintContainerT not known about
%include "cpointer.i"

%pythoncode %{
def version(HeadURL = r"$HeadURL: svn+ssh://svn.lsstcorp.org/DC2/imageproc/tickets/7/python/lsst/imageproc/imageprocLib.i $"):
    """Return a version given a HeadURL string; default: imageproc's version"""
    return guessSvnVersion(HeadURL)

%}

// Here you give the names of the functions you want to Swig

/*******************/
/* ImageSubtract.h */

// this should be in kernel.i but i need it now
%template(LinearCombinationKernelPtrTypeD) boost::shared_ptr<lsst::fw::LinearCombinationKernel<double> >;

%include "lsst/imageproc/ImageSubtract.h"
%template(DiffImContainer_D)                            lsst::imageproc::DiffImContainer<double>;
%template(vectorDiffImContainer_D)                      std::vector<lsst::imageproc::DiffImContainer<double> >;
%template(computePsfMatchingKernelForMaskedImage_FU8DD) lsst::imageproc::computePsfMatchingKernelForMaskedImage<float, 
                                                                                                                lsst::fw::maskPixelType, 
                                                                                                                double, 
                                                                                                                double>;
%template(computePsfMatchingKernelForPostageStamp_FU8D) lsst::imageproc::computePsfMatchingKernelForPostageStamp<float, 
                                                                                                                 lsst::fw::maskPixelType, 
                                                                                                                 double>;
%template(getCollectionOfFootprintsForPsfMatching_FU8)  lsst::imageproc::getCollectionOfFootprintsForPsfMatching<float, 
                                                                                                                 lsst::fw::maskPixelType>;
%template(computePcaKernelBasis_D)                      lsst::imageproc::computePcaKernelBasis<double>;
%template(computeSpatiallyVaryingPsfMatchingKernel_DD)  lsst::imageproc::computeSpatiallyVaryingPsfMatchingKernel<double, double>;
%template(generateDeltaFunctionKernelSet_D)             lsst::imageproc::generateDeltaFunctionKernelSet<double>;
%template(generateAlardLuptonKernelSet_D)               lsst::imageproc::generateAlardLuptonKernelSet<double>;
%template(checkMaskedImageForDiffim_FU8)                lsst::imageproc::checkMaskedImageForDiffim<float, lsst::fw::maskPixelType>;
%template(calculateMaskedImageResiduals_FU8)            lsst::imageproc::calculateMaskedImageResiduals<float, lsst::fw::maskPixelType>;
%template(calculateImageResiduals_F)                    lsst::imageproc::calculateImageResiduals<float>;
// getCollectionOfMaskedImagesForPsfMatching is not templated

/* ImageSubtract.h */
/*******************/
/* PCA.h */

%include "lsst/imageproc/PCA.h"
%template(computePca_F)                  lsst::imageproc::computePca<vw::math::Matrix<float>, vw::math::Vector<float> >;
%template(computePca_D)                  lsst::imageproc::computePca<vw::math::Matrix<double>, vw::math::Vector<double> >;
%template(decomposeMatrixUsingBasis_F)   lsst::imageproc::decomposeMatrixUsingBasis<vw::math::Matrix<float> >;
%template(decomposeMatrixUsingBasis_D)   lsst::imageproc::decomposeMatrixUsingBasis<vw::math::Matrix<double> >;
%template(approximateMatrixUsingBasis_F) lsst::imageproc::approximateMatrixUsingBasis<vw::math::Matrix<float> >;
%template(approximateMatrixUsingBasis_D) lsst::imageproc::approximateMatrixUsingBasis<vw::math::Matrix<double> >;

/* PCA.h */
/*******************/
/* MinimizerFunctionBase.h */

// Until we merge this in to fw trunk
%include "Minuit/FCNBase.h"
%include "lsst/fw/MinimizerFunctionBase.h"
%template(MinimizerFunctionBase1_F) lsst::fw::function::MinimizerFunctionBase1<float>;
%template(MinimizerFunctionBase1_D) lsst::fw::function::MinimizerFunctionBase1<double>;
%template(MinimizerFunctionBase2_F) lsst::fw::function::MinimizerFunctionBase2<float>;
%template(MinimizerFunctionBase2_D) lsst::fw::function::MinimizerFunctionBase2<double>;
%template(minimize_F)               lsst::fw::function::minimize<float>;
%template(minimize_D)               lsst::fw::function::minimize<double>;

/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
