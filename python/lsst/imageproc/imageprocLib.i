// -*- lsst-c++ -*-
%define imageprocLib_DOCSTRING
"
Python bindings for imageproc module
"
%enddef

%feature("autodoc", "1");
%module(docstring=imageprocLib_DOCSTRING) imageprocLib

%{
#   include <lsst/imageproc/ImageSubtract.h>
#   include <lsst/detection/Footprint.h>
%}

%inline %{
namespace lsst { namespace fw { } }
namespace lsst { namespace imageproc { } }
namespace lsst { namespace detection { } }

using namespace lsst;
using namespace lsst::imageproc;
%}

%init %{
%}

%include "lsst/mwi/p_lsstSwig.i"
%include "lsst/fw/Core/lsstImageTypes.i"     // vw and Image/Mask types and typedefs

%pythoncode %{
def version(HeadURL = r"$HeadURL: svn+ssh://svn.lsstcorp.org/DC2/imageproc/tickets/7/python/lsst/imageproc/imageprocLib.i $"):
    """Return a version given a HeadURL string; default: imageproc's version"""
    return guessSvnVersion(HeadURL)

%}


// Here you give the names of the functions you want to Swig

/*******************/
/* ImageSubtract.h */

%include "lsst/imageproc/ImageSubtract.h"
%template(computePsfMatchingKernelForMaskedImage_DIDD) lsst::imageproc::computePsfMatchingKernelForMaskedImage<double, int32, double, double>;
%template(computePsfMatchingKernelForMaskedImage_DiDD) lsst::imageproc::computePsfMatchingKernelForMaskedImage<double, int16, double, double>;
%template(computePsfMatchingKernelForMaskedImage_FIDD) lsst::imageproc::computePsfMatchingKernelForMaskedImage<float,  int32, double, double>;
%template(computePsfMatchingKernelForMaskedImage_FiDD) lsst::imageproc::computePsfMatchingKernelForMaskedImage<float,  int16, double, double>;

%template(computePsfMatchingKernelForPostageStamp_DID) lsst::imageproc::computePsfMatchingKernelForPostageStamp<double, int32, double>;
%template(computePsfMatchingKernelForPostageStamp_DiD) lsst::imageproc::computePsfMatchingKernelForPostageStamp<double, int16, double>;
%template(computePsfMatchingKernelForPostageStamp_FID) lsst::imageproc::computePsfMatchingKernelForPostageStamp<float,  int32, double>;
%template(computePsfMatchingKernelForPostageStamp_FiD) lsst::imageproc::computePsfMatchingKernelForPostageStamp<float,  int16, double>;

%template(getCollectionOfFootprintsForPsfMatching_DI)  lsst::imageproc::getCollectionOfFootprintsForPsfMatching<double, int32>;
%template(getCollectionOfFootprintsForPsfMatching_Di)  lsst::imageproc::getCollectionOfFootprintsForPsfMatching<double, int16>;
%template(getCollectionOfFootprintsForPsfMatching_FI)  lsst::imageproc::getCollectionOfFootprintsForPsfMatching<float,  int32>;
%template(getCollectionOfFootprintsForPsfMatching_Fi)  lsst::imageproc::getCollectionOfFootprintsForPsfMatching<float,  int16>;

// getCollectionOfMaskedImagesForPsfMatching is not templated

%template(computePcaKernelBasis_D)                     lsst::imageproc::computePcaKernelBasis<double>;
%template(computeSpatiallyVaryingPsfMatchingKernel_DD) lsst::imageproc::computeSpatiallyVaryingPsfMatchingKernel<double, double>;
%template(generateDeltaFunctionKernelSet_D)            lsst::imageproc::generateDeltaFunctionKernelSet<double>;
%template(generateAlardLuptonKernelSet_D)              lsst::imageproc::generateAlardLuptonKernelSet<double>;

%template(checkMaskedImageForDiffim_DI)                lsst::imageproc::checkMaskedImageForDiffim<double, int32>;
%template(checkMaskedImageForDiffim_Di)                lsst::imageproc::checkMaskedImageForDiffim<double, int16>;
%template(checkMaskedImageForDiffim_FI)                lsst::imageproc::checkMaskedImageForDiffim<float, int32>;
%template(checkMaskedImageForDiffim_Fi)                lsst::imageproc::checkMaskedImageForDiffim<float, int16>;

%template(calculateMaskedImageResiduals_DI)            lsst::imageproc::calculateMaskedImageResiduals<double, int32>;
%template(calculateMaskedImageResiduals_Di)            lsst::imageproc::calculateMaskedImageResiduals<double, int16>;
%template(calculateMaskedImageResiduals_FI)            lsst::imageproc::calculateMaskedImageResiduals<float,  int32>;
%template(calculateMaskedImageResiduals_Fi)            lsst::imageproc::calculateMaskedImageResiduals<float,  int16>;

%template(calculateImageResiduals_D)                   lsst::imageproc::calculateImageResiduals<double>;
%template(calculateImageResiduals_F)                   lsst::imageproc::calculateImageResiduals<float>;

/* ImageSubtract.h */
/*******************/
/* PCA.h */

#include "vw/Math/Matrix.h" 
#include "vw/Math/Vector.h"
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
