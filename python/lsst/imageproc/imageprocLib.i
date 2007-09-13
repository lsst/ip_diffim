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

using namespace lsst::fw;
%}

%init %{
%}

%include "lsst/mwi/p_lsstSwig.i"
%include "lsst/fw/Core/lsstImageTypes.i"     // vw and Image/Mask types and typedefs

%include "lsst/imageproc/ImageSubtract.h"
%include "lsst/imageproc/PCA.h"
// Until you merge this in to fw trunk
%include "lsst/fw/MinimizerFunctionBase.h"


/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
