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

%include "lsst/imageproc/ImageSubtract.h"
%include "lsst/imageproc/PCA.h"
// Until you merge this in to fw trunk
%include "lsst/fw/MinimizerFunctionBase.h"


/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
