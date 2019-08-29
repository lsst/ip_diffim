
from .assessSpatialKernelVisitor import *
from .buildSingleKernelVisitor import *
from .buildSpatialKernelVisitor import *
from .kernelPca import *
from .kernelSumVisitor import *

from ..deprecated import deprecate_policy as _deprecate_policy


AssessSpatialKernelVisitorF = _deprecate_policy(AssessSpatialKernelVisitorF)
BuildSingleKernelVisitorF = _deprecate_policy(BuildSingleKernelVisitorF)
BuildSpatialKernelVisitorF = _deprecate_policy(BuildSpatialKernelVisitorF)
KernelSumVisitorF = _deprecate_policy(KernelSumVisitorF)

makeKernelSumVisitor = _deprecate_policy(makeKernelSumVisitor)
