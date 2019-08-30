
from .assessSpatialKernelVisitor import *
from .buildSingleKernelVisitor import *
from .buildSpatialKernelVisitor import *
from .kernelPca import *
from .kernelSumVisitor import *

from ..deprecated import deprecate_policy as _deprecate_policy


AssessSpatialKernelVisitorF = _deprecate_policy(AssessSpatialKernelVisitorF, policy_args=[2])
BuildSingleKernelVisitorF = _deprecate_policy(BuildSingleKernelVisitorF, policy_args=[1])
BuildSpatialKernelVisitorF = _deprecate_policy(BuildSpatialKernelVisitorF, policy_args=[1])
KernelSumVisitorF = _deprecate_policy(KernelSumVisitorF)

makeKernelSumVisitor = _deprecate_policy(makeKernelSumVisitor)
