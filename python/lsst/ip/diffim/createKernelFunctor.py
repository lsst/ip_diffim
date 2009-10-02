import lsst.pex.policy as pexPolicy
import lsst.afw.math.mathLib as afwMath
import diffimLib 

def createKernelFunctor(policy):
    kCols = policy.getInt("kernelCols")
    kRows = policy.getInt("kernelRows")
    
    useAlard = policy.getBool("useAlardKernel")
    if useAlard:
        alardNGauss   = policy.getInt("alardNGauss")
        alardSigGauss = policy.getDoubleArray("alardSigGauss")
        alardDegGauss = policy.getIntArray("alardDegGauss")
        
        assert len(alardSigGauss) == alardNGauss
        assert len(alardDegGauss) == alardNGauss
        assert kCols == kRows          # square
        assert kCols % 2 == 1          # odd sized
        
        kHalfWidth = kCols//2
        basisList  = diffimLib.generateAlardLuptonKernelSet(kHalfWidth, alardNGauss, alardSigGauss, alardDegGauss)
        kFunctor   = diffimLib.PsfMatchingFunctorF(basisList)
        return kFunctor

    # else we use a delta function basis set
    basisList = diffimLib.generateDeltaFunctionKernelSet(kCols, kRows)
    useRegularization = policy.getBool("useRegularization")
    if useRegularization:
        regularizationOrder      = policy.getInt('regularizationOrder')
        regularizationBoundary   = policy.getInt('regularizationBoundary')
        regularizationDifference = policy.getInt('regularizationDifference')
        H = diffimLib.generateFiniteDifferenceRegularization(kCols, kRows,
                                                             regularizationOrder,
                                                             regularizationBoundary,
                                                             regularizationDifference)

        kFunctor = diffimLib.PsfMatchingFunctorF(basisList, H)
        return kFunctor

    kFunctor = diffimLib.PsfMatchingFunctorF(basisList)
    return kFunctor

    
