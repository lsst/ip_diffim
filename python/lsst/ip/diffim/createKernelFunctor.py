import lsst.pex.logging as pexLogging
import diffimLib 

def createKernelFunctor(policy):
    kCols = policy.getInt("kernelCols")
    kRows = policy.getInt("kernelRows")
    kType = policy.getString("kernelBasisSet")
    useRegularization = policy.getBool("useRegularization")
    
    if kType == "alard-lupton":
        alardNGauss   = policy.getInt("alardNGauss")
        alardSigGauss = policy.getDoubleArray("alardSigGauss")
        alardDegGauss = policy.getIntArray("alardDegGauss")

        if useRegularization:
            pexLogging.Trace("lsst.ip.diffim.createKernelFunctor", 1,
                             "Warning : Regularization not enabled for Alard-Lupton kernels")

        try:
            assert len(alardSigGauss) == alardNGauss
            assert len(alardDegGauss) == alardNGauss
            assert kCols == kRows          # square
            assert kCols % 2 == 1          # odd sized
        except:
            raise RuntimeError("Invalid Alard-Lupton kernel configuration")
        
        kHalfWidth = kCols//2
        basisList  = diffimLib.generateAlardLuptonBasisSet(kHalfWidth,
                                                           alardNGauss,
                                                           alardSigGauss,
                                                           alardDegGauss)
        kFunctor   = diffimLib.PsfMatchingFunctorF(basisList)
        return kFunctor

    elif kType == "delta-function":
        basisList = diffimLib.generateDeltaFunctionBasisSet(kCols, kRows)
        
        if useRegularization:
            regularizationOrder      = policy.getInt("regularizationOrder")
            regularizationBoundary   = policy.getInt("regularizationBoundary")
            regularizationDifference = policy.getInt("regularizationDifference")
            h = diffimLib.generateFiniteDifferenceRegularization(kCols, kRows,
                                                                 regularizationOrder,
                                                                 regularizationBoundary,
                                                                 regularizationDifference)

            kFunctor = diffimLib.PsfMatchingFunctorF(basisList, h)
            return kFunctor

        kFunctor = diffimLib.PsfMatchingFunctorF(basisList)
        return kFunctor

    else:
        raise RuntimeError("Invalid kernel type : %s" % (kType))
    
