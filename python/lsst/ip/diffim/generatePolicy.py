import eups
import os
import lsst.pex.policy as pexPolicy
import lsst.pex.logging as pexLog

def generateDefaultPolicy(path, modify=True, fwhm=3.5):
    diffimDir    = eups.productDir('ip_diffim')
    diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')
    policy       = pexPolicy.Policy.createPolicy(diffimPolicy)
    if modify:
        modifyKernelPolicy(policy, fwhm)
    return policy

def modifyKernelPolicy(policy, fwhm=3.5):
    # Modify the kernel policy parameters based upon the images FWHM
    
    kernelRadius = policy.getDouble("kernelRadiusFwhmScaling") * fwhm
    kernelRadius = min(kernelRadius, policy.getInt("kernelRadiusMax"))
    kernelRadius = max(kernelRadius, policy.getInt("kernelRadiusMin"))
    kernelSize   = 2 * int(kernelRadius + 0.5) + 1
    policy.set("kernelRows", kernelSize)
    policy.set("kernelCols", kernelSize)
    
    fpGrow = policy.getDouble("fpGrowFwhmScaling") * fwhm
    fpGrow = min(fpGrow, policy.getInt("fpGrowMax"))
    fpGrow = max(fpGrow, policy.getInt("fpGrowMin"))
    policy.set("fpGrowPix", int(fpGrow))

    alardSig = [x*fwhm for x in policy.getDoubleArray("alardSigFwhmScaling")]
    policy.set("alardSigGauss", alardSig[0])
    for sig in alardSig[1:]:
        policy.add("alardSigGauss", sig)
    
    pexLog.Trace("lsst.ip.diffim.generateKernelPolicy", 2,
                 "Using Psf Fwhm     : %.2f px" % (fwhm))

    pexLog.Trace("lsst.ip.diffim.generateKernelPolicy", 2,
                 "Kernel rows/cols   : %d x %d px" % (policy.getInt("kernelRows"),
                                                     policy.getInt("kernelCols")))
    
    pexLog.Trace("lsst.ip.diffim.generateKernelPolicy", 2,
                 "Footprint grow rad : %d px" % (policy.getInt("fpGrowPix")))

    outStr = ", ".join(["%.2f"%(x) for x in policy.getDoubleArray("alardSigGauss")])
    pexLog.Trace("lsst.ip.diffim.generateKernelPolicy", 2,
                 "A/L gaussian sig   : %s px" % (outStr))


