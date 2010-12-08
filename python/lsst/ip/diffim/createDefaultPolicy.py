import lsst.pex.policy as pexPolicy
import lsst.pex.logging as pexLog

import math
sigma2fwhm = 2. * math.sqrt(2. * math.log(2.))

# Here we assume we have the FWHM of the images' PSF in pixels.  We
# will generate an alard-lupton kernel policy based upon these values.
# A complication is that the "sigma" of the AL gaussians is related to
# the FWHM through a scaling of 2.35.

def createDefaultPolicy(policyPath, modify=True, fwhm=3.5):
    policy = pexPolicy.Policy.createPolicy(policyPath)
    if modify:
        modifyKernelPolicy(policy, fwhm)
    return policy

def modifyKernelPolicy(policy, fwhm=3.5):
    # Modify the kernel policy parameters based upon the images FWHM

    kernelSize  = policy.getDouble("kernelSizeFwhmScaling") * fwhm
    kernelSize  = min(kernelRadius, policy.getInt("kernelSizeMax"))
    kernelSize  = max(kernelRadius, policy.getInt("kernelSizeMin"))
    kernelSize += (1 - kernelSize % 2) # make odd sized
    policy.set("kernelSize", kernelSize)

    fpGrow = policy.getPolicy("detectionPolicy").getDouble("fpGrowFwhmScaling") * fwhm
    fpGrow = min(fpGrow, policy.getInt("fpGrowMax"))
    fpGrow = max(fpGrow, policy.getInt("fpGrowMin"))
    policy.getPolicy("detectionPolicy").set("fpGrowPix", int(fpGrow))

    sigma = fwhm / sigma2fwhm
    alardSig = [x*sigma for x in policy.getDoubleArray("alardSigFwhmScaling")]
    policy.set("alardSigGauss", alardSig[0])
    for sig in alardSig[1:]:
        policy.add("alardSigGauss", sig)
    
    pexLog.Trace("lsst.ip.diffim.createDefaultPolicy.modifyKernelPolicy", 2,
                 "Using Psf Fwhm     : %.2f px" % (fwhm))

    pexLog.Trace("lsst.ip.diffim.createDefaultPolicy.modifyKernelPolicy", 2,
                 "Kernel size   : %d x %d px" % (policy.getInt("kernelSize"),
                                                 policy.getInt("kernelSize")))
    
    pexLog.Trace("lsst.ip.diffim.createDefaultPolicy.modifyKernelPolicy", 2,
                 "Footprint grow rad : %d px" % (policy.getInt("fpGrowPix")))

    outStr = ", ".join(["%.2f" % (x) for x in policy.getDoubleArray("alardSigGauss")])
    pexLog.Trace("lsst.ip.diffim.createDefaultPolicy.modifyKernelPolicy", 2,
                 "A/L gaussian sig   : %s px" % (outStr))


