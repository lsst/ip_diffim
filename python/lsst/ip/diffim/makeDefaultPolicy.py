import lsst.pex.policy as pexPolicy
import lsst.pex.logging as pexLog

import math
sigma2fwhm = 2. * math.sqrt(2. * math.log(2.))

POLICY_DEFAULT_FWHM = 3.5

# Here we assume we have the FWHM of the images' PSF in pixels.  We
# will generate an alard-lupton kernel policy based upon these values.
# A complication is that the "sigma" of the AL gaussians is related to
# the FWHM through a scaling of 2.35.

def makeDefaultPolicy(PolicyDictName = "PsfMatchingDictionary.paf", mergePolicyPath = None,
                      fwhm = POLICY_DEFAULT_FWHM, modify = False): 
    if mergePolicyPath:
        policy = pexPolicy.Policy(mergePolicyPath)
    else:
        policy = pexPolicy.Policy() 
        
    policyFile = pexPolicy.DefaultPolicyFile("ip_diffim", PolicyDictName, "policy")
    defPolicy  = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)
    policy.mergeDefaults(defPolicy.getDictionary())

    if modify or policy.get("scaleByFwhm"):
        modifyKernelPolicy(policy, fwhm)
    return policy

def modifyKernelPolicy(policy, fwhm = POLICY_DEFAULT_FWHM):
    # Modify the kernel policy parameters based upon the images FWHM

    kernelSize  = int(policy.getDouble("kernelSizeFwhmScaling") * fwhm)
    kernelSize  = min(kernelSize, policy.getInt("kernelSizeMax"))
    kernelSize  = max(kernelSize, policy.getInt("kernelSizeMin"))
    kernelSize += (1 - kernelSize % 2) # make odd sized
    policy.set("kernelSize", kernelSize)

    fpGrow = int(policy.getPolicy("detectionPolicy").getDouble("fpGrowFwhmScaling") * fwhm)
    fpGrow = min(fpGrow, policy.getPolicy("detectionPolicy").getInt("fpGrowMax"))
    fpGrow = max(fpGrow, policy.getPolicy("detectionPolicy").getInt("fpGrowMin"))
    policy.getPolicy("detectionPolicy").set("fpGrowPix", fpGrow)

    sigma = fwhm / sigma2fwhm
    alardSig = [x*sigma for x in policy.getDoubleArray("alardSigFwhmScaling")]
    policy.set("alardSigGauss", alardSig[0])
    for sig in alardSig[1:]:
        policy.add("alardSigGauss", sig)
    
    pexLog.Trace("lsst.ip.diffim.makeDefaultPolicy.modifyKernelPolicy", 2,
                 "Using Psf Fwhm     : %.2f px" % (fwhm))

    pexLog.Trace("lsst.ip.diffim.makeDefaultPolicy.modifyKernelPolicy", 2,
                 "Kernel size   : %d x %d px" % (policy.getInt("kernelSize"),
                                                 policy.getInt("kernelSize")))
    
    pexLog.Trace("lsst.ip.diffim.makeDefaultPolicy.modifyKernelPolicy", 2,
                 "Footprint grow rad : %d px" % (policy.getPolicy("detectionPolicy").getInt("fpGrowPix")))

    outStr = ", ".join(["%.2f" % (x) for x in policy.getDoubleArray("alardSigGauss")])
    pexLog.Trace("lsst.ip.diffim.makeDefaultPolicy.modifyKernelPolicy", 2,
                 "A/L gaussian sig   : %s px" % (outStr))



def modifyForDeconvolution(policy):

    # Modify AL shapes
    degGauss = policy.getIntArray("alardDegGaussDeconvolution")
    policy.set("alardDegGauss", degGauss[0])
    for deg in degGauss[1:]:
        policy.add("alardDegGauss", deg)
        
    sigGauss = policy.getDoubleArray("alardSigGaussDeconvolution")
    policy.set("alardSigGauss", sigGauss[0])
    for sig in sigGauss[1:]:
        policy.add("alardSigGauss", sig)

    outStr = ", ".join(["%.2f" % (x) for x in policy.getDoubleArray("alardSigGauss")])
    pexLog.Trace("lsst.ip.diffim.makeDefaultPolicy.modifyForDeconvolution", 2,
                 "A/L gaussian sig   : %s px" % (outStr))
    outStr = ", ".join(["%.2f" % (x) for x in policy.getIntArray("alardDegGauss")])
    pexLog.Trace("lsst.ip.diffim.makeDefaultPolicy.modifyForDeconvolution", 2,
                 "A/L gaussian deg   : %s px" % (outStr))

    # Don't use core stats (core has large resids)
    policy.set("useCoreStats", False)

