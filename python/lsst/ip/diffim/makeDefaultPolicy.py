import lsst.pex.policy as pexPolicy
import lsst.pex.logging as pexLog
import os
import numpy as num
sigma2fwhm = 2. * num.sqrt(2. * num.log(2.))

POLICY_DEFAULT_FWHM = 3.5

# Here we assume we have the FWHM of the images' PSF in pixels.  We
# will generate an alard-lupton kernel policy based upon these values.
# A complication is that the "sigma" of the AL gaussians is related to
# the FWHM through a scaling of 2.35.

def makeDefaultPolicy(mergePolicy = None, PolicyDictName = "PsfMatchingDictionary.paf"):
    if mergePolicy and isinstance(mergePolicy, basestring):
        policy = pexPolicy.Policy(mergePolicy)
    elif mergePolicy and isinstance(mergePolicy, pexPolicy.Policy):
        policy = mergePolicy
    elif mergePolicy:
        raise RuntimeError, "Can't interpret provided policy: %s" % (policy)
    else:
        policy = pexPolicy.Policy() 
        
    policyFile = pexPolicy.DefaultPolicyFile("ip_diffim", PolicyDictName, "policy")
    defPolicy  = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)
    policy.mergeDefaults(defPolicy.getDictionary())

    return policy

def modifyForImagePsfMatch(defaultPolicy, templateFwhmPix, scienceFwhmPix, minSigma = 0.4):
    psfMatchPolicy = pexPolicy.Policy(
        os.path.join(os.getenv("IP_DIFFIM_DIR"), "policy", "ImagePsfMatchPolicy.paf")
        )
    psfMatchPolicy.mergeDefaults(defaultPolicy)

    # Modify the size of Alard Lupton kernels based upon the images FWHM
    #
    # Note the operation is : template x kernel = science
    #
    # Assuming the template and science image Psfs are Gaussians with
    # the Fwhm above, Fwhm_T **2 + Fwhm_K **2 = Fwhm_S **2
    #
    if templateFwhmPix == scienceFwhmPix:
        # Leave defaults as-is
        pass
    elif float(scienceFwhmPix) / float(templateFwhmPix) > 2.0:
        # Extreme convolution; central Gaussian is at the template
        # scale, outer Gaussian is at the scale to match the two
        # cores, central Gaussian is geometric mean.
        kernelCoreSigma   = templateFwhmPix / sigma2fwhm
        kernelOuterSigma  = num.sqrt(scienceFwhmPix**2 - templateFwhmPix**2) / sigma2fwhm
        kernelMiddleSigma = num.sqrt(kernelCoreSigma * kernelOuterSigma)
        if psfMatchPolicy.get("alardNGauss") == 3:
            psfMatchPolicy.set("alardSigGauss", kernelCoreSigma)
            psfMatchPolicy.add("alardSigGauss", kernelMiddleSigma)
            psfMatchPolicy.add("alardSigGauss", kernelOuterSigma)
        else:
            # Deal with this 
            pass
    elif scienceFwhmPix > templateFwhmPix:
        # Normal convolution; put the bulk of the power in the Gaussian
        # that matches the core Fwhms.  Outer gaussian corresponds to
        # the science image's Fwhm.  Middle is geometric mean to create
        # geometric progression of Gaussian sizes
        kernelCoreFwhm    = num.sqrt(scienceFwhmPix**2 - templateFwhmPix**2)
        kernelCoreSigma   = max(minSigma, kernelCoreFwhm / sigma2fwhm)
        kernelOuterSigma  = scienceFwhmPix / sigma2fwhm
        kernelMiddleSigma = num.sqrt(kernelCoreSigma * kernelOuterSigma)
        if psfMatchPolicy.get("alardNGauss") == 3:
            psfMatchPolicy.set("alardSigGauss", kernelCoreSigma)
            psfMatchPolicy.add("alardSigGauss", kernelMiddleSigma)
            psfMatchPolicy.add("alardSigGauss", kernelOuterSigma)
        else:
            # Deal with this 
            pass
    else:
        # Deconvolution; put the smallest Gaussian at the smallest
        # allowed scale, and define the progression of Gaussians using
        # a method used to derive a deconvolution sum-of-Gaussians
        # from its convolution counterpart.
        #
        # http://iopscience.iop.org/0266-5611/26/8/085002  Equation 40
        psfMatchPolicy    = modifyForDeconvolution(psfMatchPolicy)
        useOuter          = psfMatchPolicy.get('useOuterForDeconv')
        
        kernelDeconvSigma = minSigma
        kernelCoreSigma   = minSigma
        kernelOuterSigma  = templateFwhmPix / sigma2fwhm
        kernelMiddleSigma = 0.5 * (kernelCoreSigma + kernelOuterSigma)
        if psfMatchPolicy.get("alardNGauss") != 4:
            # Deal with this
            pass
        psfMatchPolicy.set("alardSigGauss", kernelDeconvSigma)
        for n in range(3):
            for j in range(n):
                sigma2jn  = (n - j) * kernelMiddleSigma**2
                sigma2jn += j * kernelOuterSigma**2
                sigma2jn -= (n + 1) * kernelCoreSigma**2
                sigmajn   = num.sqrt(sigma2jn)
                psfMatchPolicy.add("alardSigGauss", sigmajn)

        degGauss = psfMatchPolicy.getIntArray("alardDegGauss")
        psfMatchPolicy.set("alardDegGauss", degGauss[0])
        psfMatchPolicy.add("alardDegGauss", degGauss[1])
        psfMatchPolicy.add("alardDegGauss", degGauss[2])
        psfMatchPolicy.add("alardDegGauss", degGauss[3])

        if useOuter:
            psfMatchPolicy.add("alardSigGauss", kernelOuterSigma)
            psfMatchPolicy.add("alardDegGauss", degGauss[3])
            psfMatchPolicy.set("alardNGauss", 5)
            
        
    # Other size estimates
    maxFwhm     = max(scienceFwhmPix, templateFwhmPix)
    kernelSize  = int(psfMatchPolicy.getDouble("kernelSizeFwhmScaling") * maxFwhm)
    kernelSize  = min(kernelSize, psfMatchPolicy.getInt("kernelSizeMax"))
    kernelSize  = max(kernelSize, psfMatchPolicy.getInt("kernelSizeMin"))
    kernelSize += (1 - kernelSize % 2) # make odd sized
    psfMatchPolicy.set("kernelSize", kernelSize)
        
    fpGrow = int(psfMatchPolicy.getPolicy("detectionPolicy").getDouble("fpGrowFwhmScaling") * maxFwhm)
    fpGrow = min(fpGrow, psfMatchPolicy.getPolicy("detectionPolicy").getInt("fpGrowMax"))
    fpGrow = max(fpGrow, psfMatchPolicy.getPolicy("detectionPolicy").getInt("fpGrowMin"))
    psfMatchPolicy.getPolicy("detectionPolicy").set("fpGrowPix", fpGrow)
    
    return psfMatchPolicy

def modifyForDeconvolution(defaultPolicy):
    deconvPolicy = pexPolicy.Policy(
        os.path.join(os.getenv("IP_DIFFIM_DIR"), "policy", "DeconvolutionPolicy.paf")
        )
    deconvPolicy.mergeDefaults(defaultPolicy)
    return deconvPolicy

def modifyForModelPsfMatch(defaultPolicy):
    psfMatchPolicy = pexPolicy.Policy(
        os.path.join(os.getenv("IP_DIFFIM_DIR"), "policy", "ModelPsfMatchPolicy.paf")
        )
    psfMatchPolicy.mergeDefaults(defaultPolicy)
    return psfMatchPolicy

def modifyForSnapSubtraction(defaultPolicy):
    snapPolicy = pexPolicy.Policy(
        os.path.join(os.getenv("IP_DIFFIM_DIR"), "policy", "SnapSubtractionPolicy.paf")
        )
    snapPolicy.mergeDefaults(defaultPolicy)
    return snapPolicy

