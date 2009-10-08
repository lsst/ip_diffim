import os
import eups
import time
import lsst.afw.image.imageLib as afwImage
import lsst.afw.math.mathLib as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.policy as pexPolicy
import lsst.pex.logging as pexLogging

import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils

display = True

Verbosity = 5
pexLogging.Trace_setVerbosity("lsst.ip.diffim", Verbosity)

defDataDir   = eups.productDir("afwdata") 
imageProcDir = eups.productDir("ip_diffim")

defSciencePath  = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_1")
defTemplatePath = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_1_tmpl")
defPolicyPath   = os.path.join(imageProcDir, "pipeline", "ImageSubtractStageDictionary.paf")
defOutputPath   = "diffImage"

templateMaskedImage = afwImage.MaskedImageF(defTemplatePath)
scienceMaskedImage  = afwImage.MaskedImageF(defSciencePath)
policy              = pexPolicy.Policy.createPolicy(defPolicyPath)

footprints = ipDiffim.getCollectionOfFootprintsForPsfMatching(templateMaskedImage,
                                                              scienceMaskedImage,
                                                              policy)

kernelCellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(templateMaskedImage.getX0(),
                                                                     templateMaskedImage.getY0()),
                                                     templateMaskedImage.getWidth(),
                                                     templateMaskedImage.getHeight()),
                                       policy.getInt("sizeCellX"),
                                       policy.getInt("sizeCellY"))

for fp in footprints:
    bbox = fp.getBBox()
    xC   = 0.5 * ( bbox.getX0() + bbox.getX1() )
    yC   = 0.5 * ( bbox.getY0() + bbox.getY1() )
    tmi  = afwImage.MaskedImageF(templateMaskedImage,  bbox)
    smi  = afwImage.MaskedImageF(scienceMaskedImage, bbox)
    
    #if not goodKernelCandidate():
    #    continue
    
    cand = ipDiffim.makeKernelCandidate(xC, yC, tmi, smi)
    kernelCellSet.insertCandidate(cand)

# Do we reduce the dimensionality first by doing PCA on the DF/DFr kernels?
doPca = False
if doPca:
    nEigenComponents            = policy.getInt("nEigenComponents")
    spatialOrder                = policy.getInt("spatialOrder")
    nStarPerCell                = policy.getInt("nStarPerCell")
    nStarPerCellSpatialFit      = policy.getInt("nStarPerCellSpatialFit")
    tolerance                   = policy.getDouble("tolerance")
    reducedChi2ForPsfCandidates = policy.getDouble("reducedChi2ForPsfCandidates")
    nIterForPsf                 = policy.getInt("nIterForPsf")
        
kFunctor       = ipDiffim.createKernelFunctor(policy)
nIterForKernel = policy.getInt("maxSpatialIterations")

if doPca:
    kernel, eigenValues = ipDiffim.createPcaBasisFromCandidates(kernelCellSet, nEigenComponents,
                                                                spatialOrder, nStarPerCell)
else:
    spatialKernel, spatialBg = ipDiffim.fitSpatialKernelFromCandidates(kFunctor, kernelCellSet, policy)

nCellAll  = 0
nCellUsed = 0
nCandAll  = 0
nCandUsed = 0
for cell in kernelCellSet.getCellList():
    cellFound = False
    for cand in cell.begin(False): # False = include bad candidates
        cand = ipDiffim.cast_KernelCandidateF(cand)
        if cand.getStatus() == afwMath.SpatialCellCandidate.BAD:
            nCandAll  += 1
        elif cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
            nCandUsed += 1
            cellFound  = True
        #else:
        #    print 'Uh', cand.getStatus()
    if cellFound:
        nCellUsed += 1
    nCellAll += 1

print 'Using %d / %d Cells; %d / %d Candiates' % (nCellUsed, nCellAll, nCandUsed, nCandAll)

if nCellUsed == 0 or nCandUsed == 0:
    print 'Damn, unable to do a thing'

# Lets see what we got
if display:
    mos = displayUtils.Mosaic()

    # Inputs
    frame = 0
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(False): # False = include bad candidates
            cand  = ipDiffim.cast_KernelCandidateF(cand)
            rchi2 = cand.getChi2()
                
            try:
                im = cand.getImage()
                mos.append(im, "#%d: %.1f (%s)" % (cand.getId(), rchi2, cand.getStatus()))
            except Exception, e:
                pass
    mosaic = mos.makeMosaic()
    ds9.mtv(mosaic, frame=frame)
    mos.drawLabels(frame=frame)

    # Bases
    frame = 1
    mos.reset()
    basisList = spatialKernel.getKernelList()
    for idx in range(len(basisList)):
        kernel = basisList[idx]
        im   = afwImage.ImageD(spatialKernel.getDimensions())
        ksum = kernel.computeImage(im, False)
        mos.append(im, "K%d" % (idx))
    mosaic = mos.makeMosaic()
    ds9.mtv(mosaic, frame=frame)
    mos.drawLabels(frame=frame)
        

    # Spatial model
    frame = 2
    mos.reset()
    width = templateMaskedImage.getWidth()
    height = templateMaskedImage.getHeight()
    stamps = []; stampInfo = []
    for x in (0, width//2, width):
        for y in (0, height//2, height):
            im   = afwImage.ImageD(spatialKernel.getDimensions())
            ksum = spatialKernel.computeImage(im, False, float(x), float(y))
            mos.append(im, "x=%d y=%d kSum=%.2f" % (x, y, ksum))

    mosaic = mos.makeMosaic()
    ds9.mtv(mosaic, frame=frame)
    mos.drawLabels(frame=frame)
            

    # Background
    frame = 3
    backgroundIm = afwImage.ImageF(templateMaskedImage.getDimensions())
    ipDiffim.addSomethingToImage(backgroundIm, spatialBg)
    ds9.mtv(backgroundIm, frame=frame)

    # Diffim!
    frame = 4
    diffIm = ipDiffim.convolveAndSubtract(templateMaskedImage,
                                          scienceMaskedImage,
                                          spatialKernel,
                                          spatialBg)
    ds9.mtv(diffIm, frame=frame)


