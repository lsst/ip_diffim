from builtins import range
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.ip.diffim as ipDiffim
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.log.utils as logUtils
import pylab
import numpy as num

verbosity = 4
logUtils.traceSetAt("ip.diffim", verbosity)

# Define the diagnostic functions on the goodness of fit of the
# spatial model (k_s) to each individual kernel (k)

# import crossCorrelation
# r, d1, d2 = crossCorrelation.makeAutoCorrelation(kernelCellSet, spatialKernel, makePlot=True)
# reload(crossCorrelation)


def makeAutoCorrelation(kernelCellSet, spatialKernel, makePlot=False):
    kImage = afwImage.ImageD(spatialKernel.getDimensions())
    ksImage = afwImage.ImageD(spatialKernel.getDimensions())
    kInfo = []

    candList = []
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(True):  # only look at non-bad candidates
            if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                candList.append(cand.getId())

    r = []
    d1 = []
    d2 = []

    for i in range(len(candList)):
        cand1 = kernelCellSet.getCandidateById(candList[i])
        x1 = cand1.getXCenter()
        y1 = cand1.getYCenter()

        # Original kernel
        kImage1 = afwImage.ImageD(spatialKernel.getDimensions())
        cand1.getKernel(ipDiffim.KernelCandidateF.ORIG).computeImage(kImage1, False)

        # Spatial approximation
        ksImage1 = afwImage.ImageD(spatialKernel.getDimensions())
        spatialKernel.computeImage(ksImage1, False,
                                   afwImage.indexToPosition(int(x1)),
                                   afwImage.indexToPosition(int(y1)))

        # Turn into numarrays for dot product
        ksVector1 = ksImage1.getArray().ravel()
        kVector1 = kImage1.getArray().ravel()

        for j in range(i+1, len(candList)):
            cand2 = kernelCellSet.getCandidateById(candList[j])
            x2 = cand2.getXCenter()
            y2 = cand2.getYCenter()

            # Original kernel
            kImage2 = afwImage.ImageD(spatialKernel.getDimensions())
            cand2.getKernel(ipDiffim.KernelCandidateF.ORIG).computeImage(kImage2, False)

            # Spatial approximation
            ksImage2 = afwImage.ImageD(spatialKernel.getDimensions())
            spatialKernel.computeImage(ksImage2, False,
                                       afwImage.indexToPosition(int(x2)),
                                       afwImage.indexToPosition(int(y2)))

            # Turn into numarrays for dot product
            ksVector2 = ksImage2.getArray().ravel()
            kVector2 = kImage2.getArray().ravel()

            ###

            # Add to cross correlation functions
            r.append(num.sqrt((x1 - x2)**2 + (y1 - y2)**2))
            d1.append(correlationD1(kVector1, ksVector1, kVector2, ksVector2))
            d2.append(correlationD2(kVector1, ksVector1, kVector2, ksVector2))

    r = num.array(r)
    d1 = num.array(d1)
    d2 = num.array(d2)

    if makePlot:
        plotCrossCorrelation(r, d1, d2)

    return r, d1, d2


def plotCrossCorrelation(r, d1, d2, bsize=100):
    import pylab

    dbin = num.arange(0, max(r) + bsize, bsize)
    xbin = []
    d1av = []
    d1rms = []
    d2av = []
    d2rms = []

    for i in range(len(dbin)-1):
        idx, = num.where((r > dbin[i]) & (r <= dbin[i+1]))
        if len(idx) >= 1:
            d1av.append(num.mean(d1[idx]))
            d1rms.append(num.std(d1[idx]))
            d2av.append(num.mean(d2[idx]))
            d2rms.append(num.std(d2[idx]))
            xbin.append(num.mean(r[idx]))

    xbin = num.array(xbin)
    d1av = num.array(d1av)
    d1rms = num.array(d1rms)
    d2av = num.array(d2av)
    d2rms = num.array(d2rms)

    pylab.figure()

    pylab.plot(r, d1, 'ro', ms=2, alpha=0.25)
    pylab.plot(r, d2, 'gs', ms=2, alpha=0.25)
    pylab.errorbar(xbin, d1av, yerr=d1rms, fmt='ro', label='D1')
    pylab.errorbar(xbin, d2av, yerr=d2rms, fmt='gs', label='D2')
    pylab.xlabel('r')
    pylab.ylabel('Correlation Function')
    pylab.legend(loc=1, numpoints=1)
    pylab.axhline(y=0, color='k', linestyle='--')

    pylab.show()


def correlationD1(k1, ks1, k2, ks2):
    # Autocorrelation in residuals
    # [ (k - k_s) * (k - k_s) ] (r)
    t1 = (k1 - ks1)
    t2 = (k2 - ks2)

    return num.dot(t1, t2)


def correlationD2(k1, ks1, k2, ks2):
    # Data-residual cross correlation
    # [ k * (k - k_s) + (k - k_s) * k ] (r)
    t11 = k1
    t12 = (k2 - ks2)

    t21 = (k1 - ks1)
    t22 = k2

    return num.dot(t11, t12) + num.dot(t21, t22)

#################


def addNoise(mi):
    img = mi.getImage()
    seed = int(afwMath.makeStatistics(mi.getVariance(), afwMath.MAX).getValue())+1
    rdm = afwMath.Random(afwMath.Random.MT19937, seed)
    rdmImage = img.Factory(img.getDimensions())
    afwMath.randomGaussianImage(rdmImage, rdm)
    rdmImage *= num.sqrt(seed)
    img += rdmImage

    # and don't forget to add to the variance
    var = mi.getVariance()
    var += 1.0


def testAutoCorrelation(orderMake, orderFit, inMi=None, display=False):
    config = ipDiffim.ImagePsfMatchTask.ConfigClass()
    config.kernel.name = "AL"
    subconfig = config.kernel.active

    subconfig.fitForBackground = True

    stride = 100

    if inMi == None:
        width = 512
        height = 2048
        inMi = afwImage.MaskedImageF(afwGeom.Extent2I(width, height))
        for j in num.arange(stride//2, height, stride):
            j = int(j)
            for i in num.arange(stride//2, width, stride):
                i = int(i)
                inMi.set(i-1, j-1, (100., 0x0, 1.))
                inMi.set(i-1, j+0, (100., 0x0, 1.))
                inMi.set(i-1, j+1, (100., 0x0, 1.))
                inMi.set(i+0, j-1, (100., 0x0, 1.))
                inMi.set(i+0, j+0, (100., 0x0, 1.))
                inMi.set(i+0, j+1, (100., 0x0, 1.))
                inMi.set(i+1, j-1, (100., 0x0, 1.))
                inMi.set(i+1, j+0, (100., 0x0, 1.))
                inMi.set(i+1, j+1, (100., 0x0, 1.))

    addNoise(inMi)

    kSize = subconfig.kernelSize

    basicGaussian1 = afwMath.GaussianFunction2D(2., 2., 0.)
    basicKernel1 = afwMath.AnalyticKernel(kSize, kSize, basicGaussian1)

    basicGaussian2 = afwMath.GaussianFunction2D(5., 3., 0.5 * num.pi)
    basicKernel2 = afwMath.AnalyticKernel(kSize, kSize, basicGaussian2)

    basisList = []
    basisList.append(basicKernel1)
    basisList.append(basicKernel2)

    spatialKernelFunction = afwMath.PolynomialFunction2D(orderMake)
    spatialKernel = afwMath.LinearCombinationKernel(basisList, spatialKernelFunction)
    kCoeffs = [[1.0 for x in range(1, spatialKernelFunction.getNParameters()+1)],
               [0.01 * x for x in range(1, spatialKernelFunction.getNParameters()+1)]]
    spatialKernel.setSpatialParameters(kCoeffs)

    cMi = afwImage.MaskedImageF(inMi.getDimensions())
    afwMath.convolve(cMi, inMi, spatialKernel, True)

    if display:
        ds9.mtv(inMi.getImage(), frame=1)
        ds9.mtv(inMi.getVariance(), frame=2)

        ds9.mtv(cMi.getImage(), frame=3)
        ds9.mtv(cMi.getVariance(), frame=4)

    subconfig.spatialKernelOrder = orderFit
    subconfig.sizeCellX = stride
    subconfig.sizeCellY = stride
    psfmatch = ipDiffim.ImagePsfMatchTask(config=config)
    result = psfmatch.run(inMi, cMi, "subtractMaskedImages")

    differenceMaskedImage = result.subtractedImage
    spatialKernel = result.psfMatchingKernel
    spatialBg = result.backgroundModel
    kernelCellSet = result.kernelCellSet
    makeAutoCorrelation(kernelCellSet, spatialKernel, makePlot=True)


def doOverConstrained(inMi=None, display=False):
    testAutoCorrelation(orderMake=1, orderFit=3, inMi=inMi, display=display)


def doUnderConstrained(inMi=None, display=False):
    testAutoCorrelation(orderMake=2, orderFit=0, inMi=inMi, display=display)

if __name__ == '__main__':
    doOverConstrained(display=True)
    doUnderConstrained(display=True)
