import numpy, pylab
import lsst.ip.diffim as ipDiffim
from lsst.ip.diffim.diffimDebug import * 
from lsst.pex.logging import Trace


# Helper routines
def scaleImageTo8Bit(indata):
    scaled  = indata - indata.min()  # lowest value is 0.0
    scaled /= scaled.max()           # highest value is 1.0
    return numpy.minimum(numpy.floor(numpy.maximum(scaled * 256., 0)), 255)

def scaleImageForDisplay_nonlinear(indata, nonlinearity=2.0):
    return scaleImageTo8Bit( numpy.arcsinh(indata*nonlinearity) / nonlinearity )

def scaleImageForDisplay_log10(indata):
    return numpy.log10(indata)


# Plotting routines
def plotBackground( (backgroundFunction, cols, rows), nbins=100, outfile=None):
    print backgroundFunction.toString()
    
    pylab.figure()
    background = numpy.zeros((nbins,nbins))
    for yi in range(nbins):
        yidx = -1 + yi * 0.5 / nbins

        for xi in range(nbins):
            xidx = -1 + xi * 0.5 / nbins

            background[yi][xi] = backgroundFunction(xidx, yidx)

    sp_bg = pylab.subplot(111)
    im = sp_bg.imshow(background,
                      cmap=pylab.cm.jet,
                      origin='lower',
                      aspect=2.0,
                      extent=[-1, 1, -1, 1])
    sp_bg.plot(cols, rows, 'kx') # actual constraints
    sp_bg.axhline(y=0, c='k', linestyle='--')
    sp_bg.axvline(x=0, c='k', linestyle='--')
    pylab.colorbar(im)

    if outfile != None:
        pylab.savefig(outfile)

def makeApproxKernelPlot(difi, kernelInfo, title):
    meanKernel, eigenKernels = kernelInfo

    meanImage = imageToMatrix( meanKernel.computeNewImage(False)[0] )
    eigenImages = []
    for i in range(len(eigenKernels)):
        eigenImages.append( imageToMatrix( eigenKernels[i].computeNewImage(False)[0] ) )

    originalImage = imageToMatrix( difi.getSingleKernelPtr().computeNewImage(False)[0] )
    
    # Now make PCA-approximations of each Kernel, and watch how the
    # diffim residuals decrease as a function of # of eigenKernels.
    kResids = []
    dResids = []

    # Calculate statistics using approximated Kernel
    kModel = meanImage.copy()

    approxKernelPtr = matrixToKernelPtr(kModel)
    kResid = numpy.sum( numpy.abs( kModel - originalImage )**2 )
    kResids.append( (kModel, kResid) )
    
    approxDiffIm = ipDiffim.convolveAndSubtract(difi.getImageToConvolvePtr().get(),
                                                difi.getImageToNotConvolvePtr().get(),
                                                approxKernelPtr,
                                                difi.getSingleBackground())
    approxStatistics = ipDiffim.DifferenceImageStatisticsF(approxDiffIm)
    Trace('lsst.ip.diffim', 5,
          '%s Mean Kernel : Kernel Sum = %.2f, Diffim residuals = %.2f +/- %.2f sigma' % (
        title,
        approxKernelPtr.computeNewImage(False)[1],
        approxStatistics.getResidualMean(),
        approxStatistics.getResidualStd()
        ))
    dResids.append( ( imageToMatrix(approxDiffIm.getImage()),
                      approxStatistics.getResidualMean(),
                      approxStatistics.getResidualStd() )
                    )

    #import pdb
    #pdb.set_trace()
    
    for nKernel in range(len(eigenKernels)):
        eCoeff  = numpy.dot( numpy.ravel( (originalImage - meanImage) ).T, numpy.ravel( eigenImages[nKernel] ) )
        kModel += eCoeff * eigenImages[nKernel]

        approxKernelPtr = matrixToKernelPtr(kModel)
        kResid = numpy.sum( numpy.abs( kModel - originalImage )**2 )
        kResids.append( (kModel, kResid) )
        
        approxDiffIm = ipDiffim.convolveAndSubtract(difi.getImageToConvolvePtr().get(),
                                                    difi.getImageToNotConvolvePtr().get(),
                                                    approxKernelPtr,
                                                    difi.getSingleBackground())
        approxStatistics = ipDiffim.DifferenceImageStatisticsF(approxDiffIm)
        Trace('lsst.ip.diffim', 5,
              '%s PCA Kernel %d : Kernel Sum = %.2f, Diffim residuals = %.2f +/- %.2f sigma' % (
            title, nKernel,
            approxKernelPtr.computeNewImage(False)[1],
            approxStatistics.getResidualMean(),
            approxStatistics.getResidualStd()
            ))
        dResids.append( ( imageToMatrix(approxDiffIm.getImage()),
                          approxStatistics.getResidualMean(),
                          approxStatistics.getResidualStd() )
                        )

    return kResids, dResids
    

def approxKernelPlot(difi, kernelInfo, title, fontsize=10, outfile=None):
    pylab.figure()

    kResids, dResids = makeApproxKernelPlot(difi, kernelInfo, title)

    di0   = [0.050, 0.725, 0.150, 0.150]
    di1   = [0.200, 0.725, 0.150, 0.150]
    di2   = [0.350, 0.725, 0.150, 0.150]
    di3   = [0.050, 0.575, 0.150, 0.150]
    di4   = [0.200, 0.575, 0.150, 0.150]
    di5   = [0.350, 0.575, 0.150, 0.150]
    dspec = [0.575, 0.550, 0.375, 0.350]

    ki0  = [0.050, 0.275, 0.150, 0.150]
    ki1  = [0.200, 0.275, 0.150, 0.150]
    ki2  = [0.350, 0.275, 0.150, 0.150]
    ki3  = [0.050, 0.125, 0.150, 0.150]
    ki4  = [0.200, 0.125, 0.150, 0.150]
    ki5  = [0.350, 0.125, 0.150, 0.150]
    kspec= [0.575, 0.100, 0.375, 0.350]
    
    #
    ########
    #

    sp_di0 = pylab.axes(di0)
    sp_di0.imshow( scaleImageForDisplay_nonlinear(dResids[0][0]),
                   cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_di0.set_title('DI1', fontsize=fontsize, weight='bold')

    sp_di1 = pylab.axes(di1)
    if len(dResids) > 1:
        sp_di1.imshow( scaleImageForDisplay_nonlinear(dResids[1][0]),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_di1.set_title('DI2', fontsize=fontsize, weight='bold')

    sp_di2 = pylab.axes(di2)
    if len(dResids) > 2:
        sp_di2.imshow( scaleImageForDisplay_nonlinear(dResids[2][0]),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_di2.set_title('DI3', fontsize=fontsize, weight='bold')

    sp_di3 = pylab.axes(di3)
    if len(dResids) > 3:
        sp_di3.imshow( scaleImageForDisplay_nonlinear(dResids[3][0]),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_di3.set_xlabel('DI4', fontsize=fontsize, weight='bold')

    sp_di4 = pylab.axes(di4)
    if len(dResids) > 4:
        sp_di4.imshow( scaleImageForDisplay_nonlinear(dResids[4][0]),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_di4.set_xlabel('DI5', fontsize=fontsize, weight='bold')

    sp_di5 = pylab.axes(di5)
    if len(dResids) > 5:
        sp_di5.imshow( scaleImageForDisplay_nonlinear(dResids[5][0]),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_di5.set_xlabel('DI6', fontsize=fontsize, weight='bold')

    sp_dspec = pylab.axes(dspec)
    means = [x[1] for x in dResids]
    stds  = [x[2] for x in dResids]
    sp_dspec.plot(range(len(means)), means, linestyle='-', c='b')
    sp_dspec.plot(range(len(stds)), stds, linestyle='--', c='r')
    sp_dspec.set_xlabel('N', fontsize=fontsize+1, weight='bold')
    sp_dspec.set_ylabel('Residuals', fontsize=fontsize+1, weight='bold')
    sp_dspec.set_title(title, fontsize=fontsize+1, weight='bold')
    sp_dspec.axhline(y=means[-1], c='b', linestyle=':')
    sp_dspec.axhline(y=stds[-1],  c='r', linestyle=':')

    pylab.setp(sp_di0.get_xticklabels()+sp_di0.get_yticklabels()+
               sp_di1.get_xticklabels()+sp_di1.get_yticklabels()+
               sp_di2.get_xticklabels()+sp_di2.get_yticklabels()+
               sp_di3.get_xticklabels()+sp_di3.get_yticklabels()+
               sp_di4.get_xticklabels()+sp_di4.get_yticklabels()+
               sp_di5.get_xticklabels()+sp_di5.get_yticklabels(), visible=False)
    pylab.setp(sp_dspec.get_xticklabels()+sp_dspec.get_yticklabels(), fontsize=fontsize)

    #
    ########
    #

    sp_ki0 = pylab.axes(ki0)
    sp_ki0.imshow( scaleImageForDisplay_nonlinear(kResids[0][0]),
                   cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ki0.set_title('Mean', fontsize=fontsize, weight='bold')

    sp_ki1 = pylab.axes(ki1)
    if len(kResids) > 1:
        sp_ki1.imshow( scaleImageForDisplay_nonlinear(kResids[1][0]),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ki1.set_title('PK1', fontsize=fontsize, weight='bold')

    sp_ki2 = pylab.axes(ki2)
    if len(kResids) > 2:
        sp_ki2.imshow( scaleImageForDisplay_nonlinear(kResids[2][0]),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ki2.set_title('PK2', fontsize=fontsize, weight='bold')

    sp_ki3 = pylab.axes(ki3)
    if len(kResids) > 3:
        sp_ki3.imshow( scaleImageForDisplay_nonlinear(kResids[3][0]),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ki3.set_xlabel('PK3', fontsize=fontsize, weight='bold')

    sp_ki4 = pylab.axes(ki4)
    if len(kResids) > 4:
        sp_ki4.imshow( scaleImageForDisplay_nonlinear(kResids[4][0]),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ki4.set_xlabel('PK4', fontsize=fontsize, weight='bold')

    sp_ki5 = pylab.axes(ki5)
    if len(kResids) > 5:
        sp_ki5.imshow( scaleImageForDisplay_nonlinear(kResids[5][0]),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ki5.set_xlabel('PK5', fontsize=fontsize, weight='bold')

    sp_kspec = pylab.axes(kspec)
    msresids = [x[1] for x in kResids]
    sp_kspec.plot(range(len(msresids)), msresids, linestyle='--', c='k')
    sp_kspec.set_xlabel('N', fontsize=fontsize+1, weight='bold')
    sp_kspec.set_ylabel('M.S.E.', fontsize=fontsize+1, weight='bold')
    sp_kspec.semilogy()

    pylab.setp(sp_ki0.get_xticklabels()+sp_ki0.get_yticklabels()+
               sp_ki1.get_xticklabels()+sp_ki1.get_yticklabels()+
               sp_ki2.get_xticklabels()+sp_ki2.get_yticklabels()+
               sp_ki3.get_xticklabels()+sp_ki3.get_yticklabels()+
               sp_ki4.get_xticklabels()+sp_ki4.get_yticklabels()+
               sp_ki5.get_xticklabels()+sp_ki5.get_yticklabels(), visible=False)
    pylab.setp(sp_kspec.get_xticklabels()+sp_kspec.get_yticklabels(), fontsize=fontsize)

    if outfile != None:
        pylab.savefig(outfile)
        

def eigenKernelPlot(convInfo, deconvInfo, fontsize=10, outfile=None):
    pylab.figure()
    
    convMKernel,   convEKernels,   convEVals   = convInfo
    deconvMKernel, deconvEKernels, deconvEVals = deconvInfo
    
    ckm   = [0.050, 0.725, 0.150, 0.150]
    ck1   = [0.200, 0.725, 0.150, 0.150]
    ck2   = [0.350, 0.725, 0.150, 0.150]
    ck3   = [0.050, 0.575, 0.150, 0.150]
    ck4   = [0.200, 0.575, 0.150, 0.150]
    ck5   = [0.350, 0.575, 0.150, 0.150]
    cspec = [0.575, 0.550, 0.375, 0.350]

    dckm  = [0.050, 0.275, 0.150, 0.150]
    dck1  = [0.200, 0.275, 0.150, 0.150]
    dck2  = [0.350, 0.275, 0.150, 0.150]
    dck3  = [0.050, 0.125, 0.150, 0.150]
    dck4  = [0.200, 0.125, 0.150, 0.150]
    dck5  = [0.350, 0.125, 0.150, 0.150]
    dcspec= [0.575, 0.100, 0.375, 0.350]
    
    #
    ########
    #

    sp_ckm = pylab.axes(ckm)
    sp_ckm.imshow( scaleImageForDisplay_nonlinear( imageToMatrix( convMKernel.computeNewImage(False)[0] ) ), 
                   cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ckm.set_title('Mean', fontsize=fontsize, weight='bold')

    sp_ck1 = pylab.axes(ck1)
    if len(convEKernels) > 0:
        sp_ck1.imshow( scaleImageForDisplay_nonlinear( imageToMatrix( convEKernels[0].computeNewImage(False)[0] ) ),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ck1.set_title('PK1', fontsize=fontsize, weight='bold')

    sp_ck2 = pylab.axes(ck2)
    if len(convEKernels) > 1:
        sp_ck2.imshow( scaleImageForDisplay_nonlinear( imageToMatrix( convEKernels[1].computeNewImage(False)[0] ) ),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ck2.set_title('PK2', fontsize=fontsize, weight='bold')

    sp_ck3 = pylab.axes(ck3)
    if len(convEKernels) > 2:
        sp_ck3.imshow( scaleImageForDisplay_nonlinear( imageToMatrix( convEKernels[2].computeNewImage(False)[0] ) ),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ck3.set_xlabel('PK3', fontsize=fontsize, weight='bold')

    sp_ck4 = pylab.axes(ck4)
    if len(convEKernels) > 3:
        sp_ck4.imshow( scaleImageForDisplay_nonlinear( imageToMatrix( convEKernels[3].computeNewImage(False)[0] ) ),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ck4.set_xlabel('PK4', fontsize=fontsize, weight='bold')

    sp_ck5 = pylab.axes(ck5)
    if len(convEKernels) > 4:
        sp_ck5.imshow( scaleImageForDisplay_nonlinear( imageToMatrix( convEKernels[4].computeNewImage(False)[0] ) ),
                       cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ck5.set_xlabel('PK5', fontsize=fontsize, weight='bold')

    sp_cspec  = pylab.axes(cspec)
    sp_chist  = numpy.cumsum( convEVals )
    sp_chist /= sp_chist[-1]
    sp_cspec.plot([x+1 for x in range(len(sp_chist))], sp_chist, linestyle='--')
    sp_cspec.set_xlabel('N', fontsize=fontsize+1, weight='bold')
    sp_cspec.set_ylabel('Cumulative', fontsize=fontsize+1, weight='bold')
    sp_cspec.set_xlim( (1, len(sp_chist)) )
    sp_cspec.set_ylim( (0, 1) )

    pylab.setp(sp_ckm.get_xticklabels()+sp_ckm.get_yticklabels()+
               sp_ck1.get_xticklabels()+sp_ck1.get_yticklabels()+
               sp_ck2.get_xticklabels()+sp_ck2.get_yticklabels()+
               sp_ck3.get_xticklabels()+sp_ck3.get_yticklabels()+
               sp_ck4.get_xticklabels()+sp_ck4.get_yticklabels()+
               sp_ck5.get_xticklabels()+sp_ck5.get_yticklabels(), visible=False)
    pylab.setp(sp_cspec.get_xticklabels()+sp_cspec.get_yticklabels(), fontsize=fontsize)

    #
    ########
    #

    sp_dckm = pylab.axes(dckm)
    sp_dckm.imshow( scaleImageForDisplay_nonlinear( imageToMatrix( deconvMKernel.computeNewImage(False)[0] ) ),
                    cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_dckm.set_title('Mean', fontsize=fontsize, weight='bold')

    sp_dck1 = pylab.axes(dck1)
    if len(deconvEKernels) > 0:
        sp_dck1.imshow( scaleImageForDisplay_nonlinear( imageToMatrix( deconvEKernels[0].computeNewImage(False)[0] ) ),
                        cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_dck1.set_title('PK1', fontsize=fontsize, weight='bold')

    sp_dck2 = pylab.axes(dck2)
    if len(deconvEKernels) > 1:
        sp_dck2.imshow( scaleImageForDisplay_nonlinear( imageToMatrix( deconvEKernels[1].computeNewImage(False)[0] ) ),
                        cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_dck2.set_title('PK2', fontsize=fontsize, weight='bold')

    sp_dck3 = pylab.axes(dck3)
    if len(deconvEKernels) > 2:
        sp_dck3.imshow( scaleImageForDisplay_nonlinear( imageToMatrix( deconvEKernels[2].computeNewImage(False)[0] ) ),
                        cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_dck3.set_xlabel('PK3', fontsize=fontsize, weight='bold')

    sp_dck4 = pylab.axes(dck4)
    if len(deconvEKernels) > 3:
        sp_dck4.imshow( scaleImageForDisplay_nonlinear( imageToMatrix( deconvEKernels[3].computeNewImage(False)[0] ) ),
                        cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_dck4.set_xlabel('PK4', fontsize=fontsize, weight='bold')

    sp_dck5 = pylab.axes(dck5)
    if len(deconvEKernels) > 4:
        sp_dck5.imshow( scaleImageForDisplay_nonlinear( imageToMatrix( deconvEKernels[4].computeNewImage(False)[0] ) ),
                        cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_dck5.set_xlabel('PK5', fontsize=fontsize, weight='bold')

    sp_dcspec  = pylab.axes(dcspec)
    sp_dchist  = numpy.cumsum( deconvEVals )
    sp_dchist /= sp_dchist[-1]
    sp_dcspec.plot([x+1 for x in range(len(sp_dchist))], sp_dchist, linestyle='--')
    sp_dcspec.set_xlabel('N', fontsize=fontsize+1, weight='bold')
    sp_dcspec.set_ylabel('Cumulative', fontsize=fontsize+1, weight='bold')
    sp_dcspec.set_xlim( (1, len(sp_dchist)) )
    sp_dcspec.set_ylim( (0, 1) )

    pylab.setp(sp_dckm.get_xticklabels()+sp_dckm.get_yticklabels()+
               sp_dck1.get_xticklabels()+sp_dck1.get_yticklabels()+
               sp_dck2.get_xticklabels()+sp_dck2.get_yticklabels()+
               sp_dck3.get_xticklabels()+sp_dck3.get_yticklabels()+
               sp_dck4.get_xticklabels()+sp_dck4.get_yticklabels()+
               sp_dck5.get_xticklabels()+sp_dck5.get_yticklabels(), visible=False)
    pylab.setp(sp_dcspec.get_xticklabels()+sp_dcspec.get_yticklabels(), fontsize=fontsize)

    if outfile != None:
        pylab.savefig(outfile)


def makeSigmaHistograms(info, title, fontsize, windows):
    
    bins   = numpy.arange(-5.0, +5.0, 0.5)
    theory = []
    for b in bins:
        theory.append( numpy.exp( -0.5 * b**2 ) / numpy.sqrt(2. * numpy.pi) )

    sp_im0 = pylab.axes(windows[0])
    sp_im0.imshow( scaleImageForDisplay_nonlinear(info[0]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
    sp_im0.set_title('Convolved', fontsize=fontsize, weight='bold')
    
    sp_im1 = pylab.axes(windows[1])
    sp_im1.imshow( scaleImageForDisplay_nonlinear(info[1]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
    sp_im1.set_title('Science', fontsize=fontsize, weight='bold')

    sp_im2 = pylab.axes(windows[2])
    sp_im2.imshow( scaleImageForDisplay_nonlinear(info[2]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
    sp_im2.set_title('Difference', fontsize=fontsize, weight='bold')

    sp_ker = pylab.axes(windows[3])
    sp_ker.imshow( scaleImageForDisplay_nonlinear(info[3]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
    sp_ker.set_title('Kernel', fontsize=fontsize, weight='bold')
    
    sp_his = pylab.axes(windows[4])
    sp_his.hist(info[4], bins=bins, normed=True)
    sp_his.plot(bins, theory, 'r-')
    sp_his.set_xlabel('Sigma', fontsize=fontsize, weight='bold')
    sp_his.set_ylabel('Fraction', fontsize=fontsize, weight='bold')
    sp_his.set_title('%s : %.3f +/- %.3f' % (title, info[4].mean(), info[4].std()), fontsize=fontsize+1, weight='bold')

    pylab.setp(sp_im0.get_xticklabels()+sp_im0.get_yticklabels()+
               sp_im1.get_xticklabels()+sp_im1.get_yticklabels()+
               sp_im2.get_xticklabels()+sp_im2.get_yticklabels()+
               sp_ker.get_xticklabels()+sp_ker.get_yticklabels(), visible=False)
    pylab.setp(sp_his.get_xticklabels()+sp_his.get_yticklabels(), fontsize=fontsize)


def plotSigmaHistograms(info, title, fontsize=10, outfile=None):
    pylab.figure()
    
    # info lists have
    # template image, science image, difference image (in sigma), kernel, sigmas
    # set up plotting windows
    # left, bottom, width, height
    cim0  = [0.050, 0.700, 0.150, 0.150]
    cim1  = [0.200, 0.700, 0.150, 0.150]
    cim2  = [0.350, 0.700, 0.150, 0.150]
    cker  = [0.225, 0.550, 0.100, 0.100]
    chis  = [0.575, 0.550, 0.375, 0.350]

    dcim0 = [0.050, 0.250, 0.150, 0.150]
    dcim1 = [0.200, 0.250, 0.150, 0.150]
    dcim2 = [0.350, 0.250, 0.150, 0.150]
    dcker = [0.225, 0.100, 0.100, 0.100]
    dchis = [0.575, 0.100, 0.375, 0.350]

    makeSigmaHistograms(info[0], title, fontsize, (cim0, cim1, cim2, cker, chis))
    if len(info) == 2:
        makeSigmaHistograms(info[1], title, fontsize, (dcim0, dcim1, dcim2, dcker, dchis))        
        

#    bins   = numpy.arange(-5.0, +5.0, 0.5)
#    theory = []
#    for b in bins:
#        theory.append( numpy.exp( -0.5 * b**2 ) / numpy.sqrt(2. * numpy.pi) )
#
#    #
#    ########
#    #
#    
#    sp_cim0 = pylab.axes(cim0)
#    sp_cim0.imshow( scaleImageForDisplay_nonlinear(convInfo[0]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
#    sp_cim0.set_title('Convolved', fontsize=fontsize, weight='bold')
#    
#    sp_cim1 = pylab.axes(cim1)
#    sp_cim1.imshow( scaleImageForDisplay_nonlinear(convInfo[1]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
#    sp_cim1.set_title('Science', fontsize=fontsize, weight='bold')
#
#    sp_cim2 = pylab.axes(cim2)
#    sp_cim2.imshow( scaleImageForDisplay_nonlinear(convInfo[2]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
#    sp_cim2.set_title('Difference', fontsize=fontsize, weight='bold')
#
#    sp_cker = pylab.axes(cker)
#    sp_cker.imshow( scaleImageForDisplay_nonlinear(convInfo[3]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
#    sp_cker.set_title('Kernel', fontsize=fontsize, weight='bold')
#    
#    sp_chis = pylab.axes(chis)
#    sp_chis.hist(convInfo[4], bins=bins, normed=True)
#    sp_chis.plot(bins, theory, 'r-')
#    sp_chis.set_xlabel('Sigma', fontsize=fontsize, weight='bold')
#    sp_chis.set_ylabel('Fraction', fontsize=fontsize, weight='bold')
#    sp_chis.set_title('%s : %.3f +/- %.3f' % (title, convInfo[4].mean(), convInfo[4].std()), fontsize=fontsize+1, weight='bold')
#
#    pylab.setp(sp_cim0.get_xticklabels()+sp_cim0.get_yticklabels()+
#               sp_cim1.get_xticklabels()+sp_cim1.get_yticklabels()+
#               sp_cim2.get_xticklabels()+sp_cim2.get_yticklabels()+
#               sp_cker.get_xticklabels()+sp_cker.get_yticklabels(), visible=False)
#    pylab.setp(sp_chis.get_xticklabels()+sp_chis.get_yticklabels(), fontsize=fontsize)
#
#    #
#    ########
#    #
#    
#    sp_dcim0 = pylab.axes(dcim0)
#    sp_dcim0.imshow( scaleImageForDisplay_nonlinear(deconvInfo[0]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
#    sp_dcim0.set_title('Convolved', fontsize=fontsize, weight='bold')
#    
#    sp_dcim1 = pylab.axes(dcim1)
#    sp_dcim1.imshow( scaleImageForDisplay_nonlinear(deconvInfo[1]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
#    sp_dcim1.set_title('Science', fontsize=fontsize, weight='bold')
#    
#    sp_dcim2 = pylab.axes(dcim2)
#    sp_dcim2.imshow( scaleImageForDisplay_nonlinear(deconvInfo[2]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
#    sp_dcim2.set_title('Difference', fontsize=fontsize, weight='bold')
#    
#    sp_dcker = pylab.axes(dcker)
#    sp_dcker.imshow( scaleImageForDisplay_nonlinear(deconvInfo[3]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
#    sp_dcker.set_title('Kernel', fontsize=fontsize, weight='bold')
#    
#    sp_dchis = pylab.axes(dchis)
#    sp_dchis.hist(deconvInfo[4], bins=bins, normed=True)
#    sp_dchis.plot(bins, theory, 'r-')
#    sp_dchis.set_xlabel('Sigma', fontsize=fontsize, weight='bold')
#    sp_dchis.set_ylabel('Fraction', fontsize=fontsize, weight='bold')
#    sp_dchis.set_title('%.3f +/- %.3f' % (deconvInfo[4].mean(), deconvInfo[4].std()), fontsize=fontsize+1, weight='bold')
#
#    pylab.setp(sp_dcim0.get_xticklabels()+sp_dcim0.get_yticklabels()+
#               sp_dcim1.get_xticklabels()+sp_dcim1.get_yticklabels()+
#               sp_dcim2.get_xticklabels()+sp_dcim2.get_yticklabels()+
#               sp_dcker.get_xticklabels()+sp_dcker.get_yticklabels(), visible=False)
#    pylab.setp(sp_dchis.get_xticklabels()+sp_dchis.get_yticklabels(), fontsize=fontsize)

    if outfile != None:
        pylab.savefig(outfile)
