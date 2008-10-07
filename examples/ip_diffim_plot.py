import numpy, pylab

def scaleImageTo8Bit(indata):
    scaled  = indata - indata.min()  # lowest value is 0.0
    scaled /= scaled.max()           # highest value is 1.0
    return numpy.minimum(numpy.floor(numpy.maximum(scaled * 256., 0)), 255)

def scaleImageForDisplay_nonlinear(indata, nonlinearity=2.0):
    return scaleImageTo8Bit( numpy.arcsinh(indata*nonlinearity) / nonlinearity )

def scaleImageForDisplay_log10(indata):
    return numpy.log10(indata)

def plotBackground(backgroundFunction, cols, rows, nbins=100):
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
    sp_bg.plot(cols, rows, 'kx')
    sp_bg.axhline(y=0, c='k', linestyle='--')
    sp_bg.axvline(x=0, c='k', linestyle='--')
    pylab.colorbar(im)

def eigenKernelPlot(convInfo, deconvInfo, fontsize=10):
    pylab.figure()
    
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
    sp_ckm.imshow( scaleImageForDisplay_nonlinear(convInfo[0]), cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ckm.set_title('Mean', fontsize=fontsize, weight='bold')

    sp_ck1 = pylab.axes(ck1)
    sp_ck1.imshow( scaleImageForDisplay_nonlinear(convInfo[1]), cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ck1.set_title('PK1', fontsize=fontsize, weight='bold')

    sp_ck2 = pylab.axes(ck2)
    sp_ck2.imshow( scaleImageForDisplay_nonlinear(convInfo[2]), cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ck2.set_title('PK2', fontsize=fontsize, weight='bold')

    sp_ck3 = pylab.axes(ck3)
    sp_ck3.imshow( scaleImageForDisplay_nonlinear(convInfo[3]), cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ck3.set_xlabel('PK3', fontsize=fontsize, weight='bold')

    sp_ck4 = pylab.axes(ck4)
    sp_ck4.imshow( scaleImageForDisplay_nonlinear(convInfo[4]), cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ck4.set_xlabel('PK4', fontsize=fontsize, weight='bold')

    sp_ck5 = pylab.axes(ck5)
    sp_ck5.imshow( scaleImageForDisplay_nonlinear(convInfo[5]), cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_ck5.set_xlabel('PK5', fontsize=fontsize, weight='bold')

    sp_cspec  = pylab.axes(cspec)
    sp_chist  = numpy.cumsum( convInfo[6] )
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
    sp_dckm.imshow( scaleImageForDisplay_nonlinear(deconvInfo[0]), cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_dckm.set_title('Mean', fontsize=fontsize, weight='bold')

    sp_dck1 = pylab.axes(dck1)
    sp_dck1.imshow( scaleImageForDisplay_nonlinear(deconvInfo[1]), cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_dck1.set_title('PK1', fontsize=fontsize, weight='bold')

    sp_dck2 = pylab.axes(dck2)
    sp_dck2.imshow( scaleImageForDisplay_nonlinear(deconvInfo[2]), cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_dck2.set_title('PK2', fontsize=fontsize, weight='bold')

    sp_dck3 = pylab.axes(dck3)
    sp_dck3.imshow( scaleImageForDisplay_nonlinear(deconvInfo[3]), cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_dck3.set_xlabel('PK3', fontsize=fontsize, weight='bold')

    sp_dck4 = pylab.axes(dck4)
    sp_dck4.imshow( scaleImageForDisplay_nonlinear(deconvInfo[4]), cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_dck4.set_xlabel('PK4', fontsize=fontsize, weight='bold')

    sp_dck5 = pylab.axes(dck5)
    sp_dck5.imshow( scaleImageForDisplay_nonlinear(deconvInfo[5]), cmap=pylab.cm.gray, extent=None, aspect='equal', interpolation='sinc' )
    sp_dck5.set_xlabel('PK5', fontsize=fontsize, weight='bold')

    sp_dcspec  = pylab.axes(dcspec)
    sp_dchist  = numpy.cumsum( deconvInfo[6] )
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

def sigmaHistograms(convInfo, deconvInfo, fontsize=10):
    
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

    bins   = numpy.arange(-5.0, +5.0, 0.5)
    theory = []
    for b in bins:
        theory.append( numpy.exp( -0.5 * b**2 ) / numpy.sqrt(2. * numpy.pi) )

    #
    ########
    #
    
    sp_cim0 = pylab.axes(cim0)
    sp_cim0.imshow( scaleImageForDisplay_nonlinear(convInfo[0]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
    sp_cim0.set_title('Convolved', fontsize=fontsize, weight='bold')
    
    sp_cim1 = pylab.axes(cim1)
    sp_cim1.imshow( scaleImageForDisplay_nonlinear(convInfo[1]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
    sp_cim1.set_title('Science', fontsize=fontsize, weight='bold')

    sp_cim2 = pylab.axes(cim2)
    sp_cim2.imshow( scaleImageForDisplay_nonlinear(convInfo[2]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
    sp_cim2.set_title('Difference', fontsize=fontsize, weight='bold')

    sp_cker = pylab.axes(cker)
    sp_cker.imshow( scaleImageForDisplay_nonlinear(convInfo[3]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
    sp_cker.set_title('Kernel', fontsize=fontsize, weight='bold')
    
    sp_chis = pylab.axes(chis)
    sp_chis.hist(convInfo[4], bins=bins, normed=True)
    sp_chis.plot(bins, theory, 'r-')
    sp_chis.set_xlabel('Sigma', fontsize=fontsize+1, weight='bold')
    sp_chis.set_ylabel('Fraction', fontsize=fontsize+1, weight='bold')

    pylab.setp(sp_cim0.get_xticklabels()+sp_cim0.get_yticklabels()+
               sp_cim1.get_xticklabels()+sp_cim1.get_yticklabels()+
               sp_cim2.get_xticklabels()+sp_cim2.get_yticklabels()+
               sp_cker.get_xticklabels()+sp_cker.get_yticklabels(), visible=False)
    pylab.setp(sp_chis.get_xticklabels()+sp_chis.get_yticklabels(), fontsize=fontsize)

    #
    ########
    #
    
    sp_dcim0 = pylab.axes(dcim0)
    sp_dcim0.imshow( scaleImageForDisplay_nonlinear(deconvInfo[0]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
    sp_dcim0.set_title('Convolved', fontsize=fontsize, weight='bold')
    
    sp_dcim1 = pylab.axes(dcim1)
    sp_dcim1.imshow( scaleImageForDisplay_nonlinear(deconvInfo[1]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
    sp_dcim1.set_title('Science', fontsize=fontsize, weight='bold')
    
    sp_dcim2 = pylab.axes(dcim2)
    sp_dcim2.imshow( scaleImageForDisplay_nonlinear(deconvInfo[2]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
    sp_dcim2.set_title('Difference', fontsize=fontsize, weight='bold')
    
    sp_dcker = pylab.axes(dcker)
    sp_dcker.imshow( scaleImageForDisplay_nonlinear(deconvInfo[3]), cmap=pylab.cm.gray, extent=None, aspect='equal' )
    sp_dcker.set_title('Kernel', fontsize=fontsize, weight='bold')
    
    sp_dchis = pylab.axes(dchis)
    sp_dchis.hist(deconvInfo[4], bins=bins, normed=True)
    sp_dchis.plot(bins, theory, 'r-')
    sp_dchis.set_xlabel('Sigma', fontsize=fontsize+1, weight='bold')
    sp_dchis.set_ylabel('Fraction', fontsize=fontsize+1, weight='bold')

    pylab.setp(sp_dcim0.get_xticklabels()+sp_dcim0.get_yticklabels()+
               sp_dcim1.get_xticklabels()+sp_dcim1.get_yticklabels()+
               sp_dcim2.get_xticklabels()+sp_dcim2.get_yticklabels()+
               sp_dcker.get_xticklabels()+sp_dcker.get_yticklabels(), visible=False)
    pylab.setp(sp_dchis.get_xticklabels()+sp_dchis.get_yticklabels(), fontsize=fontsize)

