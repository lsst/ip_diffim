import numpy, pylab

def scaleImageForDisplay(indata):
    return numpy.log10(indata)

def sigmaHistograms(convInfo, deconvInfo, fontsize=10):
    # info lists have
    # template image, science image, difference image (in sigma), kernel, sigmas
    # set up plotting windows
    # left, bottom, width, height
    cim0  = [0.050, 0.700, 0.150, 0.150]
    cim1  = [0.200, 0.700, 0.150, 0.150]
    cim2  = [0.350, 0.700, 0.150, 0.150]
    cker  = [0.225, 0.550, 0.100, 0.100]
    chis  = [0.550, 0.550, 0.400, 0.350]

    dcim0 = [0.050, 0.250, 0.150, 0.150]
    dcim1 = [0.200, 0.250, 0.150, 0.150]
    dcim2 = [0.350, 0.250, 0.150, 0.150]
    dcker = [0.225, 0.100, 0.100, 0.100]
    dchis = [0.550, 0.100, 0.400, 0.350]

    bins   = numpy.arange(-5.0, +5.0, 0.5)
    theory = []
    for b in bins:
        theory.append( numpy.exp( -0.5 * b**2 ) / numpy.sqrt(2. * numpy.pi) )

    #
    ########
    #
    
    sp_cim0 = pylab.axes(cim0)
    sp_cim0.imshow( scaleImageForDisplay(convInfo[0]), cmap=pylab.cm.gray )
    sp_cim0.set_title('Convolved', fontsize=fontsize, weight='bold')
    
    sp_cim1 = pylab.axes(cim1)
    sp_cim1.imshow( scaleImageForDisplay(convInfo[1]), cmap=pylab.cm.gray )
    sp_cim1.set_title('Science', fontsize=fontsize, weight='bold')

    sp_cim2 = pylab.axes(cim2)
    sp_cim2.imshow( convInfo[2], cmap=pylab.cm.gray ) # diffim not log10
    sp_cim2.set_title('Difference', fontsize=fontsize, weight='bold')

    sp_cker = pylab.axes(cker)
    sp_cker.imshow( convInfo[3], cmap=pylab.cm.gray ) # diffim not log10
    sp_cker.set_title('Kernel', fontsize=fontsize, weight='bold')
    
    sp_chis = pylab.axes(chis)
    sp_chis.hist(convInfo[4], bins=bins, normed=True)
    sp_chis.plot(bins, theory, 'r-')
    sp_chis.set_xlabel('Sigma', fontsize=fontsize+1, weight='bold')
    sp_chis.set_ylabel('N', fontsize=fontsize+1, weight='bold')

    pylab.setp(sp_cim0.get_xticklabels()+sp_cim0.get_yticklabels()+
               sp_cim1.get_xticklabels()+sp_cim1.get_yticklabels()+
               sp_cim2.get_xticklabels()+sp_cim2.get_yticklabels()+
               sp_cker.get_xticklabels()+sp_cker.get_yticklabels(), visible=False)
    pylab.setp(sp_chis.get_xticklabels()+sp_chis.get_yticklabels(), fontsize=fontsize)

    #
    ########
    #
    
    sp_dcim0 = pylab.axes(dcim0)
    sp_dcim0.imshow( scaleImageForDisplay(deconvInfo[0]), cmap=pylab.cm.gray )
    sp_dcim0.set_title('Convolved', fontsize=fontsize, weight='bold')
    
    sp_dcim1 = pylab.axes(dcim1)
    sp_dcim1.imshow( scaleImageForDisplay(deconvInfo[1]), cmap=pylab.cm.gray )
    sp_dcim1.set_title('Science', fontsize=fontsize, weight='bold')
    
    sp_dcim2 = pylab.axes(dcim2)
    sp_dcim2.imshow( deconvInfo[2], cmap=pylab.cm.gray )
    sp_dcim2.set_title('Difference', fontsize=fontsize, weight='bold')
    
    sp_dcker = pylab.axes(dcker)
    sp_dcker.imshow( deconvInfo[3], cmap=pylab.cm.gray ) # diffim not log10
    sp_dcker.set_title('Kernel', fontsize=fontsize, weight='bold')
    
    sp_dchis = pylab.axes(dchis)
    sp_dchis.hist(deconvInfo[4], bins=bins, normed=True)
    sp_dchis.plot(bins, theory, 'r-')
    sp_dchis.set_xlabel('Sigma', fontsize=fontsize+1, weight='bold')
    sp_dchis.set_ylabel('N', fontsize=fontsize+1, weight='bold')

    pylab.setp(sp_dcim0.get_xticklabels()+sp_dcim0.get_yticklabels()+
               sp_dcim1.get_xticklabels()+sp_dcim1.get_yticklabels()+
               sp_dcim2.get_xticklabels()+sp_dcim2.get_yticklabels()+
               sp_dcker.get_xticklabels()+sp_dcker.get_yticklabels(), visible=False)
    pylab.setp(sp_dchis.get_xticklabels()+sp_dchis.get_yticklabels(), fontsize=fontsize)

    pylab.show()
