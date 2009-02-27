import sys, re, os
sys.path.append('/lsst/home/becker/python/pyfits-1.3/lib/python')
import pyfits

def splitCcd(header, data, infile, ccd, rootdir='/lsst/images/repository/input/'):
    basename = re.sub('.fits', '', os.path.basename(infile))
    basedir  = os.path.join(rootdir, basename)
    raftdir  = os.path.join(basedir, str(ccd))

    if not os.path.isdir(basedir):
        os.mkdir(basedir)
    if not os.path.isdir(raftdir):
        os.mkdir(raftdir)
    
    # We have to undo whatever imsplice has done
    #
    # HISTORY imsplice: FLIPS ver 2.0 - Elixir by CFHT - Thu Jun 5 2003 - 22:11:29
    # HISTORY imsplice: Splicing the two readouts (A&B) per CCD into a unique image
    # HISTORY imsplice: Splicing results in all detectors as if read from A amplifier
    # HISTORY imsplice: NEXTEND keyword updated (/2) = number of CCDs vs. amplifiers
    # HISTORY imsplice: EXTNAME keyword becomes `ccdxx` instead of `ampxx`
    # HISTORY imsplice: AMPNAME keyword now covers both amplifiers (eg `29a + 29b')
    # HISTORY imsplice: keyword {GAIN, RDNOISE, MAXLIN} replaced by two keywords [A,B]
    # HISTORY imsplice: DATASEC and DETSEC keywords now reflect the entire CCD
    # HISTORY imsplice: BIASSEC becomes irrelevant -> replaced by BSECA & BSECB
    # HISTORY imsplice: New keywords are DETSECA, DETSECB, DSECA, DSECB, TSECA, TSECB
    # HISTORY imsplice: New keywords are ASECA, ASECB, CSECA, CSECB

    # http://cfht.hawaii.edu/Instruments/Imaging/MegaPrime/rawdata.html
    
    # Summary:
    # The CCD images are 2k x 4k.
    # We will chop into CFHT amp images 1k x 4k.
    # And then into LSST amp images 1k x 1k.

    # Details:
    # The CCD images are 2112 by 4644.
    #
    # We ignore the overscan on top (TSA)
    # 
    # A pixels  : 1:1056     1:4612
    # B pixels  : 1057:4644  1:4612
    # note that pyfits convetion is y, x
    cfhtAmpA  = data[0:4612, 0:1056]
    cfhtAmpB  = data[0:4612, 1057:4644]

    nPixY = 1153
    for i in range(4):
        y0 = i * nPixY
        y1 = y0 + nPixY

        lsstAmpA = cfhtAmpA[y0:y1, 0:1056]
        lsstAmpB = cfhtAmpB[y0:y1, 0:1056]

        headerA  = header.copy()
        headerB  = header.copy()

        # which cfht amp?
        headerA.update('AMPLIST', 'A')
        headerB.update('AMPLIST', 'B')

        # which lsst amp?
        Aid = i + 1
        Bid = i + 5
        headerA.update('LSSTAMP',  Aid)
        headerB.update('LSSTAMP',  Bid)

        # how did we get these data?
        headerA.update('LSSTSS1', '[0:4612,0:1056]',    'pyfits substamp from original CFHT image')
        headerB.update('LSSTSS1', '[0:4612,1057:4644]', 'pyfits substamp from original CFHT image')

        # again...
        headerA.update('LSSTSS2',  '[%d:%d, 0:1056]' % (y0, y1), 'pyfits substamp from trimmed CFHT image')
        headerB.update('LSSTSS2',  '[%d:%d, 0:1056]' % (y0, y1), 'pyfits substamp from trimmed CFHT image')

        # reset the appropriate readnoise
        headerA.update('RDNOISE',  headerA['RDNOISEA'])
        headerB.update('RDNOISE',  headerB['RDNOISEB'])
        
        # reset the appropriate gain
        headerA.update('GAIN',  headerA['GAINA'])
        headerB.update('GAIN',  headerB['GAINB'])

        # set the appropriate overscan
        headerA.update('BIASSEC', '[1:32,1:%d]'      % (nPixY))
        headerB.update('BIASSEC', '[1025:1056,1:%d]' % (nPixY))

        # set the appropriate datasec
        headerA.update('DATASEC', '[33:1056,1:%d]'   % (nPixY))
        headerB.update('DATASEC', '[1:1024,1:%d]'    % (nPixY))
        
        # set the appropriate trimsec
        headerA.update('TRIMSEC',  headerA['DATASEC'])
        headerB.update('TRIMSEC',  headerB['DATASEC'])

        # set the appropriate crpix1
        headerA.update('CRPIX1',  headerA['CRPIX1'])
        headerB.update('CRPIX1',  headerB['CRPIX1'] - 1024)

        # set the appropriate crpix2
        headerA.update('CRPIX2',  headerA['CRPIX2'] - y0)
        headerB.update('CRPIX2',  headerB['CRPIX2'] - y0)

        # EMPERICAL
        # 33 PIXELS IS OVERSCAN SIZE...
        # Maybe its the top overscan that we just got rid of?
        if i > 17:
            headerA.update('CRPIX1',  headerA['CRPIX1'] - 33)
        else:
            headerB.update('CRPIX1',  headerB['CRPIX1'] - 33)

        outfileA = os.path.join(raftdir, '%s_%d_%s.fits' % (basename, ccd, Aid))
        outfileB = os.path.join(raftdir, '%s_%d_%s.fits' % (basename, ccd, Bid))
        print '# Writing', outfileA
        pyfits.PrimaryHDU(lsstAmpA, headerA).writeto(outfileA, output_verify='silentfix', clobber=True)
        print '# Writing', outfileB
        pyfits.PrimaryHDU(lsstAmpB, headerB).writeto(outfileB, output_verify='silentfix', clobber=True)
    

for infile in sys.argv[1:]:
    ptr          = pyfits.open(infile)
    commonHeader = ptr[0].header

    for i in range(1, 37):
        ccdHeader = ptr[i].header
        ccdData   = ptr[i].data

        # reference image : debugging
        #outfile = 'tmp_%d.fits' % (i)
        #pyfits.PrimaryHDU(ccdData, ccdHeader).writeto(outfile, output_verify='silentfix', clobber=True)
        
        splitCcd(ccdHeader, ccdData, infile, i)
