import sys, os, re
import pyfits
import lsst.detection.detectionLib as detectionLib
import lsst.afw.image as imageLib

maskFormat = re.compile('^\[(\d+):(\d+),(\d+):(\d+)\]$')
maskPlane  = 'BAD'

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage : createMaskForCFHT.py infile.fits'
        sys.exit(1)
        
    inMask = sys.argv[1]
    if not os.path.isfile(inMask):
        print 'ERROR: Cannot locate', inMask
        sys.exit(1)

    try:
        inPtr  = pyfits.open(inMask)
    except:
        print 'ERROR: Cannot open', inMask
        sys.exit(1)

    nExtensions = len(inPtr)

    # note that the first extension is only a header
    for nExt in range(1,nExtensions):
        header   = inPtr[nExt].header
        nColMask = header['NAXIS1']
        nRowMask = header['NAXIS2']

        # make a new Mask image
        mask = imageLib.MaskUPtr( imageLib.MaskU(nColMask, nRowMask) )
        
        # note that this will have its mask planes initialized by default
        # we want to modify BAD
        badBitMask = mask.getPlaneBitMask(maskPlane)

        # put them all in a list and do at once
        footprintList = detectionLib.FootprintContainerT()
        
        for card in header.ascardlist().keys():
            if card.startswith('MASK_'):
                maskArea = header[card]
                # the convention of the cards is
                # [col_min:row_min,col_max:row_max]
                # inclusive
                match = maskFormat.match(maskArea)
                if match == None:
                    # unable to match mask area!
                    print '# WARNING: Extn', nExt, 'unable to parse', maskArea
                    continue

                group = map(int, match.groups())
                # turn into a Footprint!
                # we have to account for 2 things
                #
                # the regions in the fits file are 1-indexed
                # while the footprints are 0-indexed
                # therefore subtract 1 from the origin
                maskBBox2i = imageLib.BBox2i(group[0]-1,             # col min
                                             group[1]-1,             # row min
                                             group[2]-group[0]+1,    # col span
                                             group[3]-group[1]+1)    # row span
                maskFootprint = detectionLib.FootprintPtrT( detectionLib.Footprint(maskBBox2i) )
                footprintList.push_back(maskFootprint)
                

        # set all the bad masks at once!
        detectionLib.setMaskFromFootprintList(mask, footprintList, badBitMask)

        # save as an output file / persist / or place on clipboard
        outputMaskFile = re.sub('.fits', '_%d_mask.fits' % (nExt), inMask)
        print '# Writing', outputMaskFile
        mask.writeFits(outputMaskFile)
