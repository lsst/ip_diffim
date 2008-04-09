import optparse
import os
import re
import sys

import lsst.afw.image as afwImage

TESTING = False
fitspat = re.compile(r"(\d+)p(_tmpl)?_(\d+)_(img|var|msk).fits")
basepat = re.compile(r"(\d+)p(_tmpl)?_(\d+)")
tailpat = re.compile(r"_(img|var|msk)\.fits$")
imgpat = re.compile(r"(\d+)p(_tmpl)?_(\d+)_img\.fits$")

usage = """usage: %%prog [-c ccd] [e exp] [-r raft] [-o dir] [-d dir] fitsdir|fits"""

cl = optparse.OptionParser(usage)
cl.add_option("-c", "--ccdId", type="int", action="store", dest="forceCCD",
              help="assign output images CCD IDs starting with ccdId",
              default=None, metavar="ccdId")
cl.add_option("-f", "--filter", type="string", action="store", 
              help="assign output images the given CCD ID", dest="forceFilt",
              default=None, metavar="exposureId")
cl.add_option("-e", "--expId", type="string", action="store", dest="forceExp",
              help="assign output images the given Exposure ID", 
              default=None, metavar="exposureId")
cl.add_option("-r", "--raftId",type="string",action="store",dest="forceRaft",
              help="assume images have the given input CCD ID",
              default=None, metavar="ccdId")
cl.add_option("-o", "--outdir", type="string", action="store", dest="outdir",
              help="write images out into the given directory",
              default=None, metavar="dir")
cl.add_option("-d", "--indir", type="string", action="store", dest="indir",
              help="look for images in the given directory", default=".",
              metavar="dir")

def splitImage(inimg, outdir, baseimg, start=1):

    if not TESTING:
        inputMaskedImage = afwImage.MaskedImageF()
        inputMaskedImage.readFits(inimg)
        inputWCS = afwImage.WCS(inputMaskedImage.getImage().getMetaData())
        inputExposure = afwImage.ExposureF(inputMaskedImage, inputWCS)

        nRowSubexposures = 4 # int(sys.argv[2]) # 4
        nColSubexposures = 2 # int(sys.argv[3]) # 2
        nRowMaskedImage = inputMaskedImage.getRows() # 4644
        nColMaskedImage = inputMaskedImage.getCols() # 2112
    
    else:
        nRowSubexposures = 4
        nColSubexposures = 2
        nRowMaskedImage  = 4644
        nColMaskedImage  = 2112

    nRowPix = int(nRowMaskedImage / nRowSubexposures)
    nColPix = int(nColMaskedImage / nColSubexposures)

    ccd = start
    for row in range(nRowSubexposures):
        for col in range(nColSubexposures):
            extn = "%d"
            if ccd < 10:
                extn = "00%d"
            elif ccd < 100:
                extn = "0%d"
            extn %= ccd

            out = os.path.join(outdir,extn);
            if os.path.exists(out) and not os.path.isdir(out):
                raise RuntimeError(out + ": exists but is not a directory")
            if not os.path.exists(out):
                if not TESTING:
                    os.makedirs(out)

            out = os.path.join(out, "".join([baseimg,"_",extn]))
            print '# Writing', "".join([baseimg,"_",extn])
        
            bbox = afwImage.BBox2i(col * nColPix,
                             row * nRowPix,
                             nColPix,
                             nRowPix)

            if not TESTING:
                outputExposure = inputExposure.getSubExposure(bbox)
                #outputExposure = inputMaskedImage.getSubImage(bbox)
                outputExposure.writeFits(out)
            
            ccd += 1

def main():

    (opts, args) = cl.parse_args();

    for arg in args:
        dir = arg
        if not os.path.isabs(dir):
            dir = os.path.join(opts.indir, arg)

        if os.path.isdir(dir):
            files = map(lambda x: tailpat.sub("", x),
                        filter(lambda x: imgpat.search(x), os.listdir(dir)))
        elif tailpat.search(dir) and os.path.exists(dir):
            files = [ tailpat.sub("", os.path.basename(dir)) ]
            dir = os.path.dirname(dir)
        elif os.path.exists(dir+"_img.fits"):
            files = [ os.path.basename(dir) ]
            dir = os.path.dirname(dir)

        if len(files) > 0 and not basepat.match(files[0]) and \
           opts.forceExp is None and opts.forceRaft is None:
            print >> sys.stderr, \
                  "%s: Unable to interpret name (use -r and -e); skipping" \
                  % files[0]
            continue
            
        for file in files:
            if opts.forceExp is None and opts.forceRaft is None:
                match = basepat.match(file)
                if match is None:
                    raise RuntimeError("programmer error: fits file " +
                                       "does not match pattern")

                base = match.group(1) + "p"
            else:
                base = opts.forceExp + "p"
            if match.group(2) is not None:
                base += match.group(2)

            out = opts.outdir;
            if opts.forceFilt is not None:
                out = os.path.join(out, opts.forceFilt)

            inb = os.path.join(dir, file)

            if opts.forceCCD is not None:
                start = opts.forceCCD
            else:
                start = int(match.group(3))
                start = ((start - 1) * 8) + 1

            splitImage(inb, out, base, start)
            
#             inb += "_tmpl"
#             base += "_tmpl"
#             if os.path.exists(inb+"_img.fits"):
#                 splitImage(inb, out, base, start)


if __name__ == "__main__":
    main()

