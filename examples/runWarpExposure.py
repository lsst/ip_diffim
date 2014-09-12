import sys
import optparse
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim
import lsst.daf.base as dafBase
from lsst.pex.logging import Log

def main():
    usage = """runWarpExposure.py refExposure towarpExposure outputExposure"""
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args()
    
    def getArg(ind):
        if ind < len(args):
            return args[ind]
    
    refWcsPath  = getArg(0)
    toWarpPath  = getArg(1)
    warpedPath  = getArg(2)

    if refWcsPath == None or toWarpPath == None or warpedPath == None:
        parser.print_help()
        sys.exit(1)
         
    print 'Reference exposure: ', refWcsPath
    print 'Exposure to be warped: ', toWarpPath
    print 'Output exposure:  ', warpedPath

    refWcsExposure = afwImage.ExposureF(refWcsPath)
    toWarpExposure = afwImage.ExposureF(toWarpPath)

    config = ipDiffim.ImagePsfMatchTask.ConfigClass()
    subconfig = config.kernel.active
    warper = afwMath.Warper.fromConfig(subconfig.warpingConfig)
    warpedExposure = warper.warpExposure(refWcsExposure.getWcs(), 
                                         toWarpExposure,
                                         destBBox = refWcsExposure.getBBox())
    warpedExposure.writeFits(warpedPath)
    
def run():
    Log.getDefaultLog()
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), 'Objects leaked:'
        print dafBase.Citizen_census(dafBase.cout, memId0)

if __name__ == '__main__':
    run()
