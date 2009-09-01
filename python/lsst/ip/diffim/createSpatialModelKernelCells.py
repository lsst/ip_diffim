# all the c++ level classes and routines
import diffimLib

# all the other diffim routines
import diffimDebug

# all the other LSST packages
import lsst.afw.math  as afwMath
import lsst.afw.image as afwImage
import lsst.afw.display.utils as displayUtils
import lsst.afw.display.ds9 as ds9

def createSpatialModelKernelCells(templateMaskedImage,
                                  scienceMaskedImage,
                                  fpInList,
                                  kFunctor,
                                  policy,
                                  cFlag='c',
                                  display=False):

    nSegmentCol = policy.get('nSegmentCol')
    nSegmentRow = policy.get('nSegmentRow')

    nSegmentColPix = int( templateMaskedImage.getWidth() / nSegmentCol )
    nSegmentRowPix = int( templateMaskedImage.getHeight() / nSegmentRow )

    spatialCells   = diffimLib.VectorSpatialModelCellF()
    
    if display:
        stamps = []; stampInfo = []
        imagePairMosaic = displayUtils.Mosaic(mode="x")
        imagePairMosaic.setGutter(2)
        imagePairMosaic.setBackground(-10)

    cellCount = 0
    for col in range(nSegmentCol):
        colMin    = max(0, col*nSegmentColPix)
        colMax    = min(templateMaskedImage.getWidth(), (col+1)*nSegmentColPix)
        colCenter = int( 0.5 * (colMin + colMax) )
        
        for row in range(nSegmentRow):
            rowMin     = max(0, row*nSegmentRowPix)
            rowMax     = min(templateMaskedImage.getHeight(), (row+1)*nSegmentRowPix)
            rowCenter  = int( 0.5 * (rowMin + rowMax) )
            label      = 'c%d' % cellCount

            modelList = diffimLib.VectorSpatialModelKernelF()

            # This is a bit dumb and could be more clever
            # Should never really have a loop within a loop within a loop
            # But we will not have *that many* Footprints...

            for fpID, fpPtr in enumerate(fpInList):
                
                fpBBox = afwImage.BBox(fpPtr.getBBox().getLLC() - templateMaskedImage.getXY0(),
                                       fpPtr.getBBox().getWidth(), fpPtr.getBBox().getHeight())
                
                fpColC = 0.5*(fpBBox.getX0() + fpBBox.getX1())
                fpRowC = 0.5*(fpBBox.getY0() + fpBBox.getY1())

                if col == 0 and row == 0:
                    if display:
                        ds9.dot("+", fpColC, fpRowC, frame=1)

                if (fpColC >= colMin) and (fpColC < colMax) and (fpRowC >= rowMin) and (fpRowC < rowMax):

                    tSubImage = afwImage.MaskedImageF(templateMaskedImage, fpBBox)
                    iSubImage = afwImage.MaskedImageF(scienceMaskedImage, fpBBox)
                    model = diffimLib.SpatialModelKernelF(fpPtr,
                                                          tSubImage,
                                                          iSubImage,
                                                          kFunctor,
                                                          policy,
                                                          False)
                    if policy.get('debugIO'):
                        diffimDebug.writeDiffImages(cFlag, fpID, model)

                    if display:
                        tmpScience  = iSubImage.Factory(iSubImage, True) # make a copy
                        tmpTemplate = tSubImage.Factory(tSubImage, True) # make a copy

                        tmpScience -=  afwMath.makeStatistics(iSubImage.getImage(), afwMath.MEDIAN).getValue()
                        tmpTemplate -= afwMath.makeStatistics(tSubImage.getImage(), afwMath.MEDIAN).getValue()
                        tmpTemplate *= afwMath.makeStatistics(tmpScience.getImage(), afwMath.MAX).getValue()/ \
                                       afwMath.makeStatistics(tmpTemplate.getImage(), afwMath.MAX).getValue()

                        if not model.isBuilt():
                            model.buildModel()

                        if model.getKernelPtr():
                            tmpKernelImage = afwImage.ImageD(model.getKernelPtr().getDimensions())
                            kSum = model.getKernelPtr().computeImage(tmpKernelImage, False)
                            tmpKernelImage = tmpScience.Factory(tmpKernelImage.convertFloat(),
                                                                afwImage.MaskU(tmpKernelImage.getDimensions()),
                                                                afwImage.ImageF(tmpKernelImage.getDimensions())
                                                                )
                            tmpKernelImage *= afwMath.makeStatistics(tmpScience.getImage(), afwMath.MAX).getValue()/ \
                                              afwMath.makeStatistics(tmpKernelImage.getImage(), afwMath.MAX).getValue()
                        else:
                            tmpKernelImage = tmpScience.Factory(tmpScience.getDimensions())
                            tmpKernelImage.set(-10)                            
                        tmpKernelImage.getMask().set(0x0)
                        
                        stamps.append(imagePairMosaic.makeMosaic([tmpTemplate, tmpScience, tmpKernelImage]))
                            
                        stampInfo.append("%s %d" % (label, fpID))

                    modelList.push_back( model )

            spatialCell = diffimLib.SpatialModelCellF(label, colCenter, rowCenter, modelList)
            spatialCells.push_back(spatialCell)

            # Formatting to the screen 
            Trace('lsst.ip.diffim.createSpatialModelKernelCells', 2, '')

            cellCount += 1

    if display and len(stamps) > 0:
        mos = displayUtils.Mosaic()
        mos.setGutter(5)
        mos.setBackground(0)

        ds9.mtv(mos.makeMosaic(stamps), frame=2)
        mos.drawLabels(stampInfo, frame=2)

    return spatialCells


