import numpy
import lsst.afw.image as afwImage

#######
# Expansions of functionality found in lsst.afw.image.testUtils
#######

def vectorFromImage(im, dtype=float):
    vec = numpy.zeros(im.getWidth()*im.getHeight(), dtype=dtype)
    idx = 0
    for row in range(im.getHeight()):
        for col in range(im.getWidth()):
            vec[idx] = im.get(col, row)
            idx     += 1
    return vec

def imageFromVector(vec, width, height, retType=afwImage.ImageF):
    im  = retType(width, height)
    idx = 0
    for row in range(height):
        for col in range(width):
            # need to cast numpy.float64 as float
            im.set(col, row, float(vec[idx]))
            idx     += 1
    return im

