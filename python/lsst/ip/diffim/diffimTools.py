import numpy
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage

# Temporary functions until we formalize this in the build system somewhere
### NOTE THESE ARE ALSO IN lsst.afw.image.testUtils

def imageToMatrix(im, dtype=float):
    arr = numpy.zeros([im.getCols(), im.getRows()], dtype=dtype)
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            arr[col, row] = im.getVal(col, row)
    return arr

def imageToVector(im, dtype=float):
    arr = numpy.zeros([im.getCols()*im.getRows()], dtype=dtype)
    n   = 0
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            arr[n] = im.getVal(col, row)
            n += 1
    return arr

def matrixToImage(arr):
    im = afwImage.ImageF(arr.shape[0], arr.shape[1])
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            im.set(col, row, arr[col, row])
    return im

def vectorToImage(arr, nCols, nRows):
    assert len(arr) == nCols * nRows
    im = afwImage.ImageF(nCols, nRows)
    n  = 0
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            im.set(col, row, arr[n])
            n += 1
    return im

def matrixToKernelPtr(arr):
    im = afwImage.ImageD(arr.shape[0], arr.shape[1])
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            im.set(col, row, arr[col, row])
    return afwMath.KernelPtr( afwMath.FixedKernel(im) ) 

def vectorToKernelPtr(arr, nCols, nRows):
    # need imageD for FixedKernels
    assert len(arr) == nCols * nRows
    im = afwImage.ImageD(nCols, nRows)
    n  = 0
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            im.set(col, row, arr[n])
            n += 1
    return afwMath.KernelPtr( afwMath.FixedKernel(im) ) 

def vectorPairToVectors(vectorPair):
    kernelVector = afwMath.vectorD()
    kernelErrorVector = afwMath.vectorD()
    for i in range(vectorPair.size()-1):
        #print i, vectorPair[i][0], vectorPair[i][1]
        kernelVector.push_back(vectorPair[i][0])
        kernelErrorVector.push_back(vectorPair[i][1])

    background = vectorPair.back()[0]
    backgroundError = vectorPair.back()[1]
    #print 'B', background, backgroundError
    return kernelVector, kernelErrorVector, background, backgroundError
