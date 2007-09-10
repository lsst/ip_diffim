#! /usr/bin/env python

from Stage import Stage
import lsst.fw.Core.fwLib as fw

class ImageSubtractStage(Stage):
    #------------------------------------------------------------------------
    def preprocess(self):
        self.activeClipboard = self.inputQueue.getNextDataset()

        #keys = self.activeClipboard.getKeys()
        keys = ['Policy', ]

        # get policy
        for key in keys:
            value = self.activeClipboard.get(key)
            print 'Python ImageSubtractStage preprocess(): stageId %i key %s value %s' % (self.stageId, key, value)

        # from policy generate kernelBasisSet
        # e.g. imageproc.generateDeltaFunctionKernelSet()

        #value = self.activeClipboard.put('BasisSet', kernelBasisSet)
        self.activeClipboard['BasisSet'] = kernelBasisSet

    #------------------------------------------------------------------------
    def process(self):
        print 'Python ImageSubtractStage process : _rank %i stageId %d' % (self._rank, self.stageId)
        self.activeClipboard = self.inputQueue.getNextDataset()

        templateImage = self.activeClipboard['TemplateImage']
        inputImage    = self.activeClipboard['InputImage']
        basisSet      = self.activeClipboard['BasisSet']
        policy        = self.activeClipboard['Policy']

        # generate kernelPtr, kernelFunctionPtr, and backgroundFunctionPtr from Policy
        
        imageproc.computePsfMatchingKernelForMaskedImage(templateImage, inputImage, basisSet,
                                                         kernelPtr, kernelFunctionPtr, backgroundFunctionPtr, policy)

        differenceImage = inputImage - fw.convolve(templateImage, kernelPtr)

        # Ptr or not Ptr, that is the question.  The answer is Not
        
        self.activeClipboard['DifferenceImage'] = differenceImage
        self.activeClipboard['OutputKernel'] = kernel
        self.activeClipboard['BackgroundModel'] = backgroundFunction
        
    #------------------------------------------------------------------------
    def postprocess(self):
        print 'Python ImageSubtractStage postprocess : stageId %d' % self.stageId
        self.outputQueue.addDataset(self.activeClipboard)

        # potentially delete 'BasisSet'
