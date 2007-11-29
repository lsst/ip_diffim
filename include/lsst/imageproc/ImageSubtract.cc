// -*- lsst-c++ -*-
/**
 * \file
 *
 * \brief Implementation of image subtraction functions declared in ImageSubtract.h
 *
 * This file is meant to be included by lsst/imageproc/ImageSubtract.h
 *
 * \author Andrew Becker
 *
 * \ingroup imageproc
 */
#include <iostream>
#include <limits>
#include <boost/timer.hpp> 

#include <vw/Math/Functions.h> 
#include <vw/Math/LinearAlgebra.h> 
#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 

#include <lsst/mwi/exceptions/Exception.h>
#include <lsst/mwi/utils/Trace.h>
#include <lsst/fw/FunctionLibrary.h>
#include <lsst/fw/KernelFunctions.h>
#include <lsst/fw/PixelAccessors.h>
#include <lsst/fw/minimize.h>

#include <lsst/imageproc/PCA.h>

//#define DEBUG_MATRIX 1
//#define DEBUG_TRACE 1

/**
 * \brief Compute spatially varying PSF matching kernel for image subtraction
 *
 * Implements the main use case of Image Subtraction.  Solves for kernel K and background B in the model
 * 
 *   T.conv.K + B = I
 * 
 * The difference image is (I - B - T.conv.K)
 *
 * \return the spatially varying PSF-matching kernel for image subtraction
 */
template <typename ImageT, typename MaskT, typename KernelT>
boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > lsst::imageproc::computePsfMatchingKernelForMaskedImage(
    boost::shared_ptr<lsst::fw::function::Function2<double> > &kernelFunctionPtr, ///< Function for spatial variation of kernel
    boost::shared_ptr<lsst::fw::function::Function2<double> > &backgroundFunctionPtr, ///< Function for spatial variation of background
    lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve, ///< Image to convolve (e.g. Template image)
    lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve, ///< Image to not convolve (e.g. Science image)
    std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelInBasisList, ///< Input set of basis kernels
    std::vector<lsst::detection::Footprint::PtrType> const &footprintList, ///< List of input footprints to do diffim on
    lsst::mwi::policy::Policy &policy ///< Policy directing the behavior
    ) {
    
    lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForMaskedImage", 1, 
                            "WARNING : Entering DEPRECATED subroutine computePsfMatchingKernelForMaskedImage");

   
    // Parse the Policy
    double maximumFootprintResidualMean = policy.getDouble("computePsfMatchingKernelForMaskedImage.maximumFootprintResidualMean");
    double maximumFootprintResidualVariance = policy.getDouble("computePsfMatchingKernelForMaskedImage.maximumFootprintResidualVariance");
    double maximumKernelSumNSigma = policy.getDouble("computePsfMatchingKernelForMaskedImage.maximumKernelSumNSigma");
    unsigned int maximumIterationsKernelSum = policy.getInt("computePsfMatchingKernelForMaskedImage.maximumIterationsKernelSum");
    bool debugIO = policy.getBool("debugIO", false);
    
    //int badMaskBit = policy.getInt("badMaskBit");
    //int edgeMaskBit = policy.getInt("edgeMaskBit");

    /* grab mask bits from the image to convolve, since that is what we'll be operating on */
    int badMaskBit = imageToConvolve.getMask()->getMaskPlane("BAD");
    int edgeMaskBit = imageToConvolve.getMask()->getMaskPlane("EDGE");

    MaskT edgePixelMask = (edgeMaskBit < 0) ? 0 : (1 << edgeMaskBit);
    MaskT badPixelMask = (badMaskBit < 0) ? 0 : (1 << badMaskBit);
    badPixelMask |= edgePixelMask;

    // Reusable view around each footprint
    typename lsst::fw::MaskedImage<ImageT, MaskT>::MaskedImagePtrT imageToConvolveStampPtr;
    typename lsst::fw::MaskedImage<ImageT, MaskT>::MaskedImagePtrT imageToNotConvolveStampPtr;
    typedef typename std::vector<lsst::detection::Footprint::PtrType>::const_iterator iFootprint;
    typedef typename std::vector<lsst::imageproc::DiffImContainer<KernelT> >::iterator iDiffImContainer;
    
    // Reusable view inside each convolved, subtracted footprint
    typename lsst::fw::MaskedImage<ImageT, MaskT>::MaskedImagePtrT convolvedImageStampSubPtr;
    
    // Information for each Footprint
    std::vector<lsst::imageproc::DiffImContainer<KernelT> > diffImContainerList;
    
    // Kernel information
    double imSum;
    /* NOTE - this assumes that all kernels are the same size */
    /*        and that this size is defined in the policy */ 
    unsigned int kernelRows = policy.getInt("kernelRows");
    unsigned int kernelCols = policy.getInt("kernelCols");
    lsst::fw::Image<KernelT> kImage(kernelCols, kernelRows);
    
    // Iterate over footprint
    int nFootprint = 0;
    for (iFootprint i = footprintList.begin(); i != footprintList.end(); ++i) {
        
        BBox2i footprintBBox = (*i)->getBBox();
        imageToConvolveStampPtr    = imageToConvolve.getSubImage(footprintBBox);
        imageToNotConvolveStampPtr = imageToNotConvolve.getSubImage(footprintBBox);
        
        if (debugIO) {
            imageToConvolveStampPtr->writeFits( (boost::format("tFits_%d") % nFootprint).str() );
            imageToNotConvolveStampPtr->writeFits( (boost::format("sFits_%d") % nFootprint).str() );
        }
        
        // Find best single kernel for this stamp
        double background;
        std::vector<double> kernelCoeffs = lsst::imageproc::computePsfMatchingKernelForPostageStamp(
            background, *imageToConvolveStampPtr, *imageToNotConvolveStampPtr, kernelInBasisList, policy);
        
        // Create a linear combination kernel from this 
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > footprintKernelPtr(
            new lsst::fw::LinearCombinationKernel<KernelT>(kernelInBasisList, kernelCoeffs)
            );
        
        // Container to carry around in imageproc.diffim
        lsst::imageproc::DiffImContainer<KernelT> diffImFootprintContainer;
        diffImFootprintContainer.id = nFootprint;
        diffImFootprintContainer.isGood = true;
        diffImFootprintContainer.diffImFootprintPtr = *i;
        diffImFootprintContainer.diffImKernelPtr = footprintKernelPtr;
        diffImFootprintContainer.background = background;
        
        /* get kernel sum */
        imSum = 0.0;
        footprintKernelPtr->computeImage(kImage, imSum, 0.0, 0.0, false);
        diffImFootprintContainer.kSum = imSum;
        
        vw::math::Vector<vw::int32, 2> ctrDoubled = footprintBBox.min() + footprintBBox.max();
        vw::math::Vector<double, 2> center;
        for (int ii = 0; ii < 2; ++ii) {
            center[ii] = static_cast<double>(ctrDoubled[ii]) / 2.0;
        }
        
        // NOTE - do we normalize the coordinates to go from -1 to 1?  If so, use directly below
        //diffImFootprintContainer.colcNorm = 2 * center[0] / imageToConvolve.getCols() - 1;
        //diffImFootprintContainer.rowcNorm = 2 * center[1] / imageToConvolve.getRows() - 1;
        
        // NOTE - else, use actual image coords
        diffImFootprintContainer.colcNorm = center[0];
        diffImFootprintContainer.rowcNorm = center[1];
        lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForMaskedImage", 4, 
                                (boost::format("Developing kernel %d at position %.3f %.3f.  Background = %.3f") 
                                 % diffImFootprintContainer.id % center[0] % center[1] % background));
        
        // QA - calculate the residual of the subtracted image here
        lsst::fw::MaskedImage<ImageT, MaskT>
            convolvedImageStamp = lsst::fw::kernel::convolve(*imageToConvolveStampPtr, 
                                                             *footprintKernelPtr, 
                                                             edgeMaskBit, 
                                                             false);
        
        // Add in background
        convolvedImageStamp += background;   
        
        if (debugIO) {
            kImage.writeFits( (boost::format("kFits_%d.fits") % nFootprint).str() );
            convolvedImageStamp.writeFits( (boost::format("cFits_%d") % nFootprint).str() );
        }
        
        // Do actual subtraction
        convolvedImageStamp -= (*imageToNotConvolveStampPtr);
        convolvedImageStamp *= -1.0;
        
        if (debugIO) {
            convolvedImageStamp.writeFits( (boost::format("dFits_%d") % nFootprint).str() );
        }
        
        /* calculate stats in the difference image */
        double meanOfResiduals = 0.0;
        double varianceOfResiduals = 0.0;
        int nGoodPixels = 0;
        lsst::imageproc::calculateMaskedImageResiduals(nGoodPixels, meanOfResiduals, varianceOfResiduals, convolvedImageStamp, badPixelMask);
        diffImFootprintContainer.footprintResidualMean = meanOfResiduals;
        diffImFootprintContainer.footprintResidualVariance = varianceOfResiduals/nGoodPixels; // Per pixel
        
        if (!(std::isfinite(diffImFootprintContainer.footprintResidualMean)) || 
            !(std::isfinite(diffImFootprintContainer.footprintResidualVariance))) {
            lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForMaskedImage", 4, 
                                    (boost::format("# Kernel %d (kSum=%.3f), analysis of footprint failed : %.3f %.3f (%d pixels)") 
                                     % diffImFootprintContainer.id 
                                     % diffImFootprintContainer.kSum
                                     % diffImFootprintContainer.footprintResidualMean 
                                     % diffImFootprintContainer.footprintResidualVariance 
                                     % nGoodPixels));
            diffImFootprintContainer.isGood = false;
        }
        else if (fabs(diffImFootprintContainer.footprintResidualMean) > maximumFootprintResidualMean) {
            lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForMaskedImage", 4, 
                                    (boost::format("# Kernel %d (kSum=%.3f), bad mean residual of footprint : %.3f (%d pixels)") 
                                     % diffImFootprintContainer.id 
                                     % diffImFootprintContainer.kSum
                                     % diffImFootprintContainer.footprintResidualMean 
                                     % nGoodPixels));
            diffImFootprintContainer.isGood = false;
        }
        else if (diffImFootprintContainer.footprintResidualVariance > maximumFootprintResidualVariance) {
            lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForMaskedImage", 4, 
                                    (boost::format("# Kernel %d (kSum=%.3f), bad residual variance of footprint : %.3f (%d pixels)") 
                                     % diffImFootprintContainer.id 
                                     % diffImFootprintContainer.kSum
                                     % diffImFootprintContainer.footprintResidualVariance 
                                     % nGoodPixels));
            diffImFootprintContainer.isGood = false;
        }
        else {
            lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForMaskedImage", 4, 
                                    (boost::format("Kernel %d (kSum=%.3f), mean and variance of residuals in footprint : %.3f %.3f (%d pixels)") 
                                     % diffImFootprintContainer.id 
                                     % diffImFootprintContainer.kSum
                                     % diffImFootprintContainer.footprintResidualMean 
                                     % diffImFootprintContainer.footprintResidualVariance 
                                     % nGoodPixels));
        }
        
        diffImContainerList.push_back(diffImFootprintContainer);
        nFootprint += 1;
    }
    
    /* 
       Run test comparing the kernel sums; things that are bad (e.g. variable
       stars, moving objects, cosmic rays) will have a kernel sum far from the
       mean.  Reject these.  
    */
    lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForMaskedImage", 3, "Rejecting kernels with deviant kSums");
    
    vw::math::Vector<double> kernelSumVector(nFootprint);
    int kernelCount;

    unsigned int nReject = 1;
    unsigned int nIter = 0;
    double kSumMean, kSumStdev;
    double kMeanMean, kMeanVariance;
    while ( (nIter < maximumIterationsKernelSum) and (nReject != 0) ) {
        kernelSumVector.set_size(nFootprint);  /* erases old data */
        kernelCount = 0;
        kMeanMean = 0.0;
        kMeanVariance = 0.0;
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {    
            if ((*i).isGood == true) {
                kernelSumVector(kernelCount) = (*i).kSum;
                kMeanMean     += (*i).footprintResidualMean;
                kMeanVariance += (*i).footprintResidualVariance;
                kernelCount   += 1;
            }
        }
        kernelSumVector.set_size(kernelCount, true); /* only keep data with entries */
        
        calculateVectorStatistics(kernelSumVector, kSumMean, kSumStdev);
        kSumStdev = sqrt(kSumStdev);
        
        kMeanMean     /= kernelCount; /* mean of the mean residual across all kernels */
        kMeanVariance /= kernelCount; /* mean of the variance in residuals across all kernels */

        nReject = 0;
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {    
            if ((*i).isGood == true) {
                if (fabs( (*i).kSum - kSumMean ) > (kSumStdev * maximumKernelSumNSigma) ) {
                    /* The kernel sum is too large/small for some reason; reject it */
                    (*i).isGood = false;
                    nReject += 1;
                    
                    lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForMaskedImage", 5, 
                                            (boost::format("# Kernel %d (kSum=%.3f) REJECTED due to bad kernel sum (mean=%.3f, std=%.3f)")
                                             % (*i).id % (*i).kSum % kSumMean % kSumStdev));
                }
            }
        }
        
        lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForMaskedImage", 4, 
                                (boost::format("Iteration %d : Kernel Sum = %.3f (%.3f); Average Mean Residual (%.3f); Average Residual Variance (%.3f)")
                                 % nIter % kSumMean % kSumStdev % kMeanMean % kMeanVariance));
        nIter += 1;
    } // End of iteration
    
    
    
    // In the end we want to test if the kernelInBasisList is Delta Functions; if so, do PCA
    // For DC2 we just do it
    std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > kernelOutBasisList =
        lsst::imageproc::computePcaKernelBasis(diffImContainerList, policy);
    
    // Compute spatial variation of the kernel 
    return lsst::imageproc::computeSpatiallyVaryingPsfMatchingKernel(
        kernelFunctionPtr,
        backgroundFunctionPtr,
        diffImContainerList,
        kernelOutBasisList, 
        policy);
}

/**
 * \brief Compute a psf-matching kernel for a postage stamp
 */
template <typename ImageT, typename MaskT, typename KernelT>
std::vector<double> lsst::imageproc::computePsfMatchingKernelForPostageStamp(
    double &background, ///< Difference in the backgrounds
    lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve, ///< Image to convolve
    lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve, ///< Image to not convolve
    std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelInBasisList, ///< Input kernel basis set
    lsst::mwi::policy::Policy &policy ///< Policy directing the behavior
    ) { 
    
    /* grab mask bits from the image to convolve, since that is what we'll be operating on */
    int edgeMaskBit = imageToConvolve.getMask()->getMaskPlane("EDGE");
    
    int nKernelParameters = 0;
    int nBackgroundParameters = 0;
    int nParameters = 0;

    boost::timer t;
    double time;
    t.restart();

    lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForPostageStamp", 3, 
                            "Entering subroutine computePsfMatchingKernelForPostageStamp");
    
    // We assume that each kernel in the Set has 1 parameter you fit for
    nKernelParameters = kernelInBasisList.size();
    // Or, we just assume that across a single kernel, background 0th order.  This quite makes sense.
    nBackgroundParameters = 1;
    // Total number of parameters
    nParameters = nKernelParameters + nBackgroundParameters;
    
    vw::math::Vector<double> B(nParameters);
    vw::math::Matrix<double> M(nParameters, nParameters);
    for (unsigned int i = nParameters; i--;) {
        B(i) = 0;
        for (unsigned int j = nParameters; j--;) {
            M(i,j) = 0;
        }
    }
    
    // convolve creates a MaskedImage, push it onto the back of the Vector
    // need to use shared pointers because MaskedImage copy does not work
    std::vector<boost::shared_ptr<lsst::fw::MaskedImage<ImageT, MaskT> > > convolvedImageList(nKernelParameters);
    // and an iterator over this
    typename std::vector<boost::shared_ptr<lsst::fw::MaskedImage<ImageT, MaskT> > >::iterator citer = convolvedImageList.begin();
    
    // Iterator for input kernel basis
    typename std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > >::const_iterator kiter = kernelInBasisList.begin();
    // Create C_ij in the formalism of Alard & Lupton
    for (; kiter != kernelInBasisList.end(); ++kiter, ++citer) {
        
        lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForPostageStamp", 6, 
                                "Convolving an Object with Basis");
        
        // NOTE : we could also *precompute* the entire template image convolved with these functions
        //        and save them somewhere to avoid this step each time.  however, our paradigm is to
        //        compute whatever is needed on the fly.  hence this step here.
        boost::shared_ptr<lsst::fw::MaskedImage<ImageT, MaskT> > imagePtr(
            new lsst::fw::MaskedImage<ImageT, MaskT>
            (lsst::fw::kernel::convolve(imageToConvolve, **kiter, edgeMaskBit, false))
            );
        
        lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForPostageStamp", 6, 
                                "Convolved an Object with Basis");
        
        *citer = imagePtr;
        
    } 
    
    // NOTE ABOUT CONVOLUTION :
    // getCtrCol:getCtrRow pixels are masked on the left:bottom side
    // getCols()-getCtrCol():getRows()-getCtrRow() masked on right/top side
    // 
    // The convolved image and the input image are by default the same size, so
    // we offset our initial pixel references by the same amount
    kiter = kernelInBasisList.begin();
    citer = convolvedImageList.begin();
    unsigned int startCol = (*kiter)->getCtrCol();
    unsigned int startRow = (*kiter)->getCtrRow();
    // NOTE - I determined I needed this +1 by eye
    unsigned int endCol   = (*citer)->getCols() - ((*kiter)->getCols() - (*kiter)->getCtrCol()) + 1;
    unsigned int endRow   = (*citer)->getRows() - ((*kiter)->getRows() - (*kiter)->getCtrRow()) + 1;
    // NOTE - we need to enforce that the input images are large enough
    // How about some multiple of the PSF FWHM?  Or second moments?
    
    // An accessor for each convolution plane
    // NOTE : MaskedPixelAccessor has no empty constructor, therefore we need to push_back()
    std::vector<lsst::fw::MaskedPixelAccessor<ImageT, MaskT> > convolvedAccessorRowList;
    for (citer = convolvedImageList.begin(); citer != convolvedImageList.end(); ++citer) {
        convolvedAccessorRowList.push_back(lsst::fw::MaskedPixelAccessor<ImageT, MaskT>(**citer));
    }
    
    // An accessor for each input image; address rows and cols separately
    lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imageToConvolveRow(imageToConvolve);
    lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imageToNotConvolveRow(imageToNotConvolve);
    
    // Take into account buffer for kernel images
    imageToConvolveRow.advance(startCol, startRow);
    imageToNotConvolveRow.advance(startCol, startRow);
    for (int ki = 0; ki < nKernelParameters; ++ki) {
        convolvedAccessorRowList[ki].advance(startCol, startRow);
    }

    for (unsigned int row = startRow; row < endRow; ++row) {
        
        // An accessor for each convolution plane
        std::vector<lsst::fw::MaskedPixelAccessor<ImageT, MaskT> > convolvedAccessorColList = convolvedAccessorRowList;
        
        // An accessor for each input image; places the col accessor at the correct row
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imageToConvolveCol = imageToConvolveRow;
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imageToNotConvolveCol = imageToNotConvolveRow;
        
        
        /* NOTE on timing here.
           
        Systematic check
        0) Baseline for each column: 1.36s
        1) Cut out accessing the iterators : 1.14s
        2) Cut out stepping accessors : 1.14s
        3) ITS THE BOOST FORMATTING AGAIN!!!!!!!
        
        
        Redo the check
        Baseline for the convolve above : 7.4s
        Baseline for each row : 0.05s * e.g. 52s = 2.7s
        Baseline for pseudoinverse : 2.1s 
        */
        
        for (unsigned int col = startCol; col < endCol; ++col) {
            
            ImageT ncCamera   = *imageToNotConvolveCol.image;
            ImageT ncVariance = *imageToNotConvolveCol.variance;
            MaskT  ncMask     = *imageToNotConvolveCol.mask;
            
            ImageT cVariance  = *imageToConvolveCol.variance;
            
            // Variance for a particlar pixel; do we use this variance of the
            // input data, or include the variance after its been convolved with
            // the basis?  For now, use the average of the input varianes.
            ImageT iVariance  = 1.0 / (cVariance + ncVariance);
            
            /* TESTING */
            //iVariance = 1.0;
            
            lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForPostageStamp", 6, 
                                    boost::format("Accessing image row %d col %d : %.3f %.3f %d") 
                                    % row % col % ncCamera % ncVariance % ncMask);
            
            // kernel index i
            typename std::vector<lsst::fw::MaskedPixelAccessor<ImageT, MaskT> >::iterator
                kiteri = convolvedAccessorColList.begin();
            
            for (int kidxi = 0; kiteri != convolvedAccessorColList.end(); ++kiteri, ++kidxi) {
                ImageT cdCamerai   = *kiteri->image;
                
                /* NOTE - These inner trace statements can ENTIRELY kill the run time */
#if defined(DEBUG_TRACE)
                ImageT cdVariancei = *kiteri->variance;
                MaskT  cdMaski     = *kiteri->mask;
                lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForPostageStamp", 6, 
                                        boost::format("Accessing convolved image %d : %.3f %.3f %d")
                                        % kidxi % cdCamerai % cdVariancei % cdMaski);
#endif
                
                // kernel index j 
                typename std::vector<lsst::fw::MaskedPixelAccessor<ImageT, MaskT> >::iterator kiterj = kiteri;
                for (int kidxj = kidxi; kiterj != convolvedAccessorColList.end(); ++kiterj, ++kidxj) {
                    ImageT cdCameraj   = *kiterj->image;
                    
                    /* NOTE - These inner trace statements can ENTIRELY kill the run time */
#if defined(DEBUG_TRACE)
                    ImageT cdVariancej = *kiterj->variance;
                    MaskT  cdMaskj     = *kiterj->mask;
                    lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForPostageStamp", 6, 
                                            boost::format("Accessing convolved image %d : %.3f %.3f %d")
                                            % kidxj % cdCameraj % cdVariancej % cdMaskj);
#endif
                    
                    M[kidxi][kidxj] += cdCamerai * cdCameraj * iVariance;
                } 
                
                B[kidxi] += ncCamera * cdCamerai * iVariance;
                
                // Constant background term; effectively j=kidxj+1
                M[kidxi][nParameters-1] += cdCamerai * iVariance;
            } 
            
            // Background term; effectively i=kidxi+1
            B[nParameters-1] += ncCamera * iVariance;
            M[nParameters-1][nParameters-1] += 1.0 * iVariance;
            
            lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForPostageStamp", 6, 
                                    boost::format("Background terms : %.3f %.3f") 
                                    % B[nParameters-1] % M[nParameters-1][nParameters-1]);
            
            // Step each accessor in column
            imageToConvolveCol.nextCol();
            imageToNotConvolveCol.nextCol();
            for (int ki = 0; ki < nKernelParameters; ++ki) {
                convolvedAccessorColList[ki].nextCol();
            }             
            
        } // col
        
        // Step each accessor in row
        imageToConvolveRow.nextRow();
        imageToNotConvolveRow.nextRow();
        for (int ki = 0; ki < nKernelParameters; ++ki) {
            convolvedAccessorRowList[ki].nextRow();
        }
        
    } // row

    // NOTE: If we are going to regularize the solution to M, this is the place to do it
    
    // Fill in rest of M
    for (int kidxi=0; kidxi < nParameters; ++kidxi) 
        for (int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) 
            M[kidxj][kidxi] = M[kidxi][kidxj];
    
#if defined(DEBUG_MATRIX)
    std::cout << "B : " << B << std::endl;
    std::cout << "M : " << M << std::endl;
#endif

    time = t.elapsed();
    lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForPostageStamp", 5, 
                            boost::format("Total compute time before matrix inversions : %.2f s") 
                            % time);


    // Invert using SVD and Pseudoinverse
    vw::math::Matrix<double> Minv;
    Minv = vw::math::pseudoinverse(M);
    //Minv = vw::math::inverse(M);
    
#if defined(DEBUG_MATRIX)
    std::cout << "Minv : " << Minv << std::endl;
#endif
    
    // Solve for x in Mx = B
    vw::math::Vector<double> Soln = Minv * B;
    
#if defined(DEBUG_MATRIX)
    std::cout << "Solution : " << Soln << std::endl;
#endif
    
    time = t.elapsed();
    lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForPostageStamp", 5, 
                            boost::format("Total compute time after matrix inversions : %.2f s") 
                            % time);

    // Translate from VW std::vectors to std std::vectors
    std::vector<double> kernelCoeffs(kernelInBasisList.size());
    for (int ki = 0; ki < nKernelParameters; ++ki) {
        kernelCoeffs[ki] = Soln[ki];
    }
    background = Soln[nParameters-1];
    return kernelCoeffs;
}

/**
 * \brief Find a list of footprints suitable for difference imaging
 *
 * \return a list of footprints suitable for difference imaging
 */
template <typename ImageT, typename MaskT>
std::vector<lsst::detection::Footprint::PtrType> lsst::imageproc::getCollectionOfFootprintsForPsfMatching(
    lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve, ///< Template image; is convolved
    lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve, ///< Science image; is not convolved
    lsst::mwi::policy::Policy &policy ///< Policy directing the behavior
    ) {
    
    // Parse the Policy
    unsigned int footprintDiffimNpixMin = policy.getInt("getCollectionOfFootprintsForPsfMatching.footprintDiffimNpixMin");
    unsigned int footprintDiffimGrow = policy.getInt("getCollectionOfFootprintsForPsfMatching.footprintDiffimGrow");
    double footprintDetectionThreshold = policy.getDouble("getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold");
    int minimumCleanFootprints = policy.getInt("getCollectionOfFootprintsForPsfMatching.minimumCleanFootprints");
    double detectionThresholdScaling = policy.getDouble("getCollectionOfFootprintsForPsfMatching.detectionThresholdScaling");
    double minimumDetectionThreshold = policy.getDouble("getCollectionOfFootprintsForPsfMatching.minimumDetectionThreshold");

    /* grab mask bits from the image to convolve, since that is what we'll be operating on */
    int badMaskBit = imageToConvolve.getMask()->getMaskPlane("BAD");
    MaskT badPixelMask = (badMaskBit < 0) ? 0 : (1 << badMaskBit);

    // Reusable view of each Footprint
    typename lsst::fw::MaskedImage<ImageT, MaskT>::MaskedImagePtrT imageToConvolveFootprintPtr;
    typename lsst::fw::MaskedImage<ImageT, MaskT>::MaskedImagePtrT imageToNotConvolveFootprintPtr;

    // Reusable list of Footprints
    std::vector<lsst::detection::Footprint::PtrType> footprintListIn;
    std::vector<lsst::detection::Footprint::PtrType> footprintListOut;

    int nCleanFootprints = 0;
    while ( (nCleanFootprints < minimumCleanFootprints) and (footprintDetectionThreshold > minimumDetectionThreshold) ) {
        footprintListIn.clear();
        footprintListOut.clear();
        
        /* Find detections */
        lsst::detection::DetectionSet<ImageT, MaskT> 
            detectionSet(imageToConvolve, lsst::detection::Threshold(footprintDetectionThreshold, lsst::detection::Threshold::VALUE));
        
        /* get the footprints */
        footprintListIn = detectionSet.getFootprints();
        
        nCleanFootprints = 0;
        for (std::vector<lsst::detection::Footprint::PtrType>::iterator i = footprintListIn.begin(); i != footprintListIn.end(); ++i) {
            /* footprint has not enough pixels */
            if (static_cast<unsigned int>((*i)->getNpix()) < footprintDiffimNpixMin) {
                continue;
            } 

            /* grab the BBox and grow it; this will eventually be overridden by a grow method on Footprint */
            BBox2i footprintBBox = (*i)->getBBox();
            footprintBBox.grow(footprintBBox.max() + vw::Vector2i(footprintDiffimGrow,footprintDiffimGrow));
            footprintBBox.grow(footprintBBox.min() - vw::Vector2i(footprintDiffimGrow,footprintDiffimGrow));

            /* grab a subimage; there is an exception if its e.g. too close to the image */
            try {
                imageToConvolveFootprintPtr = imageToConvolve.getSubImage(footprintBBox);
                imageToNotConvolveFootprintPtr = imageToNotConvolve.getSubImage(footprintBBox);
            } catch (lsst::mwi::exceptions::ExceptionStack &e) {
                continue;
            }
            
            if (lsst::imageproc::maskOk(*(imageToConvolveFootprintPtr->getMask()), badPixelMask) && 
                lsst::imageproc::maskOk(*(imageToNotConvolveFootprintPtr->getMask()), badPixelMask) ) {

                /* Create a new footprint with grow'd box */
                lsst::detection::Footprint::PtrType fpGrow(new lsst::detection::Footprint(footprintBBox));
                footprintListOut.push_back(fpGrow);
                
                nCleanFootprints += 1;
            }
        }
        
        footprintDetectionThreshold *= detectionThresholdScaling;
    }
    lsst::mwi::utils::Trace("lsst.imageproc.getCollectionOfFootprintsForPsfMatching", 3, 
                            (boost::format("Found %d clean footprints above threshold %.3f") 
                             % footprintListOut.size() % (footprintDetectionThreshold/detectionThresholdScaling)));

    return footprintListOut;
}

/**
 * \brief Test function to get a list of footprints relevant to a particular test image
 *
 * \return a list of footprints suitable for psf matching
 */
std::vector<lsst::detection::Footprint::PtrType> lsst::imageproc::getCollectionOfMaskedImagesForPsfMatching(
    ) {
    
    double radius = 20.0;
    std::vector<lsst::detection::Footprint::PtrType> footprintList;
    
    
    // SuperMACHO
    // radius = 41;
    // vw::BBox2i const region0(static_cast<int>(2119.-radius/2),
    //                          static_cast<int>(581.-radius/2),
    //                          static_cast<int>(radius),
    //                          static_cast<int>(radius));
    // lsst::detection::Footprint::PtrType fp0(
    //     new lsst::detection::Footprint(region0)
    //     );
    // footprintList.push_back(fp0);
    // 
    // return;
    
    vw::BBox2i const region1(static_cast<int>(78.654-radius/2),
                             static_cast<int>(3573.945-radius/2),
                             static_cast<int>(radius),
                             static_cast<int>(radius));
    lsst::detection::Footprint::PtrType fp1(
        new lsst::detection::Footprint(region1)
        );
    footprintList.push_back(fp1);
    
    
    vw::BBox2i const region4(static_cast<int>(367.756-radius/2),
                             static_cast<int>(3827.671-radius/2),
                             static_cast<int>(radius),
                             static_cast<int>(radius));
    lsst::detection::Footprint::PtrType fp4(
        new lsst::detection::Footprint(region4)
        );
    footprintList.push_back(fp4);
    
    vw::BBox2i const region5(static_cast<int>(381.062-radius/2),
                             static_cast<int>(3212.948-radius/2),
                             static_cast<int>(radius),
                             static_cast<int>(radius));
    lsst::detection::Footprint::PtrType fp5(
        new lsst::detection::Footprint(region5)
        );
    footprintList.push_back(fp5);
    
    vw::BBox2i const region6(static_cast<int>(404.433-radius/2),
                             static_cast<int>(573.462-radius/2),
                             static_cast<int>(radius),
                             static_cast<int>(radius));
    lsst::detection::Footprint::PtrType fp6(
        new lsst::detection::Footprint(region6)
        );
    footprintList.push_back(fp6);
    
    vw::BBox2i const region7(static_cast<int>(420.967-radius/2),
                             static_cast<int>(3306.310-radius/2),
                             static_cast<int>(radius),
                             static_cast<int>(radius));
    lsst::detection::Footprint::PtrType fp7(
        new lsst::detection::Footprint(region7)
        );
    footprintList.push_back(fp7);
    
    /* cosmic ray */
    vw::BBox2i const region9(static_cast<int>(546.657-radius/2),
                             static_cast<int>(285.079-radius/2),
                             static_cast<int>(radius),
                             static_cast<int>(radius));
    lsst::detection::Footprint::PtrType fp9(
        new lsst::detection::Footprint(region9)
        );
    footprintList.push_back(fp9);
    
    // OR
    //lsst::detection::Footprint *fp4 = new lsst::detection::Footprint(1, region4);
    //lsst::detection::Footprint::PtrType fpp4(fp4); 
    //footprintList.push_back(fpp4);
    
    return footprintList;
}

/**
 * \brief Compute a set of Principal Component Analysis (PCA) basis kernels
 *
 * \return A list of PCA basis kernels
 *
 * \throw lsst::mwi::exceptions::Runtime if no good kernels
 */
template <typename KernelT>
std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > lsst::imageproc::computePcaKernelBasis(
    std::vector<lsst::imageproc::DiffImContainer<KernelT> > &diffImContainerList, ///< List of input footprints
    lsst::mwi::policy::Policy &policy ///< Policy directing the behavior
    ) {
    std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > kernelPcaBasisList;
    
    // Parse the Policy
    unsigned int minimumNumberOfBases = policy.getInt("computePcaKernelBasis.minimumNumberOfBases");
    unsigned int maximumNumberOfBases = policy.getInt("computePcaKernelBasis.maximumNumberOfBases");
    double maximumFractionOfEigenvalues = policy.getDouble("computePcaKernelBasis.maximumFractionOfEigenvalues");
    double minimumAcceptibleEigenvalue = policy.getDouble("computePcaKernelBasis.minimumAcceptibleEigenvalue");
    unsigned int maximumIteratonsPCA = policy.getInt("computePcaKernelBasis.maximumIteratonsPCA");
    double maximumKernelResidualMean = policy.getDouble("computePcaKernelBasis.maximumKernelResidualMean");
    double maximumKernelResidualVariance = policy.getDouble("computePcaKernelBasis.maximumKernelResidualVariance");
    bool debugIO = policy.getBool("debugIO", false);
    
    // Image accessor
    typedef typename vw::ImageView<KernelT>::pixel_accessor imageAccessorType;
    // Iterator over struct
    typedef typename std::vector<lsst::imageproc::DiffImContainer<KernelT> >::iterator iDiffImContainer;
    
    // Matrix to invert.  Number of rows = number of pixels; number of columns = number of kernels
    // All calculations here are in double
    // Assigment of matrix sizes is rows, cols
    vw::math::Matrix<double> M;
    vw::math::Matrix<double> eVec;
    vw::math::Vector<double> eVal;
    vw::math::Vector<double> mMean;
    vw::math::Matrix<double> kernelCoefficientMatrix;
    
    // Initialize for while loop
    unsigned int nIter = 0;
    unsigned int nReject = 1;
    unsigned int nGood = 0;
    unsigned int nKCols = 0, nKRows = 0;
    unsigned int nPixels;
    unsigned int nCoeffToUse;
    lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 3, " ");
    lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 3, 
                            "Entering subroutine computePcaKernelBasis");
    
    /* Reuse */
    /* NOTE - this assumes that all kernels are the same size */
    /*        and that this size is defined in the policy */ 
    unsigned int kernelRows = policy.getInt("kernelRows");
    unsigned int kernelCols = policy.getInt("kernelCols");
    lsst::fw::Image<KernelT> kImage(kernelCols, kernelRows);
    
    // Iterate over PCA inputs until all are good
    while ( (nIter < maximumIteratonsPCA) and (nReject != 0) ) {
        
        nGood = 0;
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
            if ((*i).isGood == true) {
                nGood += 1;
                if (nKCols == 0) {
                    nKCols = (*i).diffImKernelPtr->getCols();
                }
                if (nKRows == 0) {
                    nKRows = (*i).diffImKernelPtr->getRows();
                }
            }
        }
        
        if (nGood == 0) {
            throw lsst::mwi::exceptions::DomainError("No good kernels for PCA");
        }
        else {
            lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 4, 
                                    (boost::format("PCA Iteration %d : Using %d kernels") % nIter % nGood));
        }

        /* nGood can only decrement; make sure you don't ask for more bases than you have input kernels */
        if (nGood < minimumNumberOfBases) {
            minimumNumberOfBases = nGood;
            maximumNumberOfBases = nGood;
        }
        if (nGood < maximumNumberOfBases) {
            maximumNumberOfBases = nGood;
        }
        
        

        
        nPixels = nKCols * nKRows;
        M.set_size(nPixels, nGood);
        eVec.set_size(nPixels, nGood);
        eVal.set_size(nGood);
        mMean.set_size(nPixels);
        
        // fill up matrix for PCA
        double imSum = 0.0;
        int ki = 0;
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
            if ((*i).isGood == false) {
                continue;
            }
            
            (*i).diffImKernelPtr->computeImage(kImage, imSum, 0.0, 0.0, false);
            /* insanity checking */
            assert(nKRows == kImage.getRows());
            assert(nKCols == kImage.getCols());
            
            int mIdx = 0;
            imageAccessorType imageAccessorCol(kImage.origin());
            for (unsigned int col = 0; col < nKCols; ++col, imageAccessorCol.next_col()) {
                
                imageAccessorType imageAccessorRow(imageAccessorCol);
                for (unsigned int row = 0; row < nKRows; ++row, ++mIdx, imageAccessorRow.next_row()) {
                    
                    // NOTE : 
                    //   arguments to matrix-related functions are given in row,col
                    //   arguments to image-related functions are given in col,row
                    // HOWEVER :
                    //   it doesn't matter which order we put the kernel elements into the PCA
                    //   since they are not correlated 
                    //   as long as we do the same thing when extracting the components
                    // UNLESS :
                    //   we want to put in some weighting/regularlization into the PCA
                    //   not sure if that is even possible...

                    M(mIdx, ki) = *imageAccessorRow;
                }
            }
            
            /* only increment when you have a good footprint */
            ki += 1;
        }
        
       
        // M is mean-subtracted if subtractMean = true
        // Eigencomponents in columns of eVec
        lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 4, "Computing pricipal components");
        lsst::imageproc::computePca(mMean, eVal, eVec, M, true);
        lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 4, "Computed pricipal components");

        /* This is probably as inefficient as you can get, but I'm no boost expert */
        boost::format eValueFormatter("Eigenvalues :");
        for (unsigned int i = 0; i < eVal.size(); ++i) {
            eValueFormatter = boost::format("%s %.3f") % eValueFormatter.str() % eVal[i];
        }
        lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 4, eValueFormatter);
        
        nCoeffToUse = minimumNumberOfBases;
        // Potentially override with larger number if the spectrum of eigenvalues requires it
        double evalSum = 0.0;
        for (unsigned int i = 0; i < eVal.size(); ++i) {
            evalSum += eVal[i];
        }
        double evalFrac = 0.0;
        for (unsigned int i = 0; i < nCoeffToUse; ++i) {
            evalFrac += eVal[i] / evalSum;
            if (eVal[i] < minimumAcceptibleEigenvalue) {
                lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 1, 
                                        (boost::format("WARNING : Using eigenvector whose eigenvalue (%.3e) is smaller than acceptible (%.3e)")
                                         % eVal[i] % minimumAcceptibleEigenvalue));
            }
            
        }
        if (evalFrac < maximumFractionOfEigenvalues) {
            for (; ( (nCoeffToUse < eVal.size()) && (nCoeffToUse < maximumNumberOfBases) ); ++nCoeffToUse) {
                evalFrac += eVal[nCoeffToUse] / evalSum;
                if (eVal[nCoeffToUse] < minimumAcceptibleEigenvalue) {
                    lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 1, 
                                            (boost::format("WARNING : Using eigenvector whose eigenvalue (%.3e) is smaller than acceptible (%.3e)")
                                             % eVal[nCoeffToUse] % minimumAcceptibleEigenvalue));
                }
                if (evalFrac > maximumFractionOfEigenvalues) {
                    break;
                }
            }
        }
        lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 4, 
                                (boost::format("Using %d basis functions (plus mean) representing %.5f of the variance") 
                                 % nCoeffToUse % evalFrac));
        
        // We now have the basis functions determined
        // Next determine the coefficients that go in front of all of the individual kernels

        /* 
           Note on matrix sizes; row x col 
           M    = nPix  x nGood
           eVec = nPix  x nGood
           coeff= nGood x nToUse

           inside decomposeMatrixUsingBasis :
              grab each column of M, representing each individual kernel (length nPix)
              then grab each column of eVec, representing each basis (length nPix)
              take their dot product, stuffing the results into coeff (length 1)
              do the latter only up to nToUse
              thus coeff is number of individual kernels by number of bases to use
        */
        kernelCoefficientMatrix.set_size(nGood, nCoeffToUse);
        lsst::imageproc::decomposeMatrixUsingBasis(kernelCoefficientMatrix, M, eVec, nCoeffToUse);
        
        // We next do quality control here; we reconstruct the input kernels with the truncated basis function set
        // Remember that M is mean-subtracted already
        /* 
           Note on matrix sizes; row x col
           approxM = nPix x nGood (same as M)
           eVec    = nPix x nGood
           coeff   = nGood x nToUse

           inside approximateMatrixUsingBasis :
              create a temporary vector of length nPix
              for each column of approxM i, representing a reconstructed kernel (length nPix)
                 set temporary vector equal to zero
                 for each column of eVec j, representing each basis up to nToUse (length nPix)
                 grab coefficient i,j and multiply by eVec j and add it to temporary vector
                 put the result in column of approxM
        */
        vw::math::Matrix<double> approxM(nPixels, nGood); 
        lsst::imageproc::approximateMatrixUsingBasis(approxM, eVec, kernelCoefficientMatrix, nCoeffToUse);
        
#if defined(DEBUG_MATRIX)
        std::cout << "EigenKernel coefficients : " << kernelCoefficientMatrix << std::endl;
#endif
        
        nReject = 0;
        unsigned int iKernel = 0;
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
            if ((*i).isGood == false) {
                continue;
            }

            /* NOTE TO SELF - DO WE LOOK AT THE KERNEL, OR THE DIFFERENCE IMAGE HERE */
            /* THE DIFFERENCE IMAGE IS EASIER TO ANALYZE */

            double kernelResidual;
            double kernelResidualVariance;
            
            vw::math::Vector<double> kernelDifference = (vw::math::select_col(M, iKernel) - vw::math::select_col(approxM, iKernel));
            calculateVectorStatistics(kernelDifference, kernelResidual, kernelResidualVariance);
            
            (*i).kernelResidual = kernelResidual;
            (*i).kernelResidualVariance = kernelResidualVariance;
            
            if (!(std::isfinite((*i).kernelResidual)) ||
                !(std::isfinite((*i).kernelResidualVariance))) {
                lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 2, 
                                        (boost::format("# Kernel %d, PCA analysis failed : %.3f %.3f") 
                                         % (*i).id % (*i).kernelResidual % (*i).kernelResidualVariance));
                (*i).isGood = false;
                nReject += 1;
            }
            else if ((*i).kernelResidual > maximumKernelResidualMean) {
                lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 2, 
                                        (boost::format("# Kernel %d, bad mean residual of PCA model : %.3f") 
                                         % (*i).id % (*i).kernelResidual));
                (*i).isGood = false;
                nReject += 1;
            }
            else if ((*i).kernelResidualVariance > maximumKernelResidualVariance) {
                lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 2, 
                                        (boost::format("# Kernel %d, bad residual variance of PCA model : %.3f") 
                                         % (*i).id % (*i).kernelResidualVariance));
                (*i).isGood = false;
                nReject += 1;
            }
            else {
                lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 5, 
                                        (boost::format("Kernel %d, mean and variance of residual from PCA decomposition : %10.3e %10.3e") 
                                         % (*i).id % (*i).kernelResidual % (*i).kernelResidualVariance));
            }
            /* only for the good footprints */
            iKernel += 1;
        }

        nIter += 1;
    } // End of iteration
    
    // Turn the Mean Image into a Kernel
    lsst::fw::Image<KernelT> meanImage(nKCols, nKRows);
    imageAccessorType imageAccessorCol(meanImage.origin());
    int mIdx = 0;
    for (unsigned int col = 0; col < nKCols; ++col) {
        imageAccessorType imageAccessorRow(imageAccessorCol);
        for (unsigned int row = 0; row < nKRows; ++row, ++mIdx) {
            *imageAccessorRow = mMean(mIdx);
            imageAccessorRow.next_row();
        }
        imageAccessorCol.next_col();
    }
    // The mean image is the first member of kernelPCABasisList
    kernelPcaBasisList.push_back(boost::shared_ptr<lsst::fw::Kernel<KernelT> > 
                                 (new lsst::fw::FixedKernel<KernelT>(meanImage)));
    
    if (debugIO) {
        meanImage.writeFits( (boost::format("mFits.fits")).str() );
    }
    
    // Turn each eVec into an Image and then into a Kernel
    // Dont use all eVec.cols(); only the number of bases that you want
    for (unsigned int ki = 0; ki < nCoeffToUse; ++ki) { 
        lsst::fw::Image<KernelT> basisImage(nKCols, nKRows);
        
        // Not sure how to bulk load information into Image
        int kIdx = 0;
        
        imageAccessorType imageAccessorCol(basisImage.origin());
        for (unsigned int col = 0; col < nKCols; ++col) {
            
            imageAccessorType imageAccessorRow(imageAccessorCol);
            for (unsigned int row = 0; row < nKRows; ++row, ++kIdx) {
                
                *imageAccessorRow = eVec(kIdx, ki);
                imageAccessorRow.next_row();
            }
            imageAccessorCol.next_col();
        }
        // Add to kernel basis
        kernelPcaBasisList.push_back(boost::shared_ptr<lsst::fw::Kernel<KernelT> > 
                                     (new lsst::fw::FixedKernel<KernelT>(basisImage)));
        
        if (debugIO) {
            basisImage.writeFits( (boost::format("eFits_%d.fits") % ki).str() );
        }
    }
    
    // Finally, create a new LinearCombinationKernel for each Footprint
    unsigned int iKernel = 0;
    for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
        if ((*i).isGood == false) {
            continue;
        }

        std::vector<double> kernelCoefficients;
        
        // Mean image
        kernelCoefficients.push_back(1); 
        for (unsigned int jKernel = 0; jKernel < nCoeffToUse; ++jKernel) {
            kernelCoefficients.push_back(kernelCoefficientMatrix(iKernel,jKernel));
        }        
        
        // Create a linear combination kernel from this 
        boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > pcaKernelPtr(
            new lsst::fw::LinearCombinationKernel<KernelT>(kernelPcaBasisList, kernelCoefficients)
            );
        
        (*i).diffImPcaKernelPtr = pcaKernelPtr;
        
        /* only for the good footprints */
        iKernel += 1;
    }
    return kernelPcaBasisList;
}

/**
 * \brief Compute a spatially varying PSF-matching kernel
 *
 * Fit the spatial variation of a spatially varying PSF-matching kernel
 * and the spatial variation of the background given a set of basis kernels
 * and an associated set of information about them.
 *
 * \return spatially varying PSF-matching kernel
 */
template <typename KernelT>
boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > lsst::imageproc::computeSpatiallyVaryingPsfMatchingKernel(
    boost::shared_ptr<lsst::fw::function::Function2<double> > &kernelFunctionPtr, ///< Function for spatial variation of kernel
    boost::shared_ptr<lsst::fw::function::Function2<double> > &backgroundFunctionPtr, ///< Function for spatial variation of background
    std::vector<lsst::imageproc::DiffImContainer<KernelT> > &diffImContainerList, ///< Information on each kernel
    std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelBasisList, ///< Basis kernel set
    lsst::mwi::policy::Policy &policy ///< Policy directing the behavior
    )
{
    // Parse the Policy
    unsigned int maximumIterationsSpatialFit = policy.getInt("computeSpatiallyVaryingPsfMatchingKernel.maximumIterationsSpatialFit");
    double maximumSpatialKernelResidualMean = policy.getDouble("computeSpatiallyVaryingPsfMatchingKernel.maximumSpatialKernelResidualMean");
    double maximumSpatialKernelResidualVariance = policy.getDouble("computeSpatiallyVaryingPsfMatchingKernel.maximumSpatialKernelResidualVariance");
    bool debugIO = policy.getBool("debugIO", false);
    
    // Container iterator
    typedef typename std::vector<lsst::imageproc::DiffImContainer<KernelT> >::iterator iDiffImContainer;
    // Kernel iterator
    typedef typename std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > >::const_iterator iKernelPtr;
    
    // Initialize for while loop
    unsigned int nIter = 0;
    unsigned int nReject = 1;
    unsigned int nGood = diffImContainerList.size();
    unsigned int nKernel = kernelBasisList.size();
    unsigned int nKernelParameters = kernelFunctionPtr->getNParameters();
    unsigned int nBgParameters = backgroundFunctionPtr->getNParameters();
    
    lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPsfMatchingKernel", 3, 
                            "Entering subroutine computeSpatiallyVaryingPsfMatchingKernel");

    /* for the kernel fit */
    std::vector<double> kernelMeas;
    std::vector<double> kernelVariances;
    std::vector<double> kernelColPos;
    std::vector<double> kernelRowPos;
    double nSigmaSq = 1.0;

    /* for the background fit */
    std::vector<double> backgroundMeas;
    std::vector<double> backgroundVariances;

    double imSum = 0.0;
    double bgSum = 0.0;
    
    /* NOTE - this assumes that all kernels are the same size */
    /*        and that this size is defined in the policy */ 
    unsigned int kernelRows = policy.getInt("kernelRows");
    unsigned int kernelCols = policy.getInt("kernelCols");
    lsst::fw::Image<KernelT> kImage(kernelCols, kernelRows);
    lsst::fw::Image<KernelT> kSpatialImage(kernelCols, kernelRows);
    
    // Set up the spatially varying kernel
    boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > spatiallyVaryingKernelPtr(
        new lsst::fw::LinearCombinationKernel<KernelT>(kernelBasisList, kernelFunctionPtr));
    
    /* Iterate over inputs until both the kernel model and the background model converge */
    while ( (nIter < maximumIterationsSpatialFit) and (nReject != 0) ) {
        /* Start out clean, nothing rejected yet */
        nReject = 0;

        /* Set up fit to kernel */
        kernelColPos.clear();
        kernelRowPos.clear();
        nGood = 0;
        bgSum = 0.0;
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
            /* those attributs that you grab once for each container; position and background */
            if ((*i).isGood == true) {
                kernelColPos.push_back((*i).colcNorm);
                kernelRowPos.push_back((*i).rowcNorm);

                bgSum += (*i).background;
                backgroundMeas.push_back((*i).background);
                backgroundVariances.push_back((*i).footprintResidualVariance);   /* approximation */

                lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForMaskedImage", 5, 
                                        (boost::format("Background %d at %.3f %.3f = %.3f") 
                                         % (*i).id % (*i).colcNorm % (*i).rowcNorm % (*i).background));

                nGood += 1;
            }
        }
        
        if (nGood == 0) {
            throw lsst::mwi::exceptions::DomainError("No good footprints for spatial kernel and background fitting");
        }
        
        /* ******************************** */
        /* ******************************** */
        /* Fit for spatially varying kernel */

        std::vector<std::vector<double> > fitParameters = spatiallyVaryingKernelPtr->getSpatialParameters();
        assert(fitParameters.size() == nKernel);
        assert(fitParameters[0].size() == nKernelParameters);
        
        /* The first kernel is the mean, which is constant across the image */
        fitParameters[0][0] = 1.0;
        for (unsigned int i = 1; i < nKernelParameters; ++i) {
            fitParameters[0][i] = 0.0;
        }

        /* For each kernel, grab each footprint's constraint on it, and fit for the spatial variation of those constraints */
        unsigned int nk = 1;
        for (iKernelPtr ik = (kernelBasisList.begin()+1); ik != kernelBasisList.end(); ++nk, ++ik) {
            
            boost::format spatialFitFormatter("Measurements of kernel PC%d :");
            spatialFitFormatter % nk;
            
            kernelMeas.clear();
            kernelVariances.clear();
            for (iDiffImContainer id = diffImContainerList.begin(); id != diffImContainerList.end(); ++id) {
                if ((*id).isGood == true) {
                    /* inefficient; would be nice to have a method to grab only parameter i */
                    std::vector<double> kernelParameters = (*id).diffImPcaKernelPtr->getKernelParameters(); 
                    kernelMeas.push_back(kernelParameters[nk]);

                    kernelVariances.push_back((*id).footprintResidualVariance); /* approximation */
                    spatialFitFormatter = boost::format("%s %.3f") % spatialFitFormatter.str() % kernelParameters[nk];
                }
            }
            lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPsfMatchingKernel", 6, spatialFitFormatter);
            
            /* NOTE - if we have fewer measurements than kernelFunctionPtr->getNParameters(), we should do something about it */
            std::vector<double> kernelParams(nKernelParameters);
            std::vector<double> kernelStepsize(nKernelParameters);
            std::fill(kernelParams.begin(), kernelParams.end(), 0.0);         /* start at value zero */
            std::fill(kernelStepsize.begin(), kernelStepsize.end(), 0.02);  /* assume order 2% contribution from eigenkernels */
            
            /* Run minimization */
            lsst::fw::function::FitResults kernelFit = lsst::fw::function::minimize(
                kernelFunctionPtr,
                kernelParams,
                kernelStepsize,
                kernelMeas,
                kernelVariances,
                kernelColPos,
                kernelRowPos,
                nSigmaSq
                );

            /* Results of spatial fit for this particular kernel */
            for (unsigned int i = 0; i < kernelFit.parameterList.size(); ++i) {
                fitParameters[nk][i] = kernelFit.parameterList[i];
                lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPsfMatchingKernel", 4, 
                                        (boost::format("Fit to kernel PC%d, spatial parameter %d : %.3f (%.3f,%.3f)\n") 
                                         % nk
                                         % i
                                         % kernelFit.parameterList[i]
                                         % kernelFit.parameterErrorList[i].first
                                         % kernelFit.parameterErrorList[i].second));
            }
        }
        /* After all kernels are fit, fill the spatially varying kernel parameters */
        spatiallyVaryingKernelPtr->setSpatialParameters(fitParameters);

        /* Check each footprint's actual kernel with the spatial model */
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
            if ((*i).isGood == false) {
                continue;
            }
            
            /* calculate spatial kernel at this location */
            spatiallyVaryingKernelPtr->computeImage(kSpatialImage, imSum, 
                                                    (*i).colcNorm, (*i).rowcNorm,
                                                    false);
            
            /* compare to the original kernel */
            (*i).diffImKernelPtr->computeImage(kImage, imSum, 
                                               0.0, 0.0, 
                                               false);
            
            kSpatialImage -= kImage;
            if (debugIO) {
                kSpatialImage.writeFits( (boost::format("ksmFits_%d.fits") % (*i).id).str() );
            }
            
            /* NOTE TO SELF - DO WE LOOK AT THE KERNEL, OR THE DIFFERENCE IMAGE HERE */
            /* THE DIFFERENCE IMAGE IS EASIER TO ANALYZE */

            /* How good does the spatial model match the kernel */
            double meanOfResiduals = 0.0;
            double varianceOfResiduals = 0.0;
            int nGoodPixels = 0;
            calculateImageResiduals(nGoodPixels, meanOfResiduals, varianceOfResiduals, kSpatialImage);
            (*i).spatialKernelResidual = meanOfResiduals;
            (*i).spatialKernelResidualVariance = varianceOfResiduals/nGoodPixels; // per pixel
            
            if (!(std::isfinite((*i).spatialKernelResidual)) || 
                !(std::isfinite((*i).spatialKernelResidualVariance))) {
                lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPsfMatchingKernel", 2, 
                                        (boost::format("# Kernel %d, spatial analysis failed : %.3f %.3f") 
                                         % (*i).id % (*i).spatialKernelResidual % (*i).spatialKernelResidualVariance));
                (*i).isGood = false;
                nReject += 1;
            }                    
            else if ((*i).spatialKernelResidual > maximumSpatialKernelResidualMean) {
                lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPsfMatchingKernel", 2, 
                                        (boost::format("# Kernel %d, bad mean residual of spatial fit : %.3f") 
                                         % (*i).id % (*i).spatialKernelResidual));
                (*i).isGood = false;
                nReject += 1;
            }
            else if ((*i).spatialKernelResidualVariance > maximumSpatialKernelResidualVariance) {
                lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPsfMatchingKernel", 2, 
                                        (boost::format("# Kernel %d, bad residual variance of spatial fit : %.3f") 
                                         % (*i).id % (*i).spatialKernelResidualVariance));
                (*i).isGood = false;
                nReject += 1;
            }
            else {
                lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPsfMatchingKernel", 5, 
                                        (boost::format("Kernel %d, mean and variance of residual from spatial fit : %10.3e %10.3e") 
                                         % (*i).id % (*i).spatialKernelResidual % (*i).spatialKernelResidualVariance));
            }
            
            /* And importantly, how does the spatial kernel work for diffim */
            /* TODO; part of the problem is that the images are not in DiffImContainer */
        }

        /* Fit for spatially varying kernel */
        /* ******************************** */
        /* ******************************** */
        /* Fit for spatially varying background */
        
        /* Initialize fit parameters at average background; higher order terms initially zero */
        std::vector<double> backgroundParameters(nBgParameters);
        std::vector<double> backgroundStepsize(nBgParameters);
        std::fill(backgroundParameters.begin(), backgroundParameters.end(), 0.0);
        std::fill(backgroundStepsize.begin(), backgroundStepsize.end(), 1);
        /* Initialize constant term at average, with 10% step size */
        backgroundParameters[0] = bgSum / nGood; 
        backgroundStepsize[0] = 0.1 * backgroundParameters[0];
        
        /* Run minimization */
        lsst::fw::function::FitResults backgroundFit = lsst::fw::function::minimize(
            backgroundFunctionPtr,
            backgroundParameters,
            backgroundStepsize,
            backgroundMeas,
            backgroundVariances,
            kernelColPos,
            kernelRowPos,
            nSigmaSq
            );
        backgroundFunctionPtr->setParameters(backgroundFit.parameterList);
        
        /* Debugging information */
        for (unsigned int i = 0; i < backgroundFit.parameterList.size(); ++i) {
            lsst::mwi::utils::Trace("lsst.imageproc.computePsfMatchingKernelForMaskedImage", 4, 
                                    (boost::format("Fit Background Parameter %d : %.3f (%.3f,%.3f)\n") 
                                     % i
                                     % backgroundFit.parameterList[i]
                                     % backgroundFit.parameterErrorList[i].first
                                     % backgroundFit.parameterErrorList[i].second));
        }

        /* Check each footprint's actual background with the spatial model */
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
            if ((*i).isGood == true) {
                /* calculate spatial background at this location */
                bgSum = (*backgroundFunctionPtr)((*i).colcNorm, (*i).rowcNorm);
                /* NOTE - ACTUALLY COMPARE THE MODEL TO THE DATA! */
                if (0) {
                    (*i).isGood = false;
                    nReject += 1;
                }
            }
        }

        /* Fit for spatially varying background */
        /* ******************************** */
        /* ******************************** */

        /* End of the loop */
        nIter += 1;
    }
    return spatiallyVaryingKernelPtr;
}

/**
 * \brief Compute a set of delta function basis kernels
 *
 * \return a set of delta function basis kernels
 * \throw lsst::mwi::exceptions::DomainError if nRows or nCols not positive
 */
template <typename KernelT>
std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > lsst::imageproc::generateDeltaFunctionKernelSet(
    unsigned int nCols, ///< Number of colunms in the basis kernels
    unsigned int nRows  ///< Number of rows in the basis kernels
    )
{
    if ((nCols < 1) || (nRows < 1)) {
        throw lsst::mwi::exceptions::DomainError("nRows and nCols must be positive");
    }
    const int signedNCols = static_cast<int>(nCols);
    const int signedNRows = static_cast<int>(nRows);
    std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > kernelBasisList;
    int colCtr = (signedNCols - 1) / 2;
    int rowCtr = (signedNRows - 1) / 2;
    for (int row = 0; row < signedNRows; ++row) {
        double y = static_cast<double>(row - rowCtr);
        
        for (int col = 0; col < signedNCols; ++col) {
            double x = static_cast<double>(col - colCtr);
            
            typename lsst::fw::Kernel<KernelT>::KernelFunctionPtrType kfuncPtr(
                new lsst::fw::function::IntegerDeltaFunction2<KernelT>(x, y)
                );
            
            boost::shared_ptr<lsst::fw::Kernel<KernelT> > kernelPtr(
                new lsst::fw::AnalyticKernel<KernelT>(kfuncPtr, nCols, nRows)
                );
            
            kernelBasisList.push_back(kernelPtr);
        }
    }
    return kernelBasisList;
}

/**
 * \brief Compute a set of Alard-Lupton basis kernels
 *
 * WARNING: Not yet implemented! Always throws DomainError!
 *
 * \return a set of Alard-Lupton basis kernels
 * \throw lsst::mwi::exceptions::DomainError if nRows or nCols not positive
 */
template <typename KernelT>
std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > lsst::imageproc::generateAlardLuptonKernelSet(
    unsigned int nRows, ///< Number of rows in the kernel basis
    unsigned int nCols, ///< Number of columns in the kernel basis
    std::vector<double> const &sigGauss, ///< Width of gaussians in basis; size = number of Gaussians
    std::vector<double> const &degGauss  ///< Degree of spatial variation within each Gaussian; size = sigGauss.size()
    )
{
    if ((nCols < 1) || (nRows < 1)) {
        throw lsst::mwi::exceptions::DomainError("nRows and nCols must be positive");
    }
    throw lsst::mwi::exceptions::DomainError("Not implemented");
    std::vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > kernelBasisList;
    return kernelBasisList;
}

/**
 * \brief Check a MaskedImage for the presence of bad pixels
 *
 * \return false if any bad pixels found, true otherwise
 */
template <typename MaskT>
bool lsst::imageproc::maskOk(
    lsst::fw::Mask<MaskT> const &inputMask,
    MaskT const badPixelMask ///< Mask value for bad data
    )
{
    typename lsst::fw::Mask<MaskT>::pixel_accessor rowAcc = inputMask.origin();
    for (unsigned int row = 0; row < inputMask.getRows(); ++row, rowAcc.next_row()) {
        typename lsst::fw::Mask<MaskT>::pixel_accessor colAcc = rowAcc;
        for (unsigned int col = 0; col < inputMask.getCols(); ++col, colAcc.next_col()) {
            //std::cout << "MASK " << (*colAcc) << " " << badPixelMask << " " << ((*colAcc) & badPixelMask) << std::endl;
            
            if (((*colAcc) & badPixelMask) != 0) {
                return false;
            }
        }
    }
    return true;
}

template <typename ImageT, typename MaskT>
void lsst::imageproc::calculateMaskedImageResiduals(
    int &nGoodPixels, ///< Number of good pixels in the image
    double &meanOfResiduals, ///< Mean residual/variance; ideally 0
    double &varianceOfResiduals, ///< Average variance of residual/variance; ideally 1
    lsst::fw::MaskedImage<ImageT, MaskT> const &inputImage, ///< Input image to be analyzed
    MaskT const badPixelMask ///< Mask for bad data
    ) {
    
    // BASICALLY A HACK FOR DC2 TESTING; SENDING THE SAME IMAGE TO DIFFIM RESULTS IN 0 VARIANCE; MESSING UP FITTING
    // Add in quadrature as noise
    const double MinVariance = 1e-20;
    
    double x2Sum=0.0, xSum=0.0, wSum=0.0;
    
//    bool didReport = false;
    nGoodPixels = 0;
    lsst::fw::MaskedPixelAccessor<ImageT, MaskT> rowAcc(inputImage);
    for (unsigned int row = 0; row < inputImage.getRows(); ++row, rowAcc.nextRow()) {
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> colAcc = rowAcc;
        for (unsigned int col = 0; col < inputImage.getCols(); ++col, colAcc.nextCol()) {
            // std::cout << "MASK " << (*rowAcc.mask) << " " << badPixelMask << " " << ((*rowAcc.mask) & badPixelMask) << std::endl;
            if (((*colAcc.mask) & badPixelMask) == 0) {
                double var = std::max(static_cast<double>(*colAcc.variance), MinVariance);
                nGoodPixels += 1;
                x2Sum += (*colAcc.image) * (*colAcc.image) / var;
                xSum  += (*colAcc.image) / var;
                wSum  += 1.0 / var;
//                if (!didReport && (std::isnan(x2Sum) || std::isnan(xSum) || std::isnan(wSum))) {
//                    std::cout << "x2Sum=" << x2Sum << "; xSum = " << xSum << "; wSum = " << wSum
//                        << "; *colAcc.image = " << *colAcc.image << "; *colAcc.variance = " << *colAcc.variance
//                        << "; var = " << var
//                        << "; row = " << row << "; col = " << col << std::endl;
//                    didReport = true;
//                }
            }
        }
    }
    
    if (nGoodPixels > 0) {
        meanOfResiduals = xSum / wSum;
    } else {
        meanOfResiduals = std::numeric_limits<double>::quiet_NaN();
    }
    if (nGoodPixels > 1) {
        varianceOfResiduals  = x2Sum / wSum - meanOfResiduals*meanOfResiduals;
        varianceOfResiduals *= nGoodPixels / (nGoodPixels - 1);
        varianceOfResiduals += MinVariance;
    } else {
        varianceOfResiduals = std::numeric_limits<double>::quiet_NaN();
    }
    
}

template <typename ImageT>
void lsst::imageproc::calculateImageResiduals(
    int &nGoodPixels, ///< Number of good pixels in the image
    double &meanOfResiduals, ///< Mean residual; nan if nGoodPixels < 1
    double &varianceOfResiduals, ///< Average variance of residual; nan if nGoodPixels < 2
    lsst::fw::Image<ImageT> const &inputImage ///< Input image to be analyzed
    ) {
    
    double x2Sum=0.0, xSum=0.0, wSum=0.0;
    
    nGoodPixels = 0;
    typedef typename vw::ImageView<ImageT>::pixel_accessor imageAccessorType;
    
    imageAccessorType imageAccessorCol(inputImage.origin());
    for (unsigned int col = 0; col < inputImage.getCols(); ++col) {
        
        imageAccessorType imageAccessorRow(imageAccessorCol);
        for (unsigned int row = 0; row < inputImage.getRows(); ++row) {
            nGoodPixels += 1;
            x2Sum       += (*imageAccessorRow) * (*imageAccessorRow);
            xSum        += (*imageAccessorRow);
            wSum        += 1;
            imageAccessorRow.next_row();
        }
        imageAccessorCol.next_col();
    }
    
    if (nGoodPixels > 0) {
        meanOfResiduals = xSum / wSum;
    } else {
        meanOfResiduals = std::numeric_limits<double>::quiet_NaN();
    }
    if (nGoodPixels > 1) {
        varianceOfResiduals  = x2Sum / wSum - meanOfResiduals*meanOfResiduals;
        varianceOfResiduals *= nGoodPixels / (nGoodPixels - 1);
    } else {
        varianceOfResiduals = std::numeric_limits<double>::quiet_NaN();
    }
}

template <typename VectorT>
void lsst::imageproc::calculateVectorStatistics(
    vw::math::Vector<VectorT> const &inputVector,
    double &mean,
    double &variance
    ) {
    
    double x2Sum = 0.0, xSum = 0.0, wSum = 0.0;
    for (unsigned int i = 0; i < inputVector.size(); ++i) {
        x2Sum += inputVector[i] * inputVector[i];
        xSum  += inputVector[i];
        wSum  += 1;
    }

    if (wSum > 0) {
        mean      = xSum / wSum;
    } else {
        mean = std::numeric_limits<double>::quiet_NaN();
    } 
    if (wSum > 1) {
        variance  = x2Sum / wSum - mean*mean;
        variance *= wSum / (wSum - 1);
    } else {
        variance = std::numeric_limits<double>::quiet_NaN();
    }
}

/**
 * Add a function to an image
 */
template <typename PixelT, typename FunctionT>
void lsst::imageproc::addFunction(
    lsst::fw::Image<PixelT> &image, ///< image
    lsst::fw::function::Function2<FunctionT> const &function ///< 2-d function
) {
    typedef typename lsst::fw::Image<PixelT>::pixel_accessor imageAccessorType;
    unsigned int numCols = image.getCols();
    unsigned int numRows = image.getRows();
    imageAccessorType imRow = image.origin();
    for (unsigned int row = 0; row < numRows; ++row, imRow.next_row()) {
        imageAccessorType imCol = imRow;
        double rowPos = lsst::fw::image::positionToIndex(row);
        for (unsigned int col = 0; col < numCols; ++col, imCol.next_col()) {
            double colPos = lsst::fw::image::positionToIndex(col);
            *imCol += static_cast<PixelT>(function(colPos, rowPos));
        }
    }
}
