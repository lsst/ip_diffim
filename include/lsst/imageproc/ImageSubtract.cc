// -*- lsst-c++ -*-
/**
 * \file
 *
 * Implementation of image subtraction
 *
 * Implementation of image subtraction
 *
 * \author Andrew Becker
 *
 * \ingroup imageproc
 */

#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/KernelFunctions.h>
#include <lsst/fw/PixelAccessors.h>
#include <lsst/fw/MinimizerFunctionBase.h>
#include <lsst/mwi/utils/Trace.h>
#include <lsst/mwi/exceptions/Exception.h>
#include <lsst/detection/Footprint.h>
#include <lsst/imageproc/PCA.h>
#include <vw/Math/Matrix.h> 
#include <vw/Math/Vector.h> 
#include <vw/Math/LinearAlgebra.h> 
#include <vw/Math/Functions.h> 
#include <boost/shared_ptr.hpp>
#include <Minuit/MnMigrad.h>
#include <Minuit/MnMinos.h>
#include <Minuit/MnPrint.h>
#include <Minuit/FunctionMinimum.h>
#include <string> // for upar.add
using namespace std;

/**
 * Computes spatially varying PSF matching kernel for image subtraction
 *
 * Note: Longer description here
 *
 * \return This describes the return value if there is one
 * \throw Any exceptions thrown must be described here
 * \throw Here too
 * \ingroup imageproc
 */
template <typename ImageT, typename MaskT, typename KernelT, typename FuncT>
void lsst::imageproc::computePSFMatchingKernelForMaskedImage(
    lsst::fw::MaskedImage<ImageT,MaskT> const &imageToConvolve, ///< Template image; is convolved
    lsst::fw::MaskedImage<ImageT,MaskT> const &imageToNotConvolve, ///< Science image; is not convolved
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelInBasisList, ///< Input set of basis kernels
    boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > &kernelPtr, ///< The output convolution kernel
    boost::shared_ptr<lsst::fw::function::Function2<FuncT> > &kernelFunctionPtr, ///< Function for spatial variation of kernel
    boost::shared_ptr<lsst::fw::function::Function2<FuncT> > &backgroundFunctionPtr ///< Function for spatial variation of background
    ) {

    // FROM THE POLICY
    unsigned int FOOTPRINT_DIFFIM_NPIX_MIN = 10;
    unsigned int FOOTPRINT_DIFFIM_GROW = 20;  // Eventually this will come from the kernel size, PSF FWHM, etc.

    // Find detections
    lsst::detection::DetectionSet<ImageT,MaskT> 
        detectionSet(imageToConvolve, lsst::detection::Threshold(15000, lsst::detection::Threshold::VALUE));

    // Reusable view of each Footprint
    typename lsst::fw::MaskedImage<ImageT,MaskT>::MaskedImagePtrT imageToConvolveFootprintPtr;
    typename lsst::fw::MaskedImage<ImageT,MaskT>::MaskedImagePtrT imageToNotConvolveFootprintPtr;
    
    // Do some culling of the Footprints; should be good in template image and science image
    // All it does now is looked for any masked pixels in the footprint
    vector<lsst::detection::Footprint::PtrType> footprintListOut;
    vector<lsst::detection::Footprint::PtrType> footprintListIn = detectionSet.getFootprints();

    for (vector<lsst::detection::Footprint::PtrType>::iterator i = footprintListIn.begin(); i != footprintListIn.end(); ) {
        if ((*i)->getNpix() < FOOTPRINT_DIFFIM_NPIX_MIN) {
            i = footprintListIn.erase(i); // SLOW; DON"T DO THIS IF YOU DON"T NEED TO
        }
        else {
            BBox2i footprintBBox = (*i)->getBBox();

            // NOTE - this will eventually be overridden by a grow method on Footprint
            footprintBBox.grow(footprintBBox.max() + vw::Vector2i(FOOTPRINT_DIFFIM_GROW,FOOTPRINT_DIFFIM_GROW));
            footprintBBox.grow(footprintBBox.min() - vw::Vector2i(FOOTPRINT_DIFFIM_GROW,FOOTPRINT_DIFFIM_GROW));
            
            try {
                imageToConvolveFootprintPtr = imageToConvolve.getSubImage(footprintBBox);
                imageToNotConvolveFootprintPtr = imageToNotConvolve.getSubImage(footprintBBox);
            } catch (lsst::mwi::exceptions::Exception &e) {
                // Likely that the grown subImage was off the image
                i = footprintListIn.erase(i);
            }

            if ( (lsst::imageproc::checkMaskedImageForDiffim(*imageToConvolveFootprintPtr) == false) || 
                 (lsst::imageproc::checkMaskedImageForDiffim(*imageToNotConvolveFootprintPtr) == false) ) {
                i = footprintListIn.erase(i);
            }
            else {

                // Create a new footprint with grow'd box
                lsst::detection::Footprint::PtrType fpGrow(
                    new lsst::detection::Footprint(footprintBBox)
                    );
                footprintListOut.push_back(fpGrow);
                
                ++i;
            }
        }
    }

    lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForMaskedImage", 4, 
                            (boost::format("Detected %d footprints above threshold") 
                             % footprintListOut.size()));

    lsst::imageproc::computePSFMatchingKernelForMaskedImage(imageToConvolve, imageToNotConvolve, 
                                                            kernelInBasisList, footprintListOut,
                                                            kernelPtr, kernelFunctionPtr, backgroundFunctionPtr);
}

template <typename ImageT, typename MaskT, typename KernelT, typename FuncT>
void lsst::imageproc::computePSFMatchingKernelForMaskedImage(
    lsst::fw::MaskedImage<ImageT,MaskT> const &imageToConvolve, ///< Template image; is convolved
    lsst::fw::MaskedImage<ImageT,MaskT> const &imageToNotConvolve, ///< Science image; is not convolved
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelInBasisList, ///< Input set of basis kernels
    vector<lsst::detection::Footprint::PtrType> const &footprintList, ///< List of input footprints to do diffim on
    boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > &kernelPtr, ///< The output convolution kernel
    boost::shared_ptr<lsst::fw::function::Function2<FuncT> > &kernelFunctionPtr, ///< Function for spatial variation of kernel
    boost::shared_ptr<lsst::fw::function::Function2<FuncT> > &backgroundFunctionPtr ///< Function for spatial variation of background
    ) {

    lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForMaskedImage", 2, 
                            "Entering subroutine computePSFMatchingKernelForMaskedImage");

    // FROM THE POLICY
    const double MAXIMUM_FOOTPRINT_RESIDUAL_MEAN = 1;     // Mean should be 0
    const double MAXIMUM_FOOTPRINT_RESIDUAL_VARIANCE = 2; // Variance should be 1
    const KernelT CONVOLVE_THRESHOLD = 0;
    const int EDGE_MASK_BIT = 1;
    bool BACKGROUND_FIT_CALCULATE_MINOS = false;

    double imSum;

    // Reusable view around each footprint
    typename lsst::fw::MaskedImage<ImageT,MaskT>::MaskedImagePtrT imageToConvolveStampPtr;
    typename lsst::fw::MaskedImage<ImageT,MaskT>::MaskedImagePtrT imageToNotConvolveStampPtr;

    // Information for each Footprint
    vector<lsst::imageproc::DiffImContainer<KernelT> > diffImContainerList;

    // Iterate over footprint

    int nFootprint = 0;
    for (vector<lsst::detection::Footprint::PtrType>::const_iterator i = footprintList.begin(); i != footprintList.end(); ++i) {
        BBox2i footprintBBox = (*i)->getBBox();
        imageToConvolveStampPtr    = imageToConvolve.getSubImage(footprintBBox);
        imageToNotConvolveStampPtr = imageToNotConvolve.getSubImage(footprintBBox);

        // DEBUGGING DEBUGGING DEBUGGING
        imageToConvolveStampPtr->writeFits( (boost::format("csFits_%d") % nFootprint).str() );
        imageToNotConvolveStampPtr->writeFits( (boost::format("ncsFits_%d") % nFootprint).str() );
        // DEBUGGING DEBUGGING DEBUGGING

        // Find best single kernel for this stamp
        double background;
        vector<double> kernelCoeffs(kernelInBasisList.size());
        lsst::imageproc::computePSFMatchingKernelForPostageStamp
            (*imageToConvolveStampPtr, *imageToNotConvolveStampPtr, kernelInBasisList, kernelCoeffs, background);

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

        vw::math::Vector<double,2> center = footprintBBox.center();
        // NOTE - verify center() returns col, row
        // Russ will not change function to span -1 to 1, so take this out for now
        //diffImFootprintContainer.colcNorm = 2 * center[0] / imageToConvolve.getCols() - 1;
        //diffImFootprintContainer.rowcNorm = 2 * center[1] / imageToConvolve.getRows() - 1;
        diffImFootprintContainer.colcNorm = center[0];
        diffImFootprintContainer.rowcNorm = center[1];
        lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForMaskedImage", 4, 
                                (boost::format("Developing kernel %d at position %f %f.  Background = %f") 
                                 % diffImFootprintContainer.id % center[0] % center[1] % background));

        // QA - calculate the residual of the subtracted image here
        lsst::fw::MaskedImage<ImageT, MaskT>
            convolvedImageStamp = lsst::fw::kernel::convolve(*imageToConvolveStampPtr, *footprintKernelPtr, CONVOLVE_THRESHOLD, EDGE_MASK_BIT);
        convolvedImageStamp -= (*imageToNotConvolveStampPtr);
        convolvedImageStamp *= -1;
        convolvedImageStamp -= background;

        double meanOfResiduals = 0;
        double varianceOfResiduals = 0;
        int nGoodPixels = 0;
        calculateMaskedImageResiduals(convolvedImageStamp, nGoodPixels, meanOfResiduals, varianceOfResiduals);
        diffImFootprintContainer.footprintResidualMean = meanOfResiduals;
        diffImFootprintContainer.footprintResidualVariance = varianceOfResiduals;

        if (fabs(meanOfResiduals) > MAXIMUM_FOOTPRINT_RESIDUAL_MEAN) {
            lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForMaskedImage", 4, 
                                    (boost::format("Kernel %d, bad mean residual of footprint : %f") 
                                     % diffImFootprintContainer.id % meanOfResiduals));
            diffImFootprintContainer.isGood = false;
        }
        if (varianceOfResiduals > MAXIMUM_FOOTPRINT_RESIDUAL_VARIANCE) {
            lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForMaskedImage", 4, 
                                    (boost::format("Kernel %d, bad residual variance of footprint : %f") 
                                     % diffImFootprintContainer.id % varianceOfResiduals));
            diffImFootprintContainer.isGood = false;
        }
        diffImContainerList.push_back(diffImFootprintContainer);

        // DEBUGGING DEBUGGING DEBUGGING
        lsst::fw::Image<KernelT> kImage = footprintKernelPtr->computeNewImage(imSum);
        kImage.writeFits( (boost::format("kFits_%d.fits") % nFootprint).str() );
        convolvedImageStamp.writeFits( (boost::format("d1Fits_%d") % nFootprint).str() );
        // DEBUGGING DEBUGGING DEBUGGING

        nFootprint += 1;
    }

    // In the end we want to test if the kernelInBasisList is Delta Functions; if so, do PCA
    // For DC2 we just do it
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > kernelOutBasisList;
    lsst::imageproc::computePcaKernelBasis(diffImContainerList, kernelOutBasisList);

    // Compute spatial variation of the kernel if requested
    if (kernelFunctionPtr != NULL) {
        computeSpatiallyVaryingPSFMatchingKernel(
            diffImContainerList,
            kernelOutBasisList, 
            kernelPtr,
            kernelFunctionPtr);
    }

    // Compute spatial variation of the background if requested
    if (backgroundFunctionPtr != NULL) {
        int nGood = 0;
        vector<double> backgrounds;
        vector<double> variances;
        vector<double> position1;
        vector<double> position2;
        double bgSum = 0;
        for (typename vector<lsst::imageproc::DiffImContainer<KernelT> >::iterator i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
            if ((*i).isGood == true) {
                nGood += 1;
                backgrounds.push_back((*i).background);
                bgSum += (*i).background;
                lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForMaskedImage", 4, 
                                        (boost::format(" Background %d = %f") 
                                         % (*i).id % (*i).background));
                
                variances.push_back((*i).footprintResidualVariance);   // approximation
                position1.push_back((*i).diffImFootprintPtr->getBBox().center()[0]);
                position2.push_back((*i).diffImFootprintPtr->getBBox().center()[1]);
            }
        }

        if (nGood == 0) {
            throw lsst::mwi::exceptions::Exception("No good footprints for background estimation");
        }

        double nSigmaSquared = 1;
        lsst::fw::function::MinimizerFunctionBase2<FuncT> 
            bgFcn(backgrounds, variances, position1, position2, nSigmaSquared, backgroundFunctionPtr);
        MnUserParameters bgPar;
        // Start 0th parameter at mean background with 10% variation
        bgPar.add("p0", bgSum / nGood, 0.1 * bgSum / nGood);
        for (unsigned int npar = 1; npar < backgroundFunctionPtr->getNParameters(); ++npar) {
            // Start other parameters at value 0 with small variation
            bgPar.add((boost::format("p%d") % npar).str().c_str(), 0, 0.1);
        }
        MnMigrad migrad(bgFcn, bgPar);
        FunctionMinimum bgFit = migrad();
        MnMinos minos(bgFcn, bgFit); 

        if (!(bgFit.isValid())) {
            lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForMaskedImage", 4, 
                                    (boost::format("WARNING : Spatial fit to background failed to converge")));
        }
            

        vector<double> bgFitParameters;
        for (unsigned int npar = 0; npar < backgroundFunctionPtr->getNParameters(); ++npar) {
            bgFitParameters.push_back(bgFit.userState().value((boost::format("p%d") % npar).str().c_str()));
            if ( (BACKGROUND_FIT_CALCULATE_MINOS == true) && (bgFit.isValid()) ) {
                std::pair<double,double> e = minos(npar);
                lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForMaskedImage", 4, 
                                        (boost::format("Fit Background Parameter %d : %f (%f,%f)\n\n") 
                                         % npar % bgFit.userState().value((boost::format("p%d") % npar).str().c_str())
                                         % e.first % e.second));
            }
            else {
                lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForMaskedImage", 4, 
                                        (boost::format("Fit Background Parameter %d : %f (-%f,%f)\n\n") 
                                         % npar % bgFit.userState().value((boost::format("p%d") % npar).str().c_str())
                                         % bgFit.userState().error((boost::format("p%d") % npar).str().c_str())
                                         % bgFit.userState().error((boost::format("p%d") % npar).str().c_str())));
            }
        }
        backgroundFunctionPtr->setParameters(bgFitParameters);
    }
    
}

/**
 * Single line description with no period
 *
 * Note: Long description
 *
 * \return This describes the return value if there is one
 * \throw Any exceptions thrown must be described here
 * \throw Here too
 * \ingroup I guess imageproc for me
 */
template <typename ImageT, typename MaskT, typename KernelT>
void lsst::imageproc::computePSFMatchingKernelForPostageStamp(
    lsst::fw::MaskedImage<ImageT, MaskT> const &imageToConvolve, ///< Goes with the code
    lsst::fw::MaskedImage<ImageT, MaskT> const &imageToNotConvolve, ///< This is for doxygen
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelInBasisList, ///< Input kernel basis set
    vector<double> &kernelCoeffs, ///< Output kernel basis coefficients
    double &background ///< Difference in the backgrounds
    ) { 

    int nKernelParameters=0, nBackgroundParameters=0, nParameters=0;
    const KernelT threshold = 0.0;
    const int edgeMaskBit = 1;

    lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 2, 
                            "Entering subroutine computePSFMatchingKernelForPostageStamp");
    
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
    vector<boost::shared_ptr<lsst::fw::MaskedImage<ImageT, MaskT> > > convolvedImageList(nKernelParameters);
    // and an iterator over this
    typename vector<boost::shared_ptr<lsst::fw::MaskedImage<ImageT, MaskT> > >::iterator citer = convolvedImageList.begin();
    
    // Iterator for input kernel basis
    typename vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > >::const_iterator kiter = kernelInBasisList.begin();
    int kId = 0; // debugging
    // Create C_ij in the formalism of Alard & Lupton
    for (; kiter != kernelInBasisList.end(); ++kiter, ++citer, ++kId) {

        lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 3, 
                                "Convolving an Object with Basis");
        
        // NOTE : we could also *precompute* the entire template image convolved with these functions
        //        and save them somewhere to avoid this step each time.  however, our paradigm is to
        //        compute whatever is needed on the fly.  hence this step here.
        boost::shared_ptr<lsst::fw::MaskedImage<ImageT, MaskT> > imagePtr(
            new lsst::fw::MaskedImage<ImageT, MaskT>
            (lsst::fw::kernel::convolve(imageToConvolve, **kiter, threshold, edgeMaskBit))
            );

        lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 3, 
                                "Convolved an Object with Basis");

        *citer = imagePtr;
        
        imagePtr->writeFits( (boost::format("cFits_%d") % kId).str() );
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
    vector<lsst::fw::MaskedPixelAccessor<ImageT, MaskT> > convolvedAccessorRowList;
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

    // integral over image's dx and dy
    for (unsigned int row = startRow; row < endRow; ++row) {

        // An accessor for each convolution plane
        vector<lsst::fw::MaskedPixelAccessor<ImageT, MaskT> > convolvedAccessorColList = convolvedAccessorRowList;

        // An accessor for each input image; places the col accessor at the correct row
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imageToConvolveCol = imageToConvolveRow;
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imageToNotConvolveCol = imageToNotConvolveRow;

        for (unsigned int col = startCol; col < endCol; ++col) {

            ImageT ncCamera = *imageToNotConvolveCol.image;
            ImageT ncVariance = *imageToNotConvolveCol.variance;
            MaskT  ncMask = *imageToNotConvolveCol.mask;
            ImageT cVariance = *imageToConvolveCol.variance;
            ImageT iVariance = 1.0 / (cVariance + ncVariance);

            // Do we skip these?
            if (ncMask & 0x1) {
                continue;
            }

            lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 5, 
                                    boost::format("Accessing image row %d col %d : %f %f") 
                                    % row % col % ncCamera % ncVariance);

            // kernel index i
            typename vector<lsst::fw::MaskedPixelAccessor<ImageT, MaskT> >::iterator
                kiteri = convolvedAccessorColList.begin();

            for (int kidxi = 0; kiteri != convolvedAccessorColList.end(); ++kiteri, ++kidxi) {
                ImageT cdCamerai = *kiteri->image;
                ImageT cdVariancei = *kiteri->variance;
                MaskT cdMaski = *kiteri->mask;
                // Do we skip these?
                if (cdMaski & 0x1) {
                    continue;
                }
                lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 6, 
                                        boost::format("Accessing convolved image %d : %f %f") 
                                        % kidxi % cdCamerai % cdVariancei);

                // Not sure what is the best variance to use
                //B[kidxi] += ncCamera * cdCamerai / (ncVariance + cdVariancei);
                B[kidxi] += ncCamera * cdCamerai * iVariance;
                
                // kernel index j 
                typename vector<lsst::fw::MaskedPixelAccessor<ImageT, MaskT> >::iterator kiterj = kiteri;
                for (int kidxj = kidxi; kiterj != convolvedAccessorColList.end(); ++kiterj, ++kidxj) {
                    ImageT cdCameraj = *kiterj->image;
                    ImageT cdVariancej = *kiterj->variance;
                    MaskT cdMaskj = *kiterj->mask;
                    // Do we skip these?
                    if (cdMaskj & 0x1) {
                        continue;
                    }
                    lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 6, 
                                            boost::format("Accessing convolved image %d : %f %f") 
                                            % kidxj % cdCameraj % cdVariancej);

                    //M[kidxi][kidxj] += cdCamerai * cdCameraj / (cdVariancei + cdVariancej);
                    M[kidxi][kidxj] += cdCamerai * cdCameraj * iVariance;
                } 
                // Constant background term; effectively kidxj+1
                // NOTE : is this cVariance or ncVariance????
                //M[kidxi][nParameters-1] += cdCamerai / (cVariance + cdVariancei);
                M[kidxi][nParameters-1] += cdCamerai * iVariance;
            } 

            // Background term; effectively kidxi+1
            //B[nParameters-1] += ncCamera / (ncVariance + cVariance); 
            //M[nParameters-1][nParameters-1] += 1.0 / (cVariance + cVariance);
            B[nParameters-1] += ncCamera * iVariance;
            M[nParameters-1][nParameters-1] += 1.0 * iVariance;
            lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 5, 
                                    boost::format("Background terms : %f %f") 
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

    // Fill in rest of M
    for (int kidxi=0; kidxi < nParameters; ++kidxi) 
        for (int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) 
            M[kidxj][kidxi] = M[kidxi][kidxj];

    cout << "B : " << B << endl;
    cout << "M : " << M << endl;

    // Invert M
    vw::math::Matrix<double> Minv;
    try {
        Minv = vw::math::inverse(M);
    } catch (lsst::mwi::exceptions::Exception &e) {
        lsst::mwi::utils::Trace("lsst.imageproc.computePSFMatchingKernelForPostageStamp", 1, 
                                ("Invert matrix failed"));
        // NOTE - lets be smarter here.
        return;
    }

    // Solve for x in Mx = B
    vw::math::Vector<double> Soln = Minv * B;

    // Worry about translating here...
    for (int ki = 0; ki < nKernelParameters; ++ki) {
        kernelCoeffs[ki] = Soln[ki];
    }
    background = Soln[nParameters-1];
}

void lsst::imageproc::getCollectionOfMaskedImagesForPSFMatching(
    vector<lsst::detection::Footprint::PtrType> &footprintList ///< Collection of footprints to build kernels around
    ) {

    double radius = 20;
    
    
    const vw::BBox2i region1(static_cast<int>(78.654-radius/2),
                             static_cast<int>(3573.945-radius/2),
                             static_cast<int>(radius),
                             static_cast<int>(radius));
    lsst::detection::Footprint::PtrType fp1(
        new lsst::detection::Footprint(region1)
        );
    footprintList.push_back(fp1);
    
    
    const vw::BBox2i region4(static_cast<int>(367.756-radius/2),
                             static_cast<int>(3827.671-radius/2),
                             static_cast<int>(radius),
                             static_cast<int>(radius));
    lsst::detection::Footprint::PtrType fp4(
        new lsst::detection::Footprint(region4)
        );
    footprintList.push_back(fp4);
    
    const vw::BBox2i region5(static_cast<int>(381.062-radius/2),
                             static_cast<int>(3212.948-radius/2),
                             static_cast<int>(radius),
                             static_cast<int>(radius));
    lsst::detection::Footprint::PtrType fp5(
        new lsst::detection::Footprint(region5)
        );
    footprintList.push_back(fp5);
    
    const vw::BBox2i region6(static_cast<int>(404.433-radius/2),
                             static_cast<int>(573.462-radius/2),
                             static_cast<int>(radius),
                             static_cast<int>(radius));
    lsst::detection::Footprint::PtrType fp6(
        new lsst::detection::Footprint(region6)
        );
    footprintList.push_back(fp6);
    
    const vw::BBox2i region7(static_cast<int>(420.967-radius/2),
                             static_cast<int>(3306.310-radius/2),
                             static_cast<int>(radius),
                             static_cast<int>(radius));
    lsst::detection::Footprint::PtrType fp7(
        new lsst::detection::Footprint(region7)
        );
    footprintList.push_back(fp7);
    
    const vw::BBox2i region9(static_cast<int>(546.657-radius/2),
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
}

template <typename KernelT>
void lsst::imageproc::computePcaKernelBasis(
    vector<lsst::imageproc::DiffImContainer<KernelT> > &diffImContainerList, ///< List of input footprints
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelPcaBasisList ///< Output principal components as kernel images
    ) {
    
    // FROM THE POLICY
    const unsigned int MINIMUM_NUMBER_OF_BASES = 3; 
    const float MAXIMUM_FRACTION_OF_EIGENVALUES = 0.99;
    const float MINIMUM_ACCEPTIBLE_EIGENVALUE = 1e-12;
    const unsigned int MAXIMUM_ITER_PCA = 5;
    const double MAXIMUM_KERNEL_RESIDUAL = 1;
    const double MAXIMUM_KERNEL_RESIDUAL_VARIANCE = 2;
    
    // Image accessor
    typedef typename vw::ImageView<KernelT>::pixel_accessor imageAccessorType;
    typedef typename vector<lsst::imageproc::DiffImContainer<KernelT> >::iterator iDiffImContainer;
    // Iterator over struct
    
    // Matrix to invert.  Number of rows = number of pixels; number of columns = number of kernels
    // All calculations here are in double
    vw::math::Matrix<double> M;
    vw::math::Matrix<double> eVec;
    vw::math::Vector<double> eVal;
    vw::math::Vector<double> mMean;
    vw::math::Matrix<double> kernelCoefficientMatrix;
    
    // Initialize for while loop
    unsigned int nIter = 0;
    unsigned int nReject = 1;
    unsigned int nGood = diffImContainerList.size();
    unsigned int nKCols = 0, nKRows = 0;
    unsigned int nPixels;
    unsigned int nCoeffToUse;
    lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 2, 
                            "Entering subroutine computePcaKernelBasis");
    
    // Iterate over PCA inputs until all are good
    while ( (nIter < MAXIMUM_ITER_PCA) and (nReject != 0) ) {
        
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
            throw lsst::mwi::exceptions::Exception("No good kernels for PCA");
        }
                    
                    
        nPixels = nKCols * nKRows;
        M.set_size(nPixels, nGood);
        
        // fill up matrix for PCA
        double imSum;
        int ki = 0;
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++ki, ++i) {
            if ((*i).isGood == false) {
                continue;
            }
            
            lsst::fw::Image<KernelT> kImage = (*i).diffImKernelPtr->computeNewImage(imSum);
            
            //assert(nKRows == kImage.getRows());
            //assert(nKCols == kImage.getCols());
            
            int mIdx = 0;
            imageAccessorType imageAccessorCol(kImage.origin());
            for (unsigned int col = 0; col < nKCols; ++col) {
                
                imageAccessorType imageAccessorRow(imageAccessorCol);
                for (unsigned int row = 0; row < nKRows; ++row, ++mIdx) {
                    
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
                    imageAccessorRow.next_row();
                }
                imageAccessorCol.next_col();
            }
        }
    
        eVec.set_size(nPixels, nGood);
        eVal.set_size(nGood);
        mMean.set_size(nPixels);
        
        // M is mean-subtracted if subtractMean = true
        // Eigencomponents in columns of eVec
        lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 4, "Computing pricipal components");
        lsst::imageproc::computePca(M, mMean, eVal, eVec, true);
        lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 4, "Computed pricipal components");

        cout << "Eigenvalues : " << eVal << endl;

        nCoeffToUse = MINIMUM_NUMBER_OF_BASES;
        // Potentially override with larger number if the spectrum of eigenvalues requires it
        double evalSum = 0;
        for (unsigned int i = 0; i < eVal.size(); ++i) {
            evalSum += eVal[i];
        }
        double evalFrac = 0;
        for (unsigned int i = 0; i < nCoeffToUse; ++i) {
            evalFrac += eVal[i] / evalSum;
            if (eVal[i] < MINIMUM_ACCEPTIBLE_EIGENVALUE) {
                lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 4, 
                                        (boost::format("WARNING : Using eigenvector whose eigenvalue (%.3e) is smaller than acceptible (%.3e)")
                                         % eVal[i] % MINIMUM_ACCEPTIBLE_EIGENVALUE));
            }
                
        }
        if (evalFrac < MAXIMUM_FRACTION_OF_EIGENVALUES) {
            for (; nCoeffToUse < eVal.size(); ++nCoeffToUse) {
                evalFrac += eVal[nCoeffToUse] / evalSum;
                if (eVal[nCoeffToUse] < MINIMUM_ACCEPTIBLE_EIGENVALUE) {
                    lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 4, 
                                            (boost::format("WARNING : Using eigenvector whose eigenvalue (%.3e) is smaller than acceptible (%.3e)")
                                             % eVal[nCoeffToUse] % MINIMUM_ACCEPTIBLE_EIGENVALUE));
                }
                if (evalFrac > MAXIMUM_FRACTION_OF_EIGENVALUES) {
                    break;
                }
            }
        }
        lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 4, 
                                (boost::format("Using %d basis functions (plus mean) representing %.5f of the variance") 
                                 % nCoeffToUse % evalFrac));

        // We now have the basis functions determined
        // Next determine the coefficients that go in front of all of the individual kernels
        kernelCoefficientMatrix.set_size(nGood, nCoeffToUse);
        lsst::imageproc::decomposeMatrixUsingBasis(M, eVec, nCoeffToUse, kernelCoefficientMatrix);
        
        // We next do quality control here; we reconstruct the input kernels with the truncated basis function set
        // Remember that M is mean-subtracted already
        vw::math::Matrix<double> approxM(nPixels, nGood); 
        lsst::imageproc::approximateMatrixUsingBasis(eVec, kernelCoefficientMatrix, nCoeffToUse, approxM);

        cout << "eigenKernel coefficients" << endl;
        cout << kernelCoefficientMatrix << endl;

        nReject = 0;
        unsigned int iKernel = 0;
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++iKernel, ++i) {
            double x2Sum = 0, xSum = 0, wSum = 0;
            double residual;
            for (unsigned int jKernel = 0; jKernel < M.rows(); ++jKernel) {
                residual = M(iKernel,jKernel) - approxM(iKernel,jKernel);
                x2Sum += residual * residual;
                xSum  += residual;
                wSum  += 1;
            }

            (*i).kernelResidual = xSum / wSum;
            (*i).kernelResidualVariance = x2Sum / wSum - (*i).kernelResidual * (*i).kernelResidual;
            if ((*i).kernelResidual > MAXIMUM_KERNEL_RESIDUAL) {
                lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 4, 
                                        (boost::format("Kernel %d, bad mean residual of PCA model : %f") 
                                         % (*i).id % (*i).kernelResidual));
                (*i).isGood = false;
                nReject += 1;
            }
            if ((*i).kernelResidualVariance > MAXIMUM_KERNEL_RESIDUAL_VARIANCE) {
                lsst::mwi::utils::Trace("lsst.imageproc.computePcaKernelBasis", 4, 
                                        (boost::format("Kernel %d, bad residual variance of PCA model : %f") 
                                         % (*i).id % (*i).kernelResidualVariance));
                (*i).isGood = false;
                nReject += 1;
            }
            cout << " Static Kernel Residual " << (*i).id << " = " << (*i).kernelResidual << " " << (*i).kernelResidualVariance << endl;
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

    // DEBUGGING DEBUGGING DEBUGGING 
    meanImage.writeFits( (boost::format("mFits.fits")).str() );
    // DEBUGGING DEBUGGING DEBUGGING 

    // Turn each eVec into an Image and then into a Kernel
    //for (unsigned int ki = 0; ki < eVec.cols(); ++ki) {
    for (unsigned int ki = 0; ki < nCoeffToUse; ++ki) { // Dont use all eVec.cols(); only the number of bases that you want
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

        // DEBUGGING DEBUGGING DEBUGGING 
        basisImage.writeFits( (boost::format("eFits_%d.fits") % ki).str() );
        // DEBUGGING DEBUGGING DEBUGGING 
    }

    // Finally, create a new LinearCombinationKernel for each Footprint
    unsigned int iKernel = 0;
    for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++iKernel, ++i) {
        if ((*i).isGood == true) {
            vector<double> kernelCoefficients;

            kernelCoefficients.push_back(1); // Mean image
            for (unsigned int jKernel = 0; jKernel < nCoeffToUse; ++jKernel) {
                kernelCoefficients.push_back(kernelCoefficientMatrix(iKernel,jKernel));
            }        
            
            // Create a linear combination kernel from this 
            boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > PcaKernelPtr(
                new lsst::fw::LinearCombinationKernel<KernelT>(kernelPcaBasisList, kernelCoefficients)
                );
            
            (*i).diffImPcaKernelPtr = PcaKernelPtr;
        }
    }
}

template <typename KernelT, typename FuncT>
void lsst::imageproc::computeSpatiallyVaryingPSFMatchingKernel(
    vector<lsst::imageproc::DiffImContainer<KernelT> > &diffImContainerList, ///< Information on each kernel
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > const &kernelOutBasisList, ///< Input basis kernel set
    boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> > &spatiallyVaryingKernelPtr, ///< Output kernel
    boost::shared_ptr<lsst::fw::function::Function2<FuncT> > &kernelFunctionPtr ///< Function for spatial variation of kernel
    )
 {
     // FROM THE POLICY
     const unsigned int MAXIMUM_ITER_SPATIAL_FIT = 5;
     const double MAXIMUM_SPATIAL_KERNEL_RESIDUAL = 1;
     const double MAXIMUM_SPATIAL_KERNEL_RESIDUAL_VARIANCE = 1;
     bool SPATIAL_KERNEL_FIT_CALCULATE_MINOS = false;

     // Container iterator
     typedef typename vector<lsst::imageproc::DiffImContainer<KernelT> >::iterator iDiffImContainer;
     // Kernel iterator
     typedef typename vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > >::const_iterator iKernelPtr;

     // Initialize for while loop
     unsigned int nIter = 0;
     unsigned int nReject = 1;
     unsigned int nGood = diffImContainerList.size();
     unsigned int nKernel = kernelOutBasisList.size();
     unsigned int nParameters = kernelFunctionPtr->getNParameters();
     
     lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPSFMatchingKernel", 2, 
                             "Entering subroutine computeSpatiallyVaryingPSFMatchingKernel");
     
     vector<double> measurements;
     vector<double> variances;
     vector<double> position1;
     vector<double> position2;

         
     // Set up the spatially varying kernel
     spatiallyVaryingKernelPtr = boost::shared_ptr<lsst::fw::LinearCombinationKernel<KernelT> >
         (new lsst::fw::LinearCombinationKernel<KernelT>(kernelOutBasisList, kernelFunctionPtr));

     // Iterate over PCA inputs until all are good
     while ( (nIter < MAXIMUM_ITER_SPATIAL_FIT) and (nReject != 0) ) {

         // NOTE - For a delta function basis set, the mean image is the first entry in kernelPcaBasisList.  
         // This should *not* vary spatially.  Should we code this up or let the fitter determine that it does not vary..?

         position1.clear();
         position2.clear();
         nGood = 0;
         for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
             if ((*i).isGood == true) {
                 position1.push_back((*i).colcNorm);
                 position2.push_back((*i).rowcNorm);
                 nGood += 1;
             }
         }
         
         //vector<double> nParList(nParameters);
         vector<vector<double> > fitParameters = spatiallyVaryingKernelPtr->getSpatialParameters();
         // DEBUGGING
         assert(fitParameters.size() == nKernel);
         assert(fitParameters[0].size() == nParameters);

         unsigned int nk = 0;
         for (iKernelPtr ik = kernelOutBasisList.begin(); ik != kernelOutBasisList.end(); ++nk, ++ik) {
             
             cout << "    Spatial fit to coefficient " << nk << ", Values are : ";
             
             measurements.clear();
             variances.clear();
             for (iDiffImContainer id = diffImContainerList.begin(); id != diffImContainerList.end(); ++id) {
                 if ((*id).isGood == true) {
                     vector<double> kernelParameters = (*id).diffImPcaKernelPtr->getKernelParameters();  // inefficient
                     measurements.push_back(kernelParameters[nk]);
                     variances.push_back((*id).footprintResidualVariance); // approximation
                     cout << kernelParameters[nk] << " ";
                 }
             }
             cout << endl;

             // NOTE - if we have fewer measurements than kernelFunctionPtr->getNParameters(), we should do something about it...

             double nSigmaSquared = 1;
             lsst::fw::function::MinimizerFunctionBase2<FuncT> 
                 kernelFcn(measurements, variances, position1, position2, nSigmaSquared, kernelFunctionPtr);
             
             MnUserParameters kernelPar;
             // Start 0th spatial parameter at 1; 10% variation
             kernelPar.add("p0", 1, 0.1);
             for (unsigned int npar = 1; npar < kernelFunctionPtr->getNParameters(); ++npar) {
                 // Start other spatial parameters at value 0 
                 kernelPar.add((boost::format("p%d") % npar).str().c_str(), 0, 0.1);
             }

             //for (unsigned int i = 0; i < nKernel; ++i) {
             //cout << " FITTING " << i << " " << position1[i] << " " << position2[i] << " " << measurements[i] << " " << variances[i] << endl;
             //}
             //cout << " FIT PAR " << kernelPar << endl;

             MnMigrad migrad(kernelFcn, kernelPar);
             FunctionMinimum kernelFit = migrad();
             MnMinos minos(kernelFcn, kernelFit); 

             if (!(kernelFit.isValid())) {
                 lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPSFMatchingKernel", 4, 
                                         (boost::format("WARNING : Spatial fit to Kernel %d failed to converge")) % nk);
             }

             for (unsigned int i = 0; i < nParameters; ++i) {
                 fitParameters[nk][i] = kernelFit.userState().value((boost::format("p%d") % i).str().c_str());
                 if ( (SPATIAL_KERNEL_FIT_CALCULATE_MINOS == true) && (kernelFit.isValid()) ) {
                     std::pair<double,double> e = minos(i);
                     lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPSFMatchingKernel", 4, 
                                             (boost::format("Fit to Kernel %d, spatial Parameter %d : %f (%f,%f)\n\n") 
                                              % nk % i % kernelFit.userState().value((boost::format("p%d") % i).str().c_str())
                                              % e.first % e.second));
                 }
                 else {
                     lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPSFMatchingKernel", 4, 
                                             (boost::format("Fit to Kernel %d, spatial Parameter %d : %f (-%f,%f)\n\n") 
                                              % nk % i % kernelFit.userState().value((boost::format("p%d") % i).str().c_str())
                                              % kernelFit.userState().error((boost::format("p%d") % i).str().c_str())
                                              % kernelFit.userState().error((boost::format("p%d") % i).str().c_str())));
                 }
             }
         }

         nReject = 0;
         spatiallyVaryingKernelPtr->setSpatialParameters(fitParameters);

         double imSum;
         for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
             if ((*i).isGood == true) {
                 lsst::fw::Image<KernelT> diffImage = spatiallyVaryingKernelPtr->computeNewImage(imSum, 
                                                                                                 (*i).colcNorm,
                                                                                                 (*i).rowcNorm,
                                                                                                 false);

                 // DEBUGGING
                 cout << "CAW" << imSum << endl;

                 vector<double> kernelParams = spatiallyVaryingKernelPtr->getKernelParameters((*i).colcNorm, (*i).rowcNorm);
                 for (unsigned int caw = 0; caw < kernelParams.size(); ++caw) {
                     cout << boost::format("Kernel param %d at col=%d, row=%d; %7.2f\n")
                         % caw % (*i).colcNorm % (*i).rowcNorm % kernelParams[caw];
                 }
                 
                 lsst::fw::kernel::printKernel(
                     *spatiallyVaryingKernelPtr, static_cast<double>((*i).colcNorm), static_cast<double>((*i).rowcNorm), false, "%7.3e ");

                 vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > kList = spatiallyVaryingKernelPtr->getKernelList();
                 for (unsigned int caw = 0; caw < kList.size(); ++caw) {
                     cout << " KERNEL " << caw << endl;
                     lsst::fw::kernel::printKernel(
                         *kList[caw], static_cast<double>((*i).colcNorm), static_cast<double>((*i).rowcNorm), false, "%7.3e ");
                 }
                 // DEBUGGING

                 diffImage -= (*i).diffImPcaKernelPtr->computeNewImage(imSum);
                 diffImage.writeFits( (boost::format("ksmFits_%d.fits") % (*i).id).str() );

                 double meanOfResiduals = 0;
                 double varianceOfResiduals = 0;
                 int nGoodPixels = 0;
                 calculateImageResiduals(diffImage, nGoodPixels, meanOfResiduals, varianceOfResiduals);
                 (*i).spatialKernelResidual = meanOfResiduals;
                 (*i).spatialKernelResidualVariance = varianceOfResiduals;
                 if ((*i).spatialKernelResidual > MAXIMUM_SPATIAL_KERNEL_RESIDUAL) {
                     lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPSFMatchingKernel", 4, 
                                             (boost::format("Kernel %d, bad mean residual of spatial fit : %f") 
                                              % (*i).id % (*i).spatialKernelResidual));
                     (*i).isGood = false;
                     nReject += 1;
                 }
                 if ((*i).spatialKernelResidualVariance > MAXIMUM_SPATIAL_KERNEL_RESIDUAL_VARIANCE) {
                     lsst::mwi::utils::Trace("lsst.imageproc.computeSpatiallyVaryingPSFMatchingKernel", 4, 
                                             (boost::format("Kernel %d, bad residual variance of spatial fit : %f") 
                                              % (*i).id % (*i).spatialKernelResidualVariance));
                     (*i).isGood = false;
                     nReject += 1;
                 }
                 cout << " Spatial Kernel Residual " << (*i).id << " = " << (*i).spatialKernelResidual << " " << (*i).spatialKernelResidualVariance << endl;
             }
         }
         nIter += 1;
     }
}

template <typename KernelT>
void lsst::imageproc::generateDeltaFunctionKernelSet(
    unsigned int const nRows, ///< Number of rows in the kernel basis
    unsigned int const nCols, ///< Number of colunms in the kernel basis
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelBasisList ///< Output kernel basis function, length nRows * nCols
    )
{
    int colCtr = (nCols - 1) / 2;
    int rowCtr = (nRows - 1) / 2;
    for (unsigned int row = 0; row < nRows; ++row) {
        int y = row - rowCtr;
        
        for (unsigned int col = 0; col < nCols; ++col) {
            int x = col - colCtr;
            
            typename lsst::fw::Kernel<KernelT>::KernelFunctionPtrType kfuncPtr(
                new lsst::fw::function::IntegerDeltaFunction2<KernelT>(x, y)
                );
            
            boost::shared_ptr<lsst::fw::Kernel<KernelT> > kernelPtr(
                new lsst::fw::AnalyticKernel<KernelT>(kfuncPtr, nCols, nRows)
                );
            
            kernelBasisList.push_back(kernelPtr);
        }
    }
}

template <typename KernelT>
void lsst::imageproc::generateAlardLuptonKernelSet(
    int const nRows, ///< Number of rows in the kernel basis
    int const nCols, ///< Number of columns in the kernel basis
    vector<double> const sigGauss, ///< Width of gaussians in basis; size = number of Gaussians
    vector<double> const degGauss, ///< Degree of spatial variation within each Gaussian; size = sigGauss.size()
    vector<boost::shared_ptr<lsst::fw::Kernel<KernelT> > > &kernelBasisList ///< Output kernel basis function, length sum_n 0.5*(deg_n+1)*(deg_n+2)
    )
{
    ; // TO DO 
}

template <typename ImageT, typename MaskT>
bool lsst::imageproc::checkMaskedImageForDiffim(
    lsst::fw::MaskedImage<ImageT,MaskT> const &inputImage
    )
{
    // FROM THE POLICY
    //const int BAD_MASK_BIT = 0x1;
    
    lsst::fw::MaskedPixelAccessor<ImageT, MaskT> accessorCol(inputImage);
    for (unsigned int col = 0; col < inputImage.getCols(); ++col) {
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> accessorRow = accessorCol;
        for (unsigned int row = 0; row < inputImage.getRows(); ++row) {
            //if (((*accessorCol.mask) & BAD_MASK_BIT) != 0) {
            if ((*accessorRow.mask) > 0) {
                return false;
            }
            accessorRow.nextRow();
        }
        accessorCol.nextCol();
    }
    return true;
}

template <typename ImageT, typename MaskT>
void lsst::imageproc::calculateMaskedImageResiduals(
    lsst::fw::MaskedImage<ImageT,MaskT> const &inputImage, ///< Input difference image to be analyzed
    int &nGoodPixels, ///< Number of good pixels in the image
    double &meanOfResiduals, ///< Mean residual/variance; ideally 0
    double &varianceOfResiduals ///< Average variance of residual/variance; ideally 1
    ) {

    // FROM THE POLICY
    const int BAD_MASK_BIT = 0x1;
    // BASICALLY A HACK FOR DC2 TESTING; SENDING THE SAME IMAGE TO DIFFIM RESULTS IN 0 VARIANCE; MESSING UP FITTING
    // Add in quadrature as noise
    double MIN_VARIANCE = 1e-20;
    
    double x2Sum=0, xSum=0, wSum=0;

    nGoodPixels = 0;
    lsst::fw::MaskedPixelAccessor<ImageT, MaskT> accessorCol(inputImage);
    for (unsigned int col = 0; col < inputImage.getCols(); ++col) {
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> accessorRow = accessorCol;
        for (unsigned int row = 0; row < inputImage.getRows(); ++row) {
            if (((*accessorRow.mask) & BAD_MASK_BIT) == 0) {
                nGoodPixels += 1;
                x2Sum += (*accessorRow.image) * (*accessorRow.image) / (*accessorRow.variance);
                xSum  += (*accessorRow.image) / (*accessorRow.variance);
                wSum  += 1. / (*accessorRow.variance);
            }
            accessorRow.nextRow();
        }
        accessorCol.nextCol();
    }
    meanOfResiduals = -1;
    varianceOfResiduals = -1;

    if (nGoodPixels > 0) {
        meanOfResiduals      = xSum / wSum;
    }
    if (nGoodPixels > 1) {
        varianceOfResiduals  = x2Sum / wSum - meanOfResiduals*meanOfResiduals;
        varianceOfResiduals /= (nGoodPixels - 1);
        varianceOfResiduals += MIN_VARIANCE;
    }
}

template <typename ImageT>
void lsst::imageproc::calculateImageResiduals(
    lsst::fw::Image<ImageT> const &inputImage, ///< Input difference image to be analyzed
    int &nGoodPixels, ///< Number of good pixels in the image
    double &meanOfResiduals, ///< Mean residual
    double &varianceOfResiduals ///< Average variance of residual
    ) {

    double x2Sum=0, xSum=0, wSum=0;

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

    meanOfResiduals = -1;
    varianceOfResiduals = -1;
    
    //cout << "CAW RESID " << nGoodPixels << " " << x2Sum << " " << xSum << " " << wSum << endl;

    if (nGoodPixels > 0) {
        meanOfResiduals      = xSum / wSum;
    }
    if (nGoodPixels > 1) {
        varianceOfResiduals  = x2Sum / wSum - meanOfResiduals*meanOfResiduals;
        varianceOfResiduals /= (nGoodPixels - 1);
    }
}
