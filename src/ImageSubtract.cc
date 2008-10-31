// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Implementation of image subtraction functions declared in ImageSubtract.h
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup diffim
 */
#include <iostream>
#include <limits>
#include <boost/timer.hpp> 

#include <vw/Math/Functions.h> 
#include <vw/Math/Vector.h> 
#include <vw/Math/Matrix.h> 
#include <vw/Math/LinearAlgebra.h> 

#define LSST_MAX_TRACE 5                // NOTE -  trace statements >= 6 can ENTIRELY kill the run time
#include <lsst/afw/image.h>
#include <lsst/afw/math.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/logging/Log.h>
#include "lsst/ip/diffim/Pca.h"
#include "lsst/ip/diffim/ImageSubtract.h"

#define DEBUG_MATRIX 0

using lsst::pex::logging::Log;
using lsst::pex::logging::Rec;

/** 
 * @brief Runs Detection on a single image for significant peaks, and checks
 * returned Footprints for Masked pixels.
 *
 * Accepts two MaskedImages, one of which is to be convolved to match the
 * other.  The Detection package is run on the image to be convolved
 * (assumed to be higher S/N than the other image).  The subimages
 * associated with each returned Footprint in both images are checked for
 * Masked pixels; Footprints containing Masked pixels are rejected.  The
 * Footprints are grown by an amount specified in the Policy.  The
 * acceptible Footprints are returned in a vector.
 *
 * @return Vector of "clean" Footprints around which Image Subtraction
 * Kernels will be built.
 *
 * @ingroup diffim
 */
template <typename ImageT, typename MaskT>
std::vector<lsst::detection::Footprint::PtrType> lsst::ip::diffim::getCollectionOfFootprintsForPsfMatching(
    lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve, ///< Template image; is convolved
    lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve, ///< Science image; is not convolved
    lsst::pex::policy::Policy &policy ///< Policy directing the behavior
    ) {
    
    /* Parse the Policy */
    unsigned int footprintDiffimNpixMin = policy.getInt("getCollectionOfFootprintsForPsfMatching.footprintDiffimNpixMin");
    unsigned int footprintDiffimGrow = policy.getInt("getCollectionOfFootprintsForPsfMatching.footprintDiffimGrow");
    int minimumCleanFootprints = policy.getInt("getCollectionOfFootprintsForPsfMatching.minimumCleanFootprints");
    double footprintDetectionThreshold = policy.getDouble("getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold");
    double detectionThresholdScaling = policy.getDouble("getCollectionOfFootprintsForPsfMatching.detectionThresholdScaling");
    double minimumDetectionThreshold = policy.getDouble("getCollectionOfFootprintsForPsfMatching.minimumDetectionThreshold");

    /* grab mask bits from the image to convolve, since that is what we'll be operating on */
    int badMaskBit = imageToConvolve.getMask()->getMaskPlane("BAD");
    MaskT badPixelMask = (badMaskBit < 0) ? 0 : (1 << badMaskBit);

    /* Reusable view of each Footprint */
    typename lsst::afw::image::MaskedImage<ImageT, MaskT>::MaskedImagePtrT imageToConvolveFootprintPtr;
    typename lsst::afw::image::MaskedImage<ImageT, MaskT>::MaskedImagePtrT imageToNotConvolveFootprintPtr;

    /* Reusable list of Footprints */
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
            vw::BBox2i footprintBBox = (*i)->getBBox();
            footprintBBox.grow(footprintBBox.max() + vw::Vector2i(footprintDiffimGrow,footprintDiffimGrow));
            footprintBBox.grow(footprintBBox.min() - vw::Vector2i(footprintDiffimGrow,footprintDiffimGrow));

            /* grab a subimage; there is an exception if its e.g. too close to the image */
            try {
                imageToConvolveFootprintPtr = imageToConvolve.getSubImage(footprintBBox);
                imageToNotConvolveFootprintPtr = imageToNotConvolve.getSubImage(footprintBBox);
            } catch (lsst::pex::exceptions::ExceptionStack &e) {
                continue;
            }
            
            if (lsst::ip::diffim::maskOk(*(imageToConvolveFootprintPtr->getMask()), badPixelMask) && 
                lsst::ip::diffim::maskOk(*(imageToNotConvolveFootprintPtr->getMask()), badPixelMask) ) {

                /* Create a new footprint with grow'd box */
                lsst::detection::Footprint::PtrType fpGrow(new lsst::detection::Footprint(footprintBBox));
                footprintListOut.push_back(fpGrow);
                
                nCleanFootprints += 1;
            }
        }
        
        footprintDetectionThreshold *= detectionThresholdScaling;
    }
    lsst::pex::logging::TTrace<3>("lsst.ip.diffim.getCollectionOfFootprintsForPsfMatching", 
                                "Found %d clean footprints above threshold %.3f",
                                footprintListOut.size(), footprintDetectionThreshold/detectionThresholdScaling);

    return footprintListOut;
}

/** 
 * \brief Computes a single Kernel (Model 1) around a single subimage.
 *
 * Accepts two MaskedImages, generally subimages of a larger image, one of
 * which is to be convolved to match the other.  The output Kernel is
 * generated using an input list of basis Kernels by finding the
 * coefficients in front of each basis.
 *
 * \return Vector of coefficients representing the relative contribution of
 * each input basis function.
 *
 * \return Differential background offset between the two images
 *
 * \ingroup diffim
 */
template <typename ImageT, typename MaskT>
std::vector<double> lsst::ip::diffim::computePsfMatchingKernelForPostageStamp(
    double &background, ///< Difference in the backgrounds
    lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToConvolve, ///< Image to convolve
    lsst::afw::image::MaskedImage<ImageT, MaskT> const &imageToNotConvolve, ///< Image to not convolve
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelInBasisList, ///< Input kernel basis set
    lsst::pex::policy::Policy &policy ///< Policy directing the behavior
    ) { 
    
    /* grab mask bits from the image to convolve, since that is what we'll be operating on */
    int edgeMaskBit = imageToConvolve.getMask()->getMaskPlane("EDGE");
    
    int nKernelParameters = 0;
    int nBackgroundParameters = 0;
    int nParameters = 0;

    boost::timer t;
    double time;
    t.restart();

    lsst::pex::logging::TTrace<3>("lsst.ip.diffim.computePsfMatchingKernelForPostageStamp", 
                                "Entering subroutine computePsfMatchingKernelForPostageStamp");
    
    /* We assume that each kernel in the Set has 1 parameter you fit for */
    nKernelParameters = kernelInBasisList.size();
    /* Or, we just assume that across a single kernel, background 0th order.  This quite makes sense. */
    nBackgroundParameters = 1;
    /* Total number of parameters */
    nParameters = nKernelParameters + nBackgroundParameters;
    
    vw::math::Vector<double> B(nParameters);
    vw::math::Matrix<double> M(nParameters, nParameters);
    for (unsigned int i = nParameters; i--;) {
        B(i) = 0;
        for (unsigned int j = nParameters; j--;) {
            M(i,j) = 0;
        }
    }
    
    /* convolve creates a MaskedImage, push it onto the back of the Vector */
    /* need to use shared pointers because MaskedImage copy does not work */
    std::vector<boost::shared_ptr<lsst::afw::image::MaskedImage<ImageT, MaskT> > > convolvedImageList(nKernelParameters);
    /* and an iterator over this */
    typename std::vector<boost::shared_ptr<lsst::afw::image::MaskedImage<ImageT, MaskT> > >::iterator citer = convolvedImageList.begin();
    
    /* Iterator for input kernel basis */
    std::vector<boost::shared_ptr<lsst::afw::math::Kernel> >::const_iterator kiter = kernelInBasisList.begin();
    /* Create C_ij in the formalism of Alard & Lupton */
    for (; kiter != kernelInBasisList.end(); ++kiter, ++citer) {
        
        lsst::pex::logging::TTrace<6>("lsst.ip.diffim.computePsfMatchingKernelForPostageStamp",
                                    "Convolving an Object with Basis");
        
        /* NOTE : we could also *precompute* the entire template image convolved with these functions */
        /*        and save them somewhere to avoid this step each time.  however, our paradigm is to */
        /*        compute whatever is needed on the fly.  hence this step here. */
        boost::shared_ptr<lsst::afw::image::MaskedImage<ImageT, MaskT> > imagePtr(
            new lsst::afw::image::MaskedImage<ImageT, MaskT>
            (lsst::afw::math::convolveNew(imageToConvolve, **kiter, edgeMaskBit, false))
            );
        
        lsst::pex::logging::TTrace<6>("lsst.ip.diffim.computePsfMatchingKernelForPostageStamp",
                                    "Convolved an Object with Basis");
        
        *citer = imagePtr;
        
    } 
    
    /* NOTE ABOUT CONVOLUTION : */
    /* getCtrCol:getCtrRow pixels are masked on the left:bottom side */
    /* getCols()-getCtrCol():getRows()-getCtrRow() masked on right/top side */
    /* */
    /* The convolved image and the input image are by default the same size, so */
    /* we offset our initial pixel references by the same amount */
    kiter = kernelInBasisList.begin();
    citer = convolvedImageList.begin();
    unsigned int startCol = (*kiter)->getCtrCol();
    unsigned int startRow = (*kiter)->getCtrRow();
    /* NOTE - I determined I needed this +1 by eye */
    unsigned int endCol   = (*citer)->getCols() - ((*kiter)->getCols() - (*kiter)->getCtrCol()) + 1;
    unsigned int endRow   = (*citer)->getRows() - ((*kiter)->getRows() - (*kiter)->getCtrRow()) + 1;
    /* NOTE - we need to enforce that the input images are large enough */
    /* How about some multiple of the PSF FWHM?  Or second moments? */
    
    /* An accessor for each convolution plane */
    /* NOTE : MaskedPixelAccessor has no empty constructor, therefore we need to push_back() */
    std::vector<lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> > convolvedAccessorRowList;
    for (citer = convolvedImageList.begin(); citer != convolvedImageList.end(); ++citer) {
        convolvedAccessorRowList.push_back(lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT>(**citer));
    }
    
    /* An accessor for each input image; address rows and cols separately */
    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> imageToConvolveRow(imageToConvolve);
    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> imageToNotConvolveRow(imageToNotConvolve);
    
    /* Take into account buffer for kernel images */
    imageToConvolveRow.advance(startCol, startRow);
    imageToNotConvolveRow.advance(startCol, startRow);
    for (int ki = 0; ki < nKernelParameters; ++ki) {
        convolvedAccessorRowList[ki].advance(startCol, startRow);
    }

    for (unsigned int row = startRow; row < endRow; ++row) {
        
        /* An accessor for each convolution plane */
        std::vector<lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> > convolvedAccessorColList = convolvedAccessorRowList;
        
        /* An accessor for each input image; places the col accessor at the correct row */
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> imageToConvolveCol = imageToConvolveRow;
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> imageToNotConvolveCol = imageToNotConvolveRow;
        
        for (unsigned int col = startCol; col < endCol; ++col) {
            
            ImageT ncCamera   = *imageToNotConvolveCol.image;
            ImageT ncVariance = *imageToNotConvolveCol.variance;
            MaskT  ncMask     = *imageToNotConvolveCol.mask;
            
            ImageT cVariance  = *imageToConvolveCol.variance;
            
            /* Variance for a particlar pixel; do we use this variance of the */
            /* input data, or include the variance after its been convolved with */
            /* the basis?  For now, use the average of the input varianes. */
            ImageT iVariance  = 1.0 / (cVariance + ncVariance);
            
            lsst::pex::logging::TTrace<7>("lsst.ip.diffim.computePsfMatchingKernelForPostageStamp",
                                        "Accessing image row %d col %d : %.3f %.3f %d",
                                        row, col, ncCamera, ncVariance, ncMask);
            
            /* kernel index i */
            typename std::vector<lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> >::iterator
                kiteri = convolvedAccessorColList.begin();
            
            for (int kidxi = 0; kiteri != convolvedAccessorColList.end(); ++kiteri, ++kidxi) {
                ImageT cdCamerai   = *kiteri->image;
                
                ImageT cdVariancei = *kiteri->variance;
                MaskT  cdMaski     = *kiteri->mask;
                lsst::pex::logging::TTrace<7>("lsst.ip.diffim.computePsfMatchingKernelForPostageStamp",
                                            "Accessing convolved image %d : %.3f %.3f %d",
                                            kidxi, cdCamerai, cdVariancei, cdMaski);
                
                /* kernel index j  */
                typename std::vector<lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> >::iterator kiterj = kiteri;
                for (int kidxj = kidxi; kiterj != convolvedAccessorColList.end(); ++kiterj, ++kidxj) {
                    ImageT cdCameraj   = *kiterj->image;
                    
                    /* NOTE - These inner trace statements can ENTIRELY kill the run time */
                    ImageT cdVariancej = *kiterj->variance;
                    MaskT  cdMaskj     = *kiterj->mask;
                    lsst::pex::logging::TTrace<7>("lsst.ip.diffim.computePsfMatchingKernelForPostageStamp",
                                                "Accessing convolved image %d : %.3f %.3f %d",
                                                kidxj, cdCameraj, cdVariancej, cdMaskj);
                    
                    M[kidxi][kidxj] += cdCamerai * cdCameraj * iVariance;
                } 
                
                B[kidxi] += ncCamera * cdCamerai * iVariance;
                
                /* Constant background term; effectively j=kidxj+1 */
                M[kidxi][nParameters-1] += cdCamerai * iVariance;
            } 
            
            /* Background term; effectively i=kidxi+1 */
            B[nParameters-1] += ncCamera * iVariance;
            M[nParameters-1][nParameters-1] += 1.0 * iVariance;
            
            lsst::pex::logging::TTrace<6>("lsst.ip.diffim.computePsfMatchingKernelForPostageStamp", 
                                        "Background terms : %.3f %.3f",
                                        B[nParameters-1], M[nParameters-1][nParameters-1]);
            
            /* Step each accessor in column */
            imageToConvolveCol.nextCol();
            imageToNotConvolveCol.nextCol();
            for (int ki = 0; ki < nKernelParameters; ++ki) {
                convolvedAccessorColList[ki].nextCol();
            }             
            
        } /* col */
        
        /* Step each accessor in row */
        imageToConvolveRow.nextRow();
        imageToNotConvolveRow.nextRow();
        for (int ki = 0; ki < nKernelParameters; ++ki) {
            convolvedAccessorRowList[ki].nextRow();
        }
        
    } /* row */

    /* NOTE: If we are going to regularize the solution to M, this is the place to do it */
    
    /* Fill in rest of M */
    for (int kidxi=0; kidxi < nParameters; ++kidxi) 
        for (int kidxj=kidxi+1; kidxj < nParameters; ++kidxj) 
            M[kidxj][kidxi] = M[kidxi][kidxj];
    
#if DEBUG_MATRIX
    std::cout << "B : " << B << std::endl;
    std::cout << "M : " << M << std::endl;
#endif

    time = t.elapsed();
    lsst::pex::logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForPostageStamp", 
                                "Total compute time before matrix inversions : %.2f s", time);

    /* Invert using SVD and Pseudoinverse */
    vw::math::Matrix<double> Minv;
    Minv = vw::math::pseudoinverse(M);
    /*Minv = vw::math::inverse(M); */
    
#if DEBUG_MATRIX
    std::cout << "Minv : " << Minv << std::endl;
#endif
    
    /* Solve for x in Mx = B */
    vw::math::Vector<double> Soln = Minv * B;
    
#if DEBUG_MATRIX
    std::cout << "Solution : " << Soln << std::endl;
#endif
    
    time = t.elapsed();
    lsst::pex::logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForPostageStamp", 
                                "Total compute time after matrix inversions : %.2f s", time);

    /* Translate from VW std::vectors to std std::vectors */
    std::vector<double> kernelCoeffs(kernelInBasisList.size());
    for (int ki = 0; ki < nKernelParameters; ++ki) {
        kernelCoeffs[ki] = Soln[ki];
    }
    background = Soln[nParameters-1];

    lsst::pex::logging::TTrace<3>("lsst.ip.diffim.computePsfMatchingKernelForPostageStamp", 
                                "Leaving subroutine computePsfMatchingKernelForPostageStamp");

    return kernelCoeffs;
}

/** 
 * @brief Computes a orthonormal Kernel basis from an ensemble of input
 * Kernels.
 *
 * Accepts a vector of input DiffImContainers.  The good Kernels associated
 * with each container are used to build an orthornormal basis using a
 * Karhunen Loeve transform, implemented using Principal Component Analysis.
 *
 * Each input Kernel is decomposed using a subset of these basis Kernels
 * (yielding Kernel Model 2).  This approximate Kernel and associated
 * difference image statistics are added to the container.  The containers
 * associated with poor subtractions have the isGood flag set to False.  The
 * process iterates until no containers are rejected.  The final basis set
 * is returned as a vector of Kernels.
 *
 * @return Vector of eigenKernels resulting from the PCA analysis.
 *
 * @throw lsst::pex::exceptions::Runtime if no good kernels
 *
 * @ingroup diffim
 */
template <typename ImageT, typename MaskT>
lsst::afw::math::KernelList<lsst::afw::math::Kernel> lsst::ip::diffim::computePcaKernelBasis(
    std::vector<lsst::ip::diffim::DiffImContainer<ImageT, MaskT> > &diffImContainerList, ///< List of input footprints
    lsst::pex::policy::Policy &policy ///< Policy directing the behavior
    ) {
    
    /* Parse the Policy */
    unsigned int minimumNumberOfBases = policy.getInt("computePcaKernelBasis.minimumNumberOfBases");
    unsigned int maximumNumberOfBases = policy.getInt("computePcaKernelBasis.maximumNumberOfBases");
    double maximumFractionOfEigenvalues = policy.getDouble("computePcaKernelBasis.maximumFractionOfEigenvalues");
    double minimumAcceptibleEigenvalue = policy.getDouble("computePcaKernelBasis.minimumAcceptibleEigenvalue");
    unsigned int maximumIteratonsPCA = policy.getInt("computePcaKernelBasis.maximumIteratonsPCA");
    /*double maximumKernelResidualMean = policy.getDouble("computePcaKernelBasis.maximumKernelResidualMean"); */
    /*double maximumKernelResidualVariance = policy.getDouble("computePcaKernelBasis.maximumKernelResidualVariance"); */
    bool debugIO = policy.getBool("debugIO", false);

    double imSum;
    
    /* Image accessor */
    typedef typename vw::ImageView<lsst::afw::math::Kernel::PixelT>::pixel_accessor imageAccessorType;
    /* Iterator over struct */
    typedef typename std::vector<lsst::ip::diffim::DiffImContainer<ImageT, MaskT> >::iterator iDiffImContainer;
    
    /* Matrix to invert.  Number of rows = number of pixels; number of columns = number of kernels */
    /* All calculations here are in double */
    /* Assigment of matrix sizes is rows, cols */
    vw::math::Matrix<double> M;
    vw::math::Matrix<double> eVec;
    vw::math::Vector<double> eVal;
    vw::math::Vector<double> mMean;
    vw::math::Matrix<double> kernelCoefficientMatrix;

    Log pcaLog(Log::getDefaultLog(), "ip.diffim.computePcaKernelBasis");

    /* Initialize for while loop */
    unsigned int nIter = 0;
    unsigned int nReject = 1;
    unsigned int nGood = 0;
    unsigned int nKCols = 0, nKRows = 0;
    unsigned int nPixels;
    unsigned int nCoeffToUse;
    lsst::pex::logging::TTrace<3>("lsst.ip.diffim.computePcaKernelBasis", " ");
    lsst::pex::logging::TTrace<3>("lsst.ip.diffim.computePcaKernelBasis", 
                                "Entering subroutine computePcaKernelBasis");
    
    /* Reuse */
    /* NOTE - this assumes that all kernels are the same size */
    /*        and that this size is defined in the policy */ 
    unsigned int kernelRows = policy.getInt("kernelRows");
    unsigned int kernelCols = policy.getInt("kernelCols");
    lsst::afw::image::Image<lsst::afw::math::Kernel::PixelT> kImage(kernelCols, kernelRows);

    lsst::afw::math::KernelList<lsst::afw::math::Kernel> kernelPcaBasisList;
    
    /* Iterate over PCA inputs until all are good */
    while ( (nIter < maximumIteratonsPCA) and (nReject != 0) ) {
        
        nGood = 0;
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
            if ((*i).isGood == true) {
                nGood += 1;
                if (nKCols == 0) {
                    nKCols = (*i).kernelList[(*i).nKernel]->getCols();
                }
                if (nKRows == 0) {
                    nKRows = (*i).kernelList[(*i).nKernel]->getRows();
                }
            }
        }
        
        if (nGood == 0) {
            throw lsst::pex::exceptions::DomainError("No good kernels for PCA");
        }
        else {
            lsst::pex::logging::TTrace<4>("lsst.ip.diffim.computePcaKernelBasis", 
                                        "PCA Iteration %d : Using %d kernels", nIter, nGood);

            pcaLog.log(Log::INFO, boost::format("PCA Iteration %d : Using %d kernels") % nIter % nGood);
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
        
        /* fill up matrix for PCA */
        int ki = 0;
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
            if ((*i).isGood == false) {
                continue;
            }

            /* Note we ALWAYS compute the PCA using the first kernel image */
            (*i).kernelList[0]->computeImage(kImage, imSum, false);

            /* insanity checking */
            assert(nKRows == kImage.getRows());
            assert(nKCols == kImage.getCols());
            
            int mIdx = 0;
            imageAccessorType imageAccessorCol(kImage.origin());
            for (unsigned int col = 0; col < nKCols; ++col, imageAccessorCol.next_col()) {
                
                imageAccessorType imageAccessorRow(imageAccessorCol);
                for (unsigned int row = 0; row < nKRows; ++row, ++mIdx, imageAccessorRow.next_row()) {

                    /*
                      NOTE : 
                        arguments to matrix-related functions are given in row,col
                        arguments to image-related functions are given in col,row
                      HOWEVER :
                        it doesn't matter which order we put the kernel elements into the PCA
                        since they are not correlated 
                        as long as we do the same thing when extracting the components
                      UNLESS :
                        we want to put in some weighting/regularlization into the PCA
                        not sure if that is even possible...
                    */

                    M(mIdx, ki) = *imageAccessorRow;
                }
            }
            
            /* only increment when you have a good footprint */
            ki += 1;
        }
        
       
        /* M is mean-subtracted if subtractMean = true */
        /* Eigencomponents in columns of eVec */
        lsst::pex::logging::TTrace<4>("lsst.ip.diffim.computePcaKernelBasis", "Computing pricipal components");
        lsst::ip::diffim::computePca(mMean, eVal, eVec, M, true);
        lsst::pex::logging::TTrace<4>("lsst.ip.diffim.computePcaKernelBasis", "Computed pricipal components");

        /* This is probably as inefficient as you can get, but I'm no boost expert */
        boost::format eValueFormatter("Eigenvalues :");
        for (unsigned int i = 0; i < eVal.size(); ++i) {
            eValueFormatter = boost::format("%s %.3f") % eValueFormatter.str() % eVal[i];
        }
        lsst::pex::logging::TTrace<4>("lsst.ip.diffim.computePcaKernelBasis", eValueFormatter.str());
        pcaLog.log(Log::INFO, eValueFormatter.str());
        
        nCoeffToUse = minimumNumberOfBases;
        /* Potentially override with larger number if the spectrum of eigenvalues requires it */
        double evalSum = 0.0;
        for (unsigned int i = 0; i < eVal.size(); ++i) {
            evalSum += eVal[i];
        }
        double evalFrac = 0.0;
        for (unsigned int i = 0; i < nCoeffToUse; ++i) {
            evalFrac += eVal[i] / evalSum;
            if (eVal[i] < minimumAcceptibleEigenvalue) {
                lsst::pex::logging::TTrace<1>("lsst.ip.diffim.computePcaKernelBasis", 
                                            "WARNING : Using eigenvector whose eigenvalue (%.3e) is smaller than acceptible (%.3e)",
                                            eVal[i], minimumAcceptibleEigenvalue);
            }
            
        }
        if (evalFrac < maximumFractionOfEigenvalues) {
            for (; ( (nCoeffToUse < eVal.size()) && (nCoeffToUse < maximumNumberOfBases) ); ++nCoeffToUse) {
                evalFrac += eVal[nCoeffToUse] / evalSum;
                if (eVal[nCoeffToUse] < minimumAcceptibleEigenvalue) {
                    lsst::pex::logging::TTrace<1>("lsst.ip.diffim.computePcaKernelBasis", 
                                                "WARNING : Using eigenvector whose eigenvalue (%.3e) is smaller than acceptible (%.3e)",
                                                eVal[nCoeffToUse], minimumAcceptibleEigenvalue);
                }
                if (evalFrac > maximumFractionOfEigenvalues) {
                    break;
                }
            }
        }
        lsst::pex::logging::TTrace<4>("lsst.ip.diffim.computePcaKernelBasis", 
                                    "Using %d basis functions (plus mean) representing %.5f of the variance",
                                    nCoeffToUse, evalFrac);
        pcaLog.log(Log::INFO, boost::format("Using %d basis functions (plus mean) representing %.5f of the variance") % nCoeffToUse % evalFrac);
        
        /* We now have the basis functions determined */
        /* Next determine the coefficients that go in front of all of the individual kernels */

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
        lsst::ip::diffim::decomposeMatrixUsingBasis(kernelCoefficientMatrix, M, eVec, nCoeffToUse);
        
        /* We next do quality control here; we reconstruct the input kernels with the truncated basis function set */
        /* Remember that M is mean-subtracted already */
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
        lsst::ip::diffim::approximateMatrixUsingBasis(approxM, eVec, kernelCoefficientMatrix, nCoeffToUse);
        
#if DEBUG_MATRIX
        std::cout << "EigenKernel coefficients : " << kernelCoefficientMatrix << std::endl;
#endif

        /* Reset basis list */
        kernelPcaBasisList.clear();
        
        /* Turn the Mean Image into a Kernel */
        lsst::afw::image::Image<lsst::afw::math::Kernel::PixelT> meanImage(nKCols, nKRows);
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
        /* The mean image is the first member of kernelPCABasisList */
        kernelPcaBasisList.push_back(boost::shared_ptr<lsst::afw::math::Kernel> 
                                     (new lsst::afw::math::FixedKernel(meanImage)));
        if (debugIO) {
            meanImage.writeFits( (boost::format("mFits%d.fits") % nIter).str());
        }
        
        /* Turn each eVec into an Image and then into a Kernel */
        /* Dont use all eVec.cols(); only the number of bases that you want */
        for (unsigned int ki = 0; ki < nCoeffToUse; ++ki) { 
            lsst::afw::image::Image<lsst::afw::math::Kernel::PixelT> basisImage(nKCols, nKRows);
            
            /* Not sure how to bulk load information into Image */
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
            /* Add to kernel basis */
            kernelPcaBasisList.push_back(boost::shared_ptr<lsst::afw::math::Kernel> 
                                         (new lsst::afw::math::FixedKernel(basisImage)));
            if (debugIO) {
                basisImage.writeFits( (boost::format("eFits%d_%d.fits") % nIter % ki).str() );
            }
        }
        
        nReject = 0;
        unsigned int iKernel = 0;
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
            if ((*i).isGood == false) {
                continue;
            }

            /* Calculate kernel from PCA model */
            std::vector<double> kernelCoefficients;
            
            /* Mean image */
            kernelCoefficients.push_back(1); 
            for (unsigned int jKernel = 0; jKernel < nCoeffToUse; ++jKernel) {
                kernelCoefficients.push_back(kernelCoefficientMatrix(iKernel,jKernel));
            }        
            
            /* Create a linear combination kernel from this  */
            boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> pcaKernelPtr(
                new lsst::afw::math::LinearCombinationKernel(kernelPcaBasisList, kernelCoefficients)
                );
            
            /* Append to vectors collectively */
            /* NOTE HERE - You will get one kernel for each iteration */
            (*i).kernelList.push_back( pcaKernelPtr );
            (*i).diffimStats.push_back( lsst::ip::diffim::MaskedImageDiffimStats() );
            (*i).backgrounds.push_back( (*i).backgrounds[(*i).nKernel] );   /* This did not change here */
            pcaKernelPtr->computeImage(kImage, imSum, false);
            (*i).kernelSums.push_back( imSum );
            (*i).nKernel += 1;

            /* Calculate image stats */
            lsst::ip::diffim::computeDiffImStats((*i), (*i).nKernel, policy);
            if ( (*i).isGood == false ) {
                nReject += 1;
            }
            
            /* only increment for the good footprints */
            iKernel += 1;
            
        }

        nIter += 1;
    } /* End of iteration */
    
    
    return kernelPcaBasisList;
}

/** 
 * @brief Compute spatially varying convolution Kernel and differential
 * background model.
 *
 * Accepts a vector of input DiffImContainers, as well as Functions for the
 * spatial variation of each basis Kernel and the differential Background.
 *
 * The Kernels associated with each container are assumed to be
 * LinearCombination Kernels built from the same basis set.  The
 * coefficients for each basis Kernel are used to fit a spatially varying
 * Function for that Kernel.  A Function is also fit to the background.  The
 * spatially varying Kernel is returned.
 *
 * @return Spatially varying convolution Kernel.
 *
 * @ingroup diffim
 */
template <typename ImageT, typename MaskT>
boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> lsst::ip::diffim::computeSpatiallyVaryingPsfMatchingKernel(
    lsst::afw::math::Function2<double> &kernelFunction, ///< Function for spatial variation of kernel
    lsst::afw::math::Function2<double> &backgroundFunction, ///< Function for spatial variation of background
    std::vector<lsst::ip::diffim::DiffImContainer<ImageT, MaskT> > &diffImContainerList, ///< Information on each kernel
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelBasisList, ///< Basis kernel set
    lsst::pex::policy::Policy &policy ///< Policy directing the behavior
    )
{
    /* Parse the Policy */
    unsigned int maximumIterationsSpatialFit = policy.getInt("computeSpatiallyVaryingPsfMatchingKernel.maximumIterationsSpatialFit");
    
    /* Container iterator */
    typedef typename std::vector<lsst::ip::diffim::DiffImContainer<ImageT, MaskT> >::iterator iDiffImContainer;
    /* Kernel iterator */
    typedef typename std::vector<boost::shared_ptr<lsst::afw::math::Kernel> >::const_iterator iKernelPtr;

    Log spatialLog(Log::getDefaultLog(), "ip.diffim.computeSpatiallyVaryingPsfMatchingKernel");

    /* Initialize for while loop */
    unsigned int nIter = 0;
    unsigned int nReject = 1;
    unsigned int nGood = diffImContainerList.size();
    unsigned int nKernel = kernelBasisList.size();
    unsigned int nKernelParameters = kernelFunction.getNParameters();
    unsigned int nBgParameters = backgroundFunction.getNParameters();
    
    lsst::pex::logging::TTrace<3>("lsst.ip.diffim.computeSpatiallyVaryingPsfMatchingKernel", 
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

    double bgSum = 0.0;
    double imSum = 0.0;
    
    /* NOTE - this assumes that all kernels are the same size */
    /*        and that this size is defined in the policy */ 
    unsigned int kernelRows = policy.getInt("kernelRows");
    unsigned int kernelCols = policy.getInt("kernelCols");
    lsst::afw::image::Image<lsst::afw::math::Kernel::PixelT> kImage(kernelCols, kernelRows);
    lsst::afw::image::Image<lsst::afw::math::Kernel::PixelT> kSpatialImage(kernelCols, kernelRows);
    
    /* Set up the spatially varying kernel */
    boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> spatiallyVaryingKernelPtr(
        new lsst::afw::math::LinearCombinationKernel(kernelBasisList, kernelFunction));
    
    /* Iterate over inputs until both the kernel model and the background model converge */
    while ( (nIter < maximumIterationsSpatialFit) and (nReject != 0) ) {
        /* Start out clean, nothing rejected yet */
        nReject = 0;

        /* Set up fit to kernel */
        kernelColPos.clear();
        kernelRowPos.clear();
        nGood = 0;

        backgroundMeas.clear();
        backgroundVariances.clear(); 
        bgSum = 0.0;
        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
            /* those attributs that you grab once for each container; position and background */
            if ((*i).isGood == true) {
                kernelColPos.push_back((*i).colcNorm);
                kernelRowPos.push_back((*i).rowcNorm);

                bgSum += (*i).backgrounds[(*i).nKernel];
                backgroundMeas.push_back((*i).backgrounds[(*i).nKernel]);
                backgroundVariances.push_back((*i).diffimStats[(*i).nKernel].footprintResidualVariance);   /* approximation */

                lsst::pex::logging::TTrace<5>("lsst.ip.diffim.computePsfMatchingKernelForMaskedImage", 
                                            "Background %d at %.3f %.3f = %.3f",
                                            (*i).id, (*i).colcNorm, (*i).rowcNorm, (*i).backgrounds[(*i).nKernel]);

                nGood += 1;
            }
        }
        
        /* NOTE - Sigma clip background here, same as kernel sum? */

        if (nGood == 0) {
            throw lsst::pex::exceptions::DomainError("No good footprints for spatial kernel and background fitting");
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
                    std::vector<double> kernelParameters = (*id).kernelList[(*id).nKernel]->getKernelParameters(); 
                    kernelMeas.push_back(kernelParameters[nk]);

                    kernelVariances.push_back((*id).diffimStats[(*id).nKernel].footprintResidualVariance); /* approximation */
                    spatialFitFormatter = boost::format("%s %.3f") % spatialFitFormatter.str() % kernelParameters[nk];
                }
            }
            lsst::pex::logging::TTrace<6>("lsst.ip.diffim.computeSpatiallyVaryingPsfMatchingKernel",
                                        spatialFitFormatter.str());
            
            /* NOTE - if we have fewer measurements than kernelFunction.getNParameters(), we should do something about it */
            std::vector<double> kernelParams(nKernelParameters);
            std::vector<double> kernelStepsize(nKernelParameters);
            std::fill(kernelParams.begin(), kernelParams.end(), 0.0);         /* start at value zero */
            std::fill(kernelStepsize.begin(), kernelStepsize.end(), 0.02);  /* assume order 2% contribution from eigenkernels */
            
            /* Run minimization */
            lsst::afw::math::FitResults kernelFit = lsst::afw::math::minimize(
                kernelFunction,
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
                lsst::pex::logging::TTrace<4>("lsst.ip.diffim.computeSpatiallyVaryingPsfMatchingKernel", 
                                            "Fit to kernel PC%d, spatial parameter %d : %.3f (%.3f,%.3f)",
                                            nk, i, kernelFit.parameterList[i], kernelFit.parameterErrorList[i].first,
                                            kernelFit.parameterErrorList[i].second);
                spatialLog.log(Log::INFO,
                               boost::format("Iteration %d, Fit to Kernel PC%d, Spatial Parameter %d : %.3f (%.3f,%.3f)") %
                               nIter % nk % i % 
                               kernelFit.parameterList[i] % 
                               kernelFit.parameterErrorList[i].first % 
                               kernelFit.parameterErrorList[i].second);
            }
        }
        /* After all kernels are fit, fill the spatially varying kernel parameters */
        spatiallyVaryingKernelPtr->setSpatialParameters(fitParameters);

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
        lsst::afw::math::FitResults backgroundFit = lsst::afw::math::minimize(
            backgroundFunction,
            backgroundParameters,
            backgroundStepsize,
            backgroundMeas,
            backgroundVariances,
            kernelColPos,
            kernelRowPos,
            nSigmaSq
            );
        backgroundFunction.setParameters(backgroundFit.parameterList);
        
        /* Debugging information */
        for (unsigned int i = 0; i < backgroundFit.parameterList.size(); ++i) {
            lsst::pex::logging::TTrace<4>("lsst.ip.diffim.computePsfMatchingKernelForMaskedImage", 
                                        "Fit Background Parameter %d : %.3f (%.3f,%.3f)\n",
                                        i, backgroundFit.parameterList[i], backgroundFit.parameterErrorList[i].first,
                                        backgroundFit.parameterErrorList[i].second);
            spatialLog.log(Log::INFO,
                           boost::format("Iteration %d, Fit to Background Parameter %d : %.3f (%.3f,%.3f)") %
                           nIter % i % backgroundFit.parameterList[i] %
                           backgroundFit.parameterErrorList[i].first %
                           backgroundFit.parameterErrorList[i].second);
        }
        
        /* Fit for spatially varying background */
        /* ******************************** */
        /* ******************************** */
        /* Check each footprint's actual kernel with the spatial model */

        for (iDiffImContainer i = diffImContainerList.begin(); i != diffImContainerList.end(); ++i) {
            if ((*i).isGood == false) {
                continue;
            }
          
            (*i).kernelList.push_back( spatiallyVaryingKernelPtr );
            (*i).diffimStats.push_back( lsst::ip::diffim::MaskedImageDiffimStats() );
            (*i).backgrounds.push_back( backgroundFunction((*i).colcNorm, (*i).rowcNorm) );
            spatiallyVaryingKernelPtr->computeImage(kSpatialImage, imSum, false, (*i).colcNorm, (*i).rowcNorm);
            (*i).kernelSums.push_back( imSum );
            (*i).nKernel += 1;

            /* Calculate image stats */
            lsst::ip::diffim::computeDiffImStats((*i), (*i).nKernel, policy);
            if ( (*i).isGood == false ) {
                nReject += 1;
            }
        }

        /* End of the loop */
        nIter += 1;

    }
    return spatiallyVaryingKernelPtr;
}

template <typename ImageT, typename MaskT>
void lsst::ip::diffim::computeDiffImStats(
    lsst::ip::diffim::DiffImContainer<ImageT, MaskT> &diffImContainer,
    int const kernelID,
    lsst::pex::policy::Policy &policy
    ) {

    double imSum;
    double meanOfResiduals = 0.0;
    double varianceOfResiduals = 0.0;
    int nGoodPixels = 0;

    Log statsLog(Log::getDefaultLog(), "ip.diffim.computeDiffImStats");

    double maximumFootprintResidualMean = policy.getDouble("computeDiffImStats.maximumFootprintResidualMean");
    double maximumFootprintResidualVariance = policy.getDouble("computeDiffImStats.maximumFootprintResidualVariance");

    if ( (kernelID <= -1) || (kernelID > diffImContainer.nKernel) ) {
        throw lsst::pex::exceptions::DomainError("Request for non-existent kernel");
    }

    int edgeMaskBit = diffImContainer.imageToConvolve->getMask()->getMaskPlane("EDGE");
    int badMaskBit  = diffImContainer.imageToConvolve->getMask()->getMaskPlane("BAD");
    MaskT edgePixelMask = (edgeMaskBit < 0) ? 0 : (1 << edgeMaskBit);
    MaskT badPixelMask = (badMaskBit < 0) ? 0 : (1 << badMaskBit);
    badPixelMask |= edgePixelMask;

    bool debugIO = policy.getBool("debugIO", false);
    unsigned int kernelRows = policy.getInt("kernelRows");
    unsigned int kernelCols = policy.getInt("kernelCols");
    lsst::afw::image::Image<lsst::afw::math::Kernel::PixelT> kImage(kernelCols, kernelRows);

    lsst::afw::image::MaskedImage<ImageT, MaskT>
        convolvedImageStamp = lsst::afw::math::convolveNew(
            *(diffImContainer.imageToConvolve),
            *(diffImContainer.kernelList[kernelID]),
            edgeMaskBit, 
            false);
    
    /* Add in background */
    convolvedImageStamp += diffImContainer.backgrounds[kernelID];
    if (debugIO) {
        convolvedImageStamp.writeFits( (boost::format("cFits%d_%d") % kernelID % diffImContainer.id).str() );
    }
    
    /* Do actual subtraction */
    convolvedImageStamp -= *(diffImContainer.imageToNotConvolve);
    convolvedImageStamp *= -1.0;
    if (debugIO) {
        convolvedImageStamp.writeFits( (boost::format("dFits%d_%d") % kernelID % diffImContainer.id).str() );
    }
    
    /* calculate stats in the difference image */
    lsst::ip::diffim::calculateMaskedImageResiduals(nGoodPixels, meanOfResiduals, varianceOfResiduals, convolvedImageStamp, badPixelMask);
    diffImContainer.diffimStats[kernelID].footprintResidualMean = meanOfResiduals;
    diffImContainer.diffimStats[kernelID].footprintResidualVariance = varianceOfResiduals;

    if (!(std::isfinite(diffImContainer.diffimStats[kernelID].footprintResidualMean)) || 
        !(std::isfinite(diffImContainer.diffimStats[kernelID].footprintResidualVariance))) {
        lsst::pex::logging::TTrace<2>("lsst.ip.diffim.computeDiffImStats", 
                                    "# Kernel %d (kSum=%.3f), analysis of footprint failed : %.3f %.3f (%d pixels)",
                                    diffImContainer.id,
                                    diffImContainer.kernelSums[kernelID],
                                    diffImContainer.diffimStats[kernelID].footprintResidualMean,
                                    diffImContainer.diffimStats[kernelID].footprintResidualVariance, 
                                    nGoodPixels);

        diffImContainer.isGood = false;
    }
    else if (fabs(diffImContainer.diffimStats[kernelID].footprintResidualMean) > maximumFootprintResidualMean) {
        lsst::pex::logging::TTrace<2>("lsst.ip.diffim.computeDiffImStats", 
                                    "# Kernel %d (kSum=%.3f), bad mean residual of footprint : %.3f (%d pixels)",
                                    diffImContainer.id,
                                    diffImContainer.kernelSums[kernelID],
                                    diffImContainer.diffimStats[kernelID].footprintResidualMean,
                                    nGoodPixels);

        diffImContainer.isGood = false;
    }
    else if (diffImContainer.diffimStats[kernelID].footprintResidualVariance > maximumFootprintResidualVariance) {
        lsst::pex::logging::TTrace<2>("lsst.ip.diffim.computeDiffImStats", 
                                    "# Kernel %d (kSum=%.3f), bad residual variance of footprint : %.3f (%d pixels)",
                                    diffImContainer.id,
                                    diffImContainer.kernelSums[kernelID],
                                    diffImContainer.diffimStats[kernelID].footprintResidualVariance,
                                    nGoodPixels);

        diffImContainer.isGood = false;
    }
    else {
        lsst::pex::logging::TTrace<5>("lsst.ip.diffim.computeDiffImStats", 
                                    "Kernel %d (kSum=%.3f), mean and variance of residuals in footprint : %.3f %.3f (%d pixels)",
                                    diffImContainer.id,
                                    diffImContainer.kernelSums[kernelID],
                                    diffImContainer.diffimStats[kernelID].footprintResidualMean,
                                    diffImContainer.diffimStats[kernelID].footprintResidualVariance,
                                    nGoodPixels);
    }
    if (diffImContainer.isGood == false) {
        statsLog.log(Log::INFO, 
                     boost::format("#Kernel %d (kSum=%.3f), mean and variance of residuals in footprint : %.3f %.3f (%d pixels)") %
                     diffImContainer.id %
                     diffImContainer.kernelSums[kernelID] %
                     diffImContainer.diffimStats[kernelID].footprintResidualMean %
                     diffImContainer.diffimStats[kernelID].footprintResidualVariance %
                     nGoodPixels);
    }
    else {
        statsLog.log(Log::INFO, 
                     boost::format("Kernel %d (kSum=%.3f), mean and variance of residuals in footprint : %.3f %.3f (%d pixels)") %
                     diffImContainer.id %
                     diffImContainer.kernelSums[kernelID] %
                     diffImContainer.diffimStats[kernelID].footprintResidualMean %
                     diffImContainer.diffimStats[kernelID].footprintResidualVariance %
                     nGoodPixels);
    }

    if (debugIO) {
        diffImContainer.kernelList[kernelID]->computeImage(kImage, imSum, false);
        if (diffImContainer.isGood == true) {
            kImage.writeFits( (boost::format("kFits%dG_%d.fits") % kernelID % diffImContainer.id).str() );
        }
        else {
            kImage.writeFits( (boost::format("kFits%dB_%d.fits") % kernelID % diffImContainer.id).str() );
        }
    }
}

/** 
 * @brief Generate a basis set of delta function Kernels.
 *
 * Generates a vector of Kernels sized nCols * nRows, where each Kernel has
 * a unique pixel set to value 1.0 with the other pixels valued 0.0.  This
 * is the "delta function" basis set.
 * 
 * @return Vector of orthonormal delta function Kernels.
 *
 * @throw lsst::pex::exceptions::DomainError if nRows or nCols not positive
 *
 * @ingroup diffim
 */
lsst::afw::math::KernelList<lsst::afw::math::Kernel>
lsst::ip::diffim::generateDeltaFunctionKernelSet(
    unsigned int nCols, ///< Number of colunms in the basis kernels
    unsigned int nRows  ///< Number of rows in the basis kernels
    )
{
    if ((nCols < 1) || (nRows < 1)) {
        throw lsst::pex::exceptions::DomainError("nRows and nCols must be positive");
    }
    const int signedNCols = static_cast<int>(nCols);
    const int signedNRows = static_cast<int>(nRows);
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> kernelBasisList;
    for (int row = 0; row < signedNRows; ++row) {
        for (int col = 0; col < signedNCols; ++col) {
            boost::shared_ptr<lsst::afw::math::Kernel> kernelPtr(
                new lsst::afw::math::DeltaFunctionKernel(col, row, nCols, nRows)
                );
            
            kernelBasisList.push_back(kernelPtr);
        }
    }
    return kernelBasisList;
}

/** 
 * @brief Generate an Alard-Lupton basis set of Kernels.
 *
 * Not implemented.
 * 
 * @return Vector of Alard-Lupton Kernels.
 *
 * @throw lsst::pex::exceptions::DomainError if nRows or nCols not positive
 * @throw lsst::pex::exceptions::DomainError until implemented
 *
 * @ingroup diffim
 */
lsst::afw::math::KernelList<lsst::afw::math::Kernel>
lsst::ip::diffim::generateAlardLuptonKernelSet(
    unsigned int nRows, ///< Number of rows in the kernel basis
    unsigned int nCols, ///< Number of columns in the kernel basis
    std::vector<double> const &sigGauss, ///< Width of gaussians in basis; size = number of Gaussians
    std::vector<double> const &degGauss  ///< Degree of spatial variation within each Gaussian; size = sigGauss.size()
    )
{
    if ((nCols < 1) || (nRows < 1)) {
        throw lsst::pex::exceptions::DomainError("nRows and nCols must be positive");
    }
    throw lsst::pex::exceptions::DomainError("Not implemented");

    lsst::afw::math::KernelList<lsst::afw::math::Kernel> kernelBasisList;
    return kernelBasisList;
}

/** 
 * @brief Checks a Mask image to see if a particular Mask plane is set.
 *
 * @return True if the mask is *not* set, False if it is.
 *
 * @ingroup diffim
 */
template <typename MaskT>
bool lsst::ip::diffim::maskOk(
    lsst::afw::image::Mask<MaskT> const &inputMask,
    MaskT const badPixelMask ///< Mask value for bad data
    )
{
    typename lsst::afw::image::Mask<MaskT>::pixel_accessor rowAcc = inputMask.origin();
    for (unsigned int row = 0; row < inputMask.getRows(); ++row, rowAcc.next_row()) {
        typename lsst::afw::image::Mask<MaskT>::pixel_accessor colAcc = rowAcc;
        for (unsigned int col = 0; col < inputMask.getCols(); ++col, colAcc.next_col()) {
            /*std::cout << "MASK " << (*colAcc) << " " << badPixelMask << " " << ((*colAcc) & badPixelMask) << std::endl; */
            
            if (((*colAcc) & badPixelMask) != 0) {
                return false;
            }
        }
    }
    return true;
}

/** 
 * @brief Calculates mean and variance of values (normalized by the sqrt of
 * the image variance) in a MaskedImage
 *
 * The pixel values in the image portion are normalized by the sqrt of the
 * variance.  The mean and variance of this distribution are calculated.  If
 * the MaskedImage is a difference image, the results should follow a
 * normal(0,1) distribution.
 *
 * @return Number of unmasked pixels in the image, and the mean and variance
 * of the residuals divided by the sqrt(variance).
 *
 * @ingroup diffim
 */
template <typename ImageT, typename MaskT>
void lsst::ip::diffim::calculateMaskedImageResiduals(
    int &nGoodPixels, ///< Number of good pixels in the image
    double &meanOfResiduals, ///< Mean residual/variance; ideally 0
    double &varianceOfResiduals, ///< Average variance of residual/variance; ideally 1
    lsst::afw::image::MaskedImage<ImageT, MaskT> const &inputImage, ///< Input image to be analyzed
    MaskT const badPixelMask ///< Mask for bad data
    ) {
    
    double x2Sum=0.0, xSum=0.0;
    
    nGoodPixels = 0;
    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> rowAcc(inputImage);
    for (unsigned int row = 0; row < inputImage.getRows(); ++row, rowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> colAcc = rowAcc;
        for (unsigned int col = 0; col < inputImage.getCols(); ++col, colAcc.nextCol()) {
            if (((*colAcc.mask) & badPixelMask) == 0) {
                xSum  += (*colAcc.image) / sqrt(*colAcc.variance);
                x2Sum += (*colAcc.image) * (*colAcc.image) / (*colAcc.variance);

                nGoodPixels += 1;
            }
        }
    }
    
    if (nGoodPixels > 0) {
        meanOfResiduals = xSum / nGoodPixels;
    } else {
        meanOfResiduals = std::numeric_limits<double>::quiet_NaN();
    }
    if (nGoodPixels > 1) {
        varianceOfResiduals  = x2Sum / nGoodPixels - meanOfResiduals*meanOfResiduals;
        varianceOfResiduals *= nGoodPixels / (nGoodPixels - 1);
    } else {
        varianceOfResiduals = std::numeric_limits<double>::quiet_NaN();
    }
    
}

/** 
 * @brief Calculates mean and variance of the values in an Image
 *
 * @return Number of pixels in the image, and the mean and variance of the
 * values.
 *
 * @ingroup diffim
 */
template <typename ImageT>
void lsst::ip::diffim::calculateImageResiduals(
    int &nGoodPixels, ///< Number of good pixels in the image
    double &meanOfResiduals, ///< Mean residual; nan if nGoodPixels < 1
    double &varianceOfResiduals, ///< Average variance of residual; nan if nGoodPixels < 2
    lsst::afw::image::Image<ImageT> const &inputImage ///< Input image to be analyzed
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
        varianceOfResiduals /= nGoodPixels;
    } else {
        varianceOfResiduals = std::numeric_limits<double>::quiet_NaN();
    }
}

/** 
 * @brief Calculates mean and variance of the values in a vector
 *
 * @return The mean and variance of the values.
 *
 * @ingroup diffim
 */
template <typename VectorT>
void lsst::ip::diffim::calculateVectorStatistics(
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
 * @brief Add a Function to an Image
 *
 * @ingroup diffim
 */
template <typename PixelT, typename FunctionT>
void lsst::ip::diffim::addFunction(
    lsst::afw::image::Image<PixelT> &image, ///< image
    lsst::afw::math::Function2<FunctionT> const &function ///< 2-d function
) {
    typedef typename lsst::afw::image::Image<PixelT>::pixel_accessor imageAccessorType;
    unsigned int numCols = image.getCols();
    unsigned int numRows = image.getRows();
    imageAccessorType imRow = image.origin();
    for (unsigned int row = 0; row < numRows; ++row, imRow.next_row()) {
        imageAccessorType imCol = imRow;
        double rowPos = lsst::afw::image::positionToIndex(row);
        for (unsigned int col = 0; col < numCols; ++col, imCol.next_col()) {
            double colPos = lsst::afw::image::positionToIndex(col);
            *imCol += static_cast<PixelT>(function(colPos, rowPos));
        }
    }
}

/** 
 * @brief Subtract a Function from an Image
 *
 * @ingroup diffim
 */
template <typename PixelT, typename FunctionT>
void lsst::ip::diffim::subtractFunction(
    lsst::afw::image::Image<PixelT> &image, ///< image
    lsst::afw::math::Function2<FunctionT> const &function ///< 2-d function
) {
    typedef typename lsst::afw::image::Image<PixelT>::pixel_accessor imageAccessorType;
    unsigned int numCols = image.getCols();
    unsigned int numRows = image.getRows();
    imageAccessorType imRow = image.origin();
    for (unsigned int row = 0; row < numRows; ++row, imRow.next_row()) {
        imageAccessorType imCol = imRow;
        double rowPos = lsst::afw::image::positionToIndex(row);
        for (unsigned int col = 0; col < numCols; ++col, imCol.next_col()) {
            double colPos = lsst::afw::image::positionToIndex(col);
            *imCol -= static_cast<PixelT>(function(colPos, rowPos));
        }
    }
}

/************************************************************************************************************/
/* Explicit instantiations */

template
std::vector<double> lsst::ip::diffim::computePsfMatchingKernelForPostageStamp(
    double &background,
    lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> const &imageToConvolve,
    lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> const &imageToNotConvolve,
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelInBasisList,
    lsst::pex::policy::Policy &policy);

template
std::vector<double> lsst::ip::diffim::computePsfMatchingKernelForPostageStamp(
    double &background,
    lsst::afw::image::MaskedImage<double, lsst::afw::image::maskPixelType> const &imageToConvolve,
    lsst::afw::image::MaskedImage<double, lsst::afw::image::maskPixelType> const &imageToNotConvolve,
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelInBasisList,
    lsst::pex::policy::Policy &policy);

template
std::vector<lsst::detection::Footprint::PtrType> lsst::ip::diffim::getCollectionOfFootprintsForPsfMatching(
    lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> const &imageToConvolve,
    lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> const &imageToNotConvolve,
    lsst::pex::policy::Policy &policy);

template
std::vector<lsst::detection::Footprint::PtrType> lsst::ip::diffim::getCollectionOfFootprintsForPsfMatching(
    lsst::afw::image::MaskedImage<double, lsst::afw::image::maskPixelType> const &imageToConvolve,
    lsst::afw::image::MaskedImage<double, lsst::afw::image::maskPixelType> const &imageToNotConvolve,
    lsst::pex::policy::Policy &policy);

template 
lsst::afw::math::KernelList<lsst::afw::math::Kernel> lsst::ip::diffim::computePcaKernelBasis(
    std::vector<lsst::ip::diffim::DiffImContainer<float, lsst::afw::image::maskPixelType> > &diffImContainerList,
    lsst::pex::policy::Policy &policy);

template 
lsst::afw::math::KernelList<lsst::afw::math::Kernel> lsst::ip::diffim::computePcaKernelBasis(
    std::vector<lsst::ip::diffim::DiffImContainer<double, lsst::afw::image::maskPixelType> > &diffImContainerList,
    lsst::pex::policy::Policy &policy);

template 
boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> lsst::ip::diffim::computeSpatiallyVaryingPsfMatchingKernel(
    lsst::afw::math::Function2<double> &kernelFunction,
    lsst::afw::math::Function2<double> &backgroundFunction,
    std::vector<lsst::ip::diffim::DiffImContainer<float, lsst::afw::image::maskPixelType> > &diffImContainerList,
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelBasisList,
    lsst::pex::policy::Policy &policy);
template 
boost::shared_ptr<lsst::afw::math::LinearCombinationKernel> lsst::ip::diffim::computeSpatiallyVaryingPsfMatchingKernel(
    lsst::afw::math::Function2<double> &kernelFunction,
    lsst::afw::math::Function2<double> &backgroundFunction,
    std::vector<lsst::ip::diffim::DiffImContainer<double, lsst::afw::image::maskPixelType> > &diffImContainerList,
    lsst::afw::math::KernelList<lsst::afw::math::Kernel> const &kernelBasisList,
    lsst::pex::policy::Policy &policy);

template
bool lsst::ip::diffim::maskOk(
    lsst::afw::image::Mask<lsst::afw::image::maskPixelType> const &inputMask,
    lsst::afw::image::maskPixelType const badPixelMask);

template
void lsst::ip::diffim::calculateMaskedImageResiduals(
    int &nGoodPixels,
    double &meanOfResiduals,
    double &varianceOfResiduals,
    lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> const &inputImage,
    lsst::afw::image::maskPixelType const badPixelMask);

template
void lsst::ip::diffim::calculateMaskedImageResiduals(
    int &nGoodPixels,
    double &meanOfResiduals,
    double &varianceOfResiduals,
    lsst::afw::image::MaskedImage<double, lsst::afw::image::maskPixelType> const &inputImage,
    lsst::afw::image::maskPixelType const badPixelMask);

template
void lsst::ip::diffim::calculateImageResiduals(
    int &nGoodPixels,
    double &meanOfResiduals,
    double &varianceOfResiduals,
    lsst::afw::image::Image<float> const &inputImage);

template
void lsst::ip::diffim::calculateImageResiduals(
    int &nGoodPixels,
    double &meanOfResiduals,
    double &varianceOfResiduals,
    lsst::afw::image::Image<double> const &inputImage);

template
void lsst::ip::diffim::addFunction(
    lsst::afw::image::Image<float>&,
    lsst::afw::math::Function2<float> const&);

template
void lsst::ip::diffim::addFunction(
    lsst::afw::image::Image<float>&,
    lsst::afw::math::Function2<double> const&);

template
void lsst::ip::diffim::addFunction(
    lsst::afw::image::Image<double>&,
    lsst::afw::math::Function2<float> const&);

template
void lsst::ip::diffim::addFunction(
    lsst::afw::image::Image<double>&,
    lsst::afw::math::Function2<double> const&);
