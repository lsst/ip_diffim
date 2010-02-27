// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Implementation of KernelCandidateDetection class
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include "lsst/afw/image/Image.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/pex/exceptions/Exception.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/ip/diffim/ImageSubtract.h"
#include "lsst/ip/diffim/KernelCandidateDetection.h"

namespace afwImage  = lsst::afw::image;
namespace afwDetect = lsst::afw::detection;
namespace pexLog    = lsst::pex::logging; 
namespace pexExcept = lsst::pex::exceptions; 

namespace lsst { 
namespace ip { 
namespace diffim {

    
    template <typename PixelT>
    KernelCandidateDetection<PixelT>::KernelCandidateDetection(
        lsst::pex::policy::Policy const& policy
        ) :
        _policy(policy), 
        _badBitMask(0), 
        _footprints(std::vector<lsst::afw::detection::Footprint::Ptr>()) {

        std::cout << "CAW" << std::endl;

        std::vector<std::string> detBadMaskPlanes = _policy.getStringArray("detBadMaskPlanes");
        for (std::vector<std::string>::iterator mi = detBadMaskPlanes.begin(); 
             mi != detBadMaskPlanes.end(); ++mi){
            try {
                _badBitMask |= afwImage::Mask<afwImage::MaskPixel>::getPlaneBitMask(*mi);
            } catch (pexExcept::Exception& e) {
                pexLog::TTrace<6>("lsst.ip.diffim.KernelCandidateDetection",
                                  "Cannot update bad bit mask with %s", (*mi).c_str());
                pexLog::TTrace<7>("lsst.ip.diffim.KernelCandidateDetection",
                                  e.what());
            }
        }
        pexLog::TTrace<4>("lsst.ip.diffim.KernelCandidateDetection", 
                          "Using bad bit mask %d", _badBitMask);
    }
        
        

    /** 
     * @brief Runs Detection on a single image for significant peaks, and checks
     * returned Footprints for Masked pixels.
     *
     * @note Accepts two MaskedImages, one of which is to be convolved to match the
     * other.  The Detection package is run on either the image to be convolved
     * (assumed to be higher S/N than the other image), or the image to not be
     * convolved (assumed lower S/N; however if you run detection on a very deep
     * template, you might not have significant S/N objects in the science image).
     * The subimages associated with each returned Footprint in both images are
     * checked for Masked pixels; Footprints containing Masked pixels are rejected.
     * The Footprints are grown by an amount specified in the Policy.  The
     * acceptible Footprints are returned in a vector.
     *
     * @return Vector of "clean" Footprints around which Image Subtraction
     * Kernels will be built.
     *
     */
    template <typename PixelT>
    void KernelCandidateDetection<PixelT>::apply(
        MaskedImagePtr const& miToConvolvePtr,
        MaskedImagePtr const& miToNotConvolvePtr
        ) {
        
        // Parse the Policy
        int fpNpixMin                = _policy.getInt("fpNpixMin");
        int fpGrowPix                = _policy.getInt("fpGrowPix");
        
        bool detOnTemplate           = _policy.getBool("detOnTemplate");
        double detThreshold          = _policy.getDouble("detThreshold");
        std::string detThresholdType = _policy.getString("detThresholdType");

        /* reset private variables */
        _footprints.clear();

        // List of Footprints
        std::vector<afwDetect::Footprint::Ptr> footprintListIn;
        
        // Find detections
        afwDetect::Threshold threshold = 
            afwDetect::createThreshold(detThreshold, detThresholdType);
        
        if (detOnTemplate == true) {
            afwDetect::FootprintSet<PixelT> footprintSet(
                *(miToConvolvePtr), 
                threshold,
                "",
                fpNpixMin);
            // Get the associated footprints
            footprintListIn = footprintSet.getFootprints();
            pexLog::TTrace<4>("lsst.ip.diffim.KernelCandidateDetection.apply", 
                              "Found %d total footprints in template above %.3f %s",
                              footprintListIn.size(), detThreshold, detThresholdType.c_str());
        }
        else {
            afwDetect::FootprintSet<PixelT> footprintSet(
                *(miToNotConvolvePtr), 
                threshold,
                "",
                fpNpixMin);
            
            footprintListIn = footprintSet.getFootprints();
            pexLog::TTrace<4>("lsst.ip.diffim.KernelCandidateDetection.apply", 
                              "Found %d total footprints in science image above %.3f %s",
                              footprintListIn.size(), detThreshold, detThresholdType.c_str());
        }    
        
        // Iterate over footprints, look for "good" ones
        for (std::vector<afwDetect::Footprint::Ptr>::iterator i = footprintListIn.begin(); 
             i != footprintListIn.end(); ++i) {
            
            pexLog::TTrace<6>("lsst.ip.diffim.KernelCandidateDetection.apply", 
                              "Processing candidate %d", (*i)->getId());
            growCandidate((*i), fpGrowPix, miToConvolvePtr, miToNotConvolvePtr);
        }
        
        if (_footprints.size() == 0) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              "Unable to find any footprints for Psf matching");
        }
        
        pexLog::TTrace<1>("lsst.ip.diffim.KernelCandidateDetection.apply", 
                          "Found %d clean footprints above threshold %.3f",
                          _footprints.size(), detThreshold);
        
    }
    
    template <typename PixelT>
    bool KernelCandidateDetection<PixelT>::growCandidate(
        afwDetect::Footprint::Ptr fp, 
        int fpGrowPix,
        MaskedImagePtr const& miToConvolvePtr,
        MaskedImagePtr const& miToNotConvolvePtr
        ) {

        /* Grow the Footprint by the requested (optimal) amount; if you happen
         * to find a masked pixel within the grown Footprint, try again growing
         * by the minimal acceptable amount.  We don't want to be too sensitive
         * to masked pixels in the wings of the stamp.  If it fails this minimal
         * grow just don't use it, its too close to a masked pixel 
         */

        int fpGrowMin = _policy.getInt("fpGrowMin");
        int fpNpixMax = _policy.getInt("fpNpixMax");

        /* Functor to search through the images for masked pixels within *
         * candidate footprints.  Might want to consider changing the default
         * mask planes it looks through.
         */
        FindSetBits<afwImage::Mask<afwImage::MaskPixel> > fsb;
        
        afwImage::BBox fpBBox = fp->getBBox();
        /* Failure Condition 1) 
         * 
         * Footprint has too many pixels off the bat.  We don't want to throw
         * away these guys, they have alot of signal!  Lets just use the core of
         * it.
         * 
         */
        if (fp->getNpix() > fpNpixMax) {
            pexLog::TTrace<6>("lsst.ip.diffim.KernelCandidateDetection.apply", 
                              "Footprint has too many pix: %d (max =%d)", 
                              fp->getNpix(), fpNpixMax);
            
            int xc = int(0.5 * (fpBBox.getX0() + fpBBox.getX1()));
            int yc = int(0.5 * (fpBBox.getY0() + fpBBox.getY1()));
            afwDetect::Footprint::Ptr fpCore(
                new afwDetect::Footprint(afwImage::BBox(afwImage::PointI(xc, yc)))
                );
            return growCandidate(fpCore, fpGrowPix, miToConvolvePtr, miToNotConvolvePtr);
        } 
        
        /* Grow the footprint
         * flag true  = isotropic grow   = slow
         * flag false = 'manhattan grow' = fast
         * 
         * The manhattan masks are rotated 45 degree w.r.t. the coordinate
         * system.  They intersect the vertices of the rectangle that would
         * connect pixels (X0,Y0) (X1,Y0), (X0,Y1), (X1,Y1).
         * 
         * The isotropic masks do take considerably longer to grow and are
         * basically elliptical.  X0, X1, Y0, Y1 delimit the extent of the
         * ellipse.
         * 
         * In both cases, since the masks aren't rectangles oriented with
         * the image coordinate system, when we DO extract such rectangles
         * as subimages for kernel fitting, some corner pixels can be found
         * in multiple subimages.
         * 
         */
        afwDetect::Footprint::Ptr fpGrow = 
            afwDetect::growFootprint(fp, fpGrowPix, false);
        
        /* Next we look at the image within this Footprint.  To do this we
         * create a subimage using a BBox.  When this happens, the BBox
         * addresses the image in its own coordinate system, not the parent
         * coordinate system.  Therefore you need to shift the BBox by its XY0.
         */
        afwImage::BBox fpGrowBBox = fpGrow->getBBox();
        pexLog::TTrace<8>("lsst.ip.diffim.KernelCandidateDetection.apply", 
                          "Footprint in parent : %d,%d -> %d,%d -> %d,%d",
                          fpGrowBBox.getX0(), fpGrowBBox.getY0(), 
                          int(0.5 * (fpGrowBBox.getX0() + fpGrowBBox.getX1())),
                          int(0.5 * (fpGrowBBox.getY0() + fpGrowBBox.getY1())),
                          fpGrowBBox.getX1(), fpGrowBBox.getY1());

        fpGrowBBox.shift(-miToConvolvePtr->getX0(), -miToConvolvePtr->getY0());
        pexLog::TTrace<8>("lsst.ip.diffim.KernelCandidateDetection.apply", 
                          "Footprint in image : %d,%d -> %d,%d -> %d,%d",
                          fpGrowBBox.getX0(), fpGrowBBox.getY0(), 
                          int(0.5 * (fpGrowBBox.getX0() + fpGrowBBox.getX1())),
                          int(0.5 * (fpGrowBBox.getY0() + fpGrowBBox.getY1())),
                          fpGrowBBox.getX1(), fpGrowBBox.getY1());

        /* Failure Condition 2) 
         * Grown off the image
         */
        bool belowOriginX = fpGrowBBox.getX0() < 0;
        bool belowOriginY = fpGrowBBox.getY0() < 0;
        bool offImageX    = fpGrowBBox.getX1() > (miToConvolvePtr->getWidth() - 1);
        bool offImageY    = fpGrowBBox.getY1() > (miToConvolvePtr->getHeight() - 1);
        if (belowOriginX || belowOriginY || offImageX || offImageY) {
            pexLog::TTrace<6>("lsst.ip.diffim.KernelCandidateDetection.apply", 
                              "Footprint grown off image"); 

            if (fpGrowPix == fpGrowMin) {
                return false;
            }
            else {
                return growCandidate(fp, fpGrowMin, miToConvolvePtr, miToNotConvolvePtr);
            }
        }
        
        /* Grab subimages; report any exception */
        bool subimageHasFailed = false;
        try {
            afwImage::MaskedImage<PixelT> subImageToConvolve(*(miToConvolvePtr), fpGrowBBox);
            afwImage::MaskedImage<PixelT> subImageToNotConvolve(*(miToNotConvolvePtr), fpGrowBBox);
            
            // Search for any masked pixels within the footprint
            fsb.apply(*(subImageToConvolve.getMask()));
            if (fsb.getBits() & _badBitMask) {
                pexLog::TTrace<6>("lsst.ip.diffim.KernelCandidateDetection.apply", 
                                  "Footprint has masked pix (vals=%d) in image to convolve", 
                                  fsb.getBits()); 
                subimageHasFailed = true;
            }
            
            fsb.apply(*(subImageToNotConvolve.getMask()));
            if (fsb.getBits() & _badBitMask) {
                pexLog::TTrace<6>("lsst.ip.diffim.KernelCandidateDetection.apply", 
                                  "Footprint has masked pix (vals=%d) in image not to convolve", 
                                  fsb.getBits());
                subimageHasFailed = true;
            }
            
        } catch (pexExcept::Exception& e) {
            pexLog::TTrace<6>("lsst.ip.diffim.KernelCandidateDetection.apply",
                              "Exception caught extracting Footprint");
            pexLog::TTrace<7>("lsst.ip.diffim.KernelCandidateDetection.apply",
                              e.what());
            subimageHasFailed = true;
        }
        if (subimageHasFailed) {
            if (fpGrowPix == fpGrowMin) {
                return false;
            }
            else {
                return growCandidate(fp, fpGrowMin, miToConvolvePtr, miToNotConvolvePtr);
            }
        } else {
            /* We have a good candidate */
            _footprints.push_back(fpGrow);
            return true;
        }
    }

/***********************************************************************************************************/
//
// Explicit instantiations
//
    typedef float PixelT;
    template class KernelCandidateDetection<PixelT>;

}}} // end of namespace lsst::ip::diffim

