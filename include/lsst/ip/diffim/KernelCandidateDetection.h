// -*- lsst-c++ -*-
/**
 * @file KernelCandidateDetection.h
 *
 * @brief Detect candidates for kernels within 2 images 
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_KERNELCANDIDATEDETECTION_H
#define LSST_IP_DIFFIM_KERNELCANDIDATEDETECTION_H

#include "lsst/afw/image/Image.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/pex/policy/Policy.h"

namespace lsst { 
namespace ip { 
namespace diffim {

    /**
     * @brief Search through images for Footprints with no masked pixels
     *
     * @note Runs detection on the template; searches through both images for masked pixels
     *
     * @param templateMaskedImage  MaskedImage that will be convolved with kernel
     * @param scienceMaskedImage   MaskedImage to subtract convolved template from
     * @param policy  Policy for operations; in particular object detection
     *
     * @ingroup ip_diffim
     */    
    template <typename PixelT>
    class KernelCandidateDetection {
    public:
        typedef std::shared_ptr<KernelCandidateDetection> Ptr;
        typedef std::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > MaskedImagePtr;

        KernelCandidateDetection(lsst::pex::policy::Policy const& policy);

        virtual ~KernelCandidateDetection() {};

        void apply(MaskedImagePtr const& templateMaskedImage,
                   MaskedImagePtr const& scienceMaskedImage);
        
        bool growCandidate(lsst::afw::detection::Footprint::Ptr fp,
                           int fpGrowPix,
                           MaskedImagePtr const& templateMaskedImage,
                           MaskedImagePtr const& scienceMaskedImage);

        std::vector<lsst::afw::detection::Footprint::Ptr> getFootprints() {return _footprints;};

    private:
        lsst::pex::policy::Policy _policy;
        lsst::afw::image::MaskPixel _badBitMask;
        std::vector<lsst::afw::detection::Footprint::Ptr> _footprints;
    };


}}} // end of namespace lsst::ip::diffim
        
#endif
