// -*- lsst-c++ -*-
/**
 * @file KernelPca.cc
 *
 * @brief Implementation of KernelPca and KernelPcaVisitor 
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/log/Log.h"
#include "lsst/pex/exceptions/Runtime.h"

#include "lsst/ip/diffim/KernelCandidate.h"
#include "lsst/ip/diffim/KernelPca.h"

namespace afwMath        = lsst::afw::math;
namespace afwImage       = lsst::afw::image;
namespace pexExcept      = lsst::pex::exceptions; 

namespace lsst { 
namespace ip { 
namespace diffim {
namespace detail {

    /**
     * @class KernelPcaVisitor
     *
     * @ingroup ip_diffim
     *
     * @brief A class to run a PCA on all candidate kernels (represented as
     * Images).
     *
     * @note Templated on the pixel types of the MaskedImages it will be
     * visiting (typically float).
     *
     * @note Works in concert with a afwMath::SpatialCellSet and ip::Diffim
     * KernelPca to create a Karhunen-Loeve basis from all the good
     * KernelCandidates.  This class adds the extra functionality to subtract
     * off the mean kernel from all entries, which makes the resulting basis
     * more compact.  The user needs to manually add this mean image into the
     * resulting basis list after imagePca.analyze() is called.
     * 
     * @note KernelPca (and base class afwImage::ImagePca) weight objects of
     * different brightness differently.  However we don't necessarily want
     * images with larger kernel sums to have more weight.  Each kernel should
     * have constant weight in the Pca.  For simplicity we scale them to have
     * the same kernel sum, 1.0, and send to ImagePca that the flux (weight) is
     * 1.0.
     * 
     */
    template<typename PixelT>
    KernelPcaVisitor<PixelT>::KernelPcaVisitor(
        std::shared_ptr<KernelPca<ImageT> > imagePca ///< Set of Images to initialise
        ) :
        afwMath::CandidateVisitor(),
        _imagePca(imagePca),
        _mean() 
    {};

    template<typename PixelT>
    lsst::afw::math::KernelList KernelPcaVisitor<PixelT>::getEigenKernels() {
        afwMath::KernelList kernelList;

        std::vector<typename ImageT::Ptr> eigenImages = _imagePca->getEigenImages();
        int ncomp = eigenImages.size();

        if (_mean) {
            kernelList.push_back(afwMath::Kernel::Ptr(
                                     new afwMath::FixedKernel(
                                         afwImage::Image<afwMath::Kernel::Pixel>
                                         (*_mean, true))));
        }
        for (int i = 0; i < ncomp; i++) {
            afwImage::Image<afwMath::Kernel::Pixel> img = 
                afwImage::Image<afwMath::Kernel::Pixel>(*eigenImages[i], true);
            kernelList.push_back(afwMath::Kernel::Ptr(
                                     new afwMath::FixedKernel(img)
                                     ));
        }

        return kernelList;
    }

    template<typename PixelT>
    void KernelPcaVisitor<PixelT>::processCandidate(lsst::afw::math::SpatialCellCandidate *candidate) {
        
        KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
        if (kCandidate == NULL) {
            throw LSST_EXCEPT(pexExcept::LogicError,
                              "Failed to cast SpatialCellCandidate to KernelCandidate");
        }
        LOGL_DEBUG("TRACE5.ip.diffim.SetPcaImageVisitor.processCandidate",
                   "Processing candidate %d", kCandidate->getId());
        
        try {
            /* Normalize to unit sum */
            PTR(ImageT) kImage = kCandidate->getKernelSolution(
                KernelCandidate<PixelT>::ORIG)->makeKernelImage();
            *kImage           /= kCandidate->getKernelSolution(
                KernelCandidate<PixelT>::ORIG)->getKsum();
            /* Tell imagePca they have the same weighting in the Pca */
            _imagePca->addImage(kImage, 1.0);
        } catch(pexExcept::Exception &e) {
            return;
        }
    }

    template<typename PixelT>
    void KernelPcaVisitor<PixelT>::subtractMean() {
        /* 
           If we don't subtract off the mean before we do the Pca, the
           subsequent terms carry less of the power than if you do subtract
           off the mean.  Explicit example:
           
           With mean subtraction:
             DEBUG: Eigenvalue 0 : 0.010953 (0.373870 %)
             DEBUG: Eigenvalue 1 : 0.007927 (0.270604 %)
             DEBUG: Eigenvalue 2 : 0.001393 (0.047542 %)
             DEBUG: Eigenvalue 3 : 0.001092 (0.037261 %)
             DEBUG: Eigenvalue 4 : 0.000829 (0.028283 %)
           
           Without mean subtraction:
             DEBUG: Eigenvalue 0 : 0.168627 (0.876046 %)
             DEBUG: Eigenvalue 1 : 0.007935 (0.041223 %)
             DEBUG: Eigenvalue 2 : 0.006049 (0.031424 %)
             DEBUG: Eigenvalue 3 : 0.001188 (0.006173 %)
             DEBUG: Eigenvalue 4 : 0.001050 (0.005452 %)

           After the first term above, which basically represents the mean,
           the remaining terms carry less of the power than if you do
           subtract off the mean.  (0.041223/(1-0.876046) < 0.373870).
         */
        LOGL_DEBUG("TRACE5.ip.diffim.KernelPcaVisitor.subtractMean",
                   "Subtracting mean feature before Pca");
        
        _mean = _imagePca->getMean();
        KernelPca<ImageT>::ImageList imageList = _imagePca->getImageList();
        for (typename KernelPca<ImageT>::ImageList::const_iterator ptr = imageList.begin(), 
                 end = imageList.end(); ptr != end; ++ptr) {
            **ptr -= *_mean;
        }
    }

    /**
     * @class KernelPca
     *
     * @ingroup ip_diffim
     *
     * @brief Overrides the analyze method of base class afwImage::ImagePca
     *
     * @note Templated on the Image types it is running on (typically
     * [exclusively?] afwMath::Kernel::Pixel, which is double)
     *
     * @note This override normalizes the resulting eigenImages to have peak
     * value of 1.0.
     *
     */
    template <typename ImageT>
    void KernelPca<ImageT>::analyze()
    {
        Super::analyze();
        
        typename Super::ImageList const &eImageList = this->getEigenImages();
        typename Super::ImageList::const_iterator iter = eImageList.begin(), end = eImageList.end();
        for (size_t i = 0; iter != end; ++i, ++iter) {
            PTR(ImageT) eImage = *iter;
            
            /*
             * Normalise eigenImages to have a maximum of 1.0.  For n > 0 they
             * (should) have mean == 0, so we can't use that to normalize
             */
            afwMath::Statistics stats = afwMath::makeStatistics(*eImage, (afwMath::MIN | afwMath::MAX));
            double const min = stats.getValue(afwMath::MIN);
            double const max = stats.getValue(afwMath::MAX);
            
            double const extreme = (fabs(min) > max) ? min :max;
            if (extreme != 0.0) {
                *eImage /= extreme;
            }
        }
    }


    typedef float PixelT;
    template class KernelPcaVisitor<PixelT>;
    template class KernelPca<afwImage::Image<afwMath::Kernel::Pixel> >;

}}}} // end of namespace lsst::ip::diffim::detail
