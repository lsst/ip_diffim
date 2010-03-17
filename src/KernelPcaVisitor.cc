// -*- lsst-c++ -*-
/**
 * @file KernelPcaVisitor.h
 *
 * @brief Implementation of KernelPcaVisitor 
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/ip/diffim/KernelCandidate.h"
#include "lsst/ip/diffim/KernelPcaVisitor.h"

namespace afwMath        = lsst::afw::math;
namespace afwImage       = lsst::afw::image;
namespace pexLogging     = lsst::pex::logging; 
namespace pexExcept      = lsst::pex::exceptions; 

namespace lsst { 
namespace ip { 
namespace diffim {
namespace detail {

    /**
     * @class KernelPcaVisitor
     * @ingroup ip_diffim
     *
     * @brief A class to run a PCA on all candidate kernels (represented as Images)
     *
     * @code
        afwImage::ImagePca<ImageT> imagePca;
        detail::KernelPcaVisitor<PixelT> importStarVisitor(&imagePca);
        kernelCells.visitCandidates(&importStarVisitor, nStarPerCell);
        importStarVisitor.subtractMean();
        imagePca.analyze();
        std::vector<typename ImageT::Ptr> eigenImages = imagePca.getEigenImages();
        afwMath::KernelList kernelListRaw;
        kernelListRaw.push_back(afwMath::Kernel::Ptr(
                                    new afwMath::FixedKernel(
                                        afwImage::Image<afwMath::Kernel::Pixel>
                                        (*(importStarVisitor.returnMean()), true))));
        int const ncomp = static_cast<int>(eigenImages.size()) - 1; // -1 since we have subtracted mean
        for (int j = 0; j != ncomp; ++j) {
            kernelListRaw.push_back(afwMath::Kernel::Ptr(
                                        new afwMath::FixedKernel(
                                            afwImage::Image<afwMath::Kernel::Pixel>
                                            (*eigenImages[j], true))));
        }
     * @endcode
     *
     * @note Works in concert with a afwMath::SpatialCellSet and afwImage::ImagePca
     * to create a Karhunen-Loeve basis from all the good KernelCandidates.  This
     * class adds the extra functionality to subtract off the mean kernel from all
     * entries in afwImage::ImagePca, which makes the resulting basis more compact.
     * The user needs to manually add this mean image into the resulting basis list
     * after imagePca.analyze() is called.
     *
     * @note afwImage::ImagePca weights objects of different brightness differently.
     * However we don't necessarily want images with larger kernel sums to have more
     * weight.  Each kernel should have constant weight in the Pca.  For simplicity
     * we scale them to have the same kernel sum, 1.0, and send to ImagePca that the
     * flux (weight) is 1.0.
     * 
     */
    template<typename PixelT>
    KernelPcaVisitor<PixelT>::KernelPcaVisitor(
        lsst::afw::image::ImagePca<ImageT> *imagePca ///< Set of Images to initialise
        ) :
        afwMath::CandidateVisitor(),
        _imagePca(imagePca),
        _mean() 
    {};

    template<typename PixelT>
    void KernelPcaVisitor<PixelT>::processCandidate(lsst::afw::math::SpatialCellCandidate *candidate) {
        
        KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
        if (kCandidate == NULL) {
            throw LSST_EXCEPT(pexExcept::LogicErrorException,
                              "Failed to cast SpatialCellCandidate to KernelCandidate");
        }
        pexLogging::TTrace<6>("lsst.ip.diffim.SetPcaImageVisitor.processCandidate", 
                              "Processing candidate %d", kCandidate->getId());
        
        try {
            /* Normalize to unit sum */
            ImageT::Ptr kImage = kCandidate->getKernelSolution(
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
        pexLogging::TTrace<6>("lsst.ip.diffim.SetPcaImageVisitor.subtractMean", 
                              "Subtracting mean feature before Pca");
        
        _mean = _imagePca->getMean();
        afwImage::ImagePca<ImageT>::ImageList imageList = _imagePca->getImageList();
        for (typename afwImage::ImagePca<ImageT>::ImageList::const_iterator ptr = imageList.begin(), 
                 end = imageList.end(); ptr != end; ++ptr) {
            **ptr -= *_mean;
        }
    }

    typedef float PixelT;
    template class KernelPcaVisitor<PixelT>;

}}}} // end of namespace lsst::ip::diffim::detail
