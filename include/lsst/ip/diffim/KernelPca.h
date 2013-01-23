// -*- lsst-c++ -*-
/**
 * @file KernelPca.h
 *
 * @brief Declaration of KernelPca and KernelPcaVisitor
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_KERNELPCA_H
#define LSST_IP_DIFFIM_KERNELPCA_H

#include "lsst/afw/image.h"
#include "lsst/afw/math.h"

namespace lsst { 
namespace ip { 
namespace diffim { 
namespace detail {

    template <typename ImageT>
    class KernelPca : public lsst::afw::image::ImagePca<ImageT> {
        typedef typename lsst::afw::image::ImagePca<ImageT> Super; ///< Base class
    public:
        typedef typename boost::shared_ptr<KernelPca<ImageT> > Ptr;
        using lsst::afw::image::ImagePca<ImageT>::getEigenImages;
        using lsst::afw::image::ImagePca<ImageT>::getEigenValues;
        using lsst::afw::image::ImagePca<ImageT>::addImage;

        /// Ctor
        explicit KernelPca(bool constantWeight=true) : Super(constantWeight) {}
        
        /// Generate eigenimages that are normalised 
        virtual void analyze();
    };
    
    template<typename PixelT>
    class KernelPcaVisitor : public lsst::afw::math::CandidateVisitor {
    public:
        typedef lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel> ImageT;
        typedef boost::shared_ptr<KernelPcaVisitor<PixelT> > Ptr;
        
        KernelPcaVisitor(boost::shared_ptr<KernelPca<ImageT> > imagePca);
        virtual ~KernelPcaVisitor() {};
        
        lsst::afw::math::KernelList getEigenKernels();
        void processCandidate(lsst::afw::math::SpatialCellCandidate *candidate);
        void subtractMean();
        PTR(ImageT) returnMean() {return _mean;}
    private:
        boost::shared_ptr<KernelPca<ImageT> > _imagePca;  ///< Structure to fill with images
        PTR(ImageT) _mean;                                ///< Mean image calculated before Pca
    };

    template<typename PixelT>
    boost::shared_ptr<KernelPcaVisitor<PixelT> >
    makeKernelPcaVisitor(boost::shared_ptr<KernelPca<typename KernelPcaVisitor<PixelT>::ImageT> > imagePca) {
        return typename KernelPcaVisitor<PixelT>::Ptr(new KernelPcaVisitor<PixelT>(imagePca));
    };

}}}} // end of namespace lsst::ip::diffim::detail

#endif
