// -*- lsst-c++ -*-
/**
 * @file KernelPcaVisitor.h
 *
 * @brief Declaration of KernelPcaVisitor 
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_KERNELPCAVISITOR_H
#define LSST_IP_DIFFIM_KERNELPCAVISITOR_H

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
        
        KernelPcaVisitor(KernelPca<ImageT> *imagePca);
        virtual ~KernelPcaVisitor() {};
        
        lsst::afw::math::KernelList getEigenKernels();
        void processCandidate(lsst::afw::math::SpatialCellCandidate *candidate);
        void subtractMean();
        PTR(ImageT) returnMean() {return _mean;}
    private:
        KernelPca<ImageT> *_imagePca;  ///< Structure to fill with images
        ImageT::Ptr _mean;             ///< Mean image calculated before Pca
    };

    template<typename PixelT>
    boost::shared_ptr<KernelPcaVisitor<PixelT> >
    makeKernelPcaVisitor(KernelPca<typename KernelPcaVisitor<PixelT>::ImageT> *imagePca) {
        return typename KernelPcaVisitor<PixelT>::Ptr(new KernelPcaVisitor<PixelT>(imagePca));
    };

}}}} // end of namespace lsst::ip::diffim::detail

#endif
