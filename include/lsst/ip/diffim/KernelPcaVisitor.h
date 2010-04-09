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
    
    template<typename PixelT>
    class KernelPcaVisitor : public lsst::afw::math::CandidateVisitor {
    public:
        typedef lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel> ImageT;
        typedef boost::shared_ptr<KernelPcaVisitor<PixelT> > Ptr;
        
        KernelPcaVisitor(lsst::afw::image::ImagePca<ImageT> *imagePca);
        virtual ~KernelPcaVisitor() {};
        
        //getEigenKernels();
        void processCandidate(lsst::afw::math::SpatialCellCandidate *candidate);
        void subtractMean();
        ImageT::Ptr returnMean() {return _mean;}
    private:
        lsst::afw::image::ImagePca<ImageT> *_imagePca; ///< Structure to fill with images
        ImageT::Ptr _mean;                             ///< Mean image calculated before Pca
    };

    template<typename PixelT>
    boost::shared_ptr<KernelPcaVisitor<PixelT> >
    makeKernelPcaVisitor(lsst::afw::image::ImagePca<typename KernelPcaVisitor<PixelT>::ImageT> *imagePca) {
        return typename KernelPcaVisitor<PixelT>::Ptr(new KernelPcaVisitor<PixelT>(imagePca));
    }
    
}}}} // end of namespace lsst::ip::diffim::detail

#endif
