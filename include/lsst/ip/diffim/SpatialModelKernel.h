// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Class used by SpatialModelCell for spatial Kernel fitting
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup afw
 */

#ifndef LSST_IP_DIFFIM_SPATIALMODELKERNEL_H
#define LSST_IP_DIFFIM_SPATIALMODELKERNEL_H

#include <boost/shared_ptr.hpp>
#include <Eigen/Core>

#include <lsst/afw/math/SpatialCell.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/math/KernelFunctions.h>
#include <lsst/pex/policy/Policy.h>
#include <lsst/afw/detection/Footprint.h>
#include <lsst/sdqa/SdqaRating.h>

#include <lsst/ip/diffim/ImageSubtract.h>

namespace lsst {
namespace ip {
namespace diffim {

    /** 
     * 
     * @brief Class stored in SpatialCells for spatial Kernel fitting
     * 
     * KernelCandidate is a single Kernel derived around a source.  We'll assign
     * them to sets of SpatialCells; these sets will then be used to fit a
     * spatial model to the PSF.
     */    

    template <typename PixelT>
    class KernelCandidate : public lsst::afw::math::SpatialCellImageCandidate<lsst::afw::image::Image<lsst::afw::math::Kernel::PixelT> > {
    public: 
        typedef lsst::afw::image::Image<lsst::afw::math::Kernel::PixelT> ImageT;

    private:
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::_image;

    public:
        typedef boost::shared_ptr<KernelCandidate> Ptr;
        typedef boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > MaskedImagePtr;

        /** Constructor
         *
         * @param xCenter Col position of object
         * @param yCenter Row position of object
         * @param miToConvolvePtr  Pointer to template image
         * @param miToNotConvolvePtr  Pointer to science image
         */
        KernelCandidate(float const xCenter,
                        float const yCenter, 
                        MaskedImagePtr const& miToConvolvePtr,
                        MaskedImagePtr const& miToNotConvolvePtr) :
            lsst::afw::math::SpatialCellImageCandidate<ImageT>(xCenter, yCenter),
            _miToConvolvePtr(miToConvolvePtr),
            _miToNotConvolvePtr(miToNotConvolvePtr),
            _haveImage(false),
            _haveKernel(false) {
        }
        
        /// Destructor
        ~KernelCandidate() {};

        /**
         * Return Cell rating
         * 
         * @note Required method for use by SpatialCell; e.g. total flux
         */
        double getCandidateRating() const { return 100; }

        MaskedImagePtr getMiToConvolvePtr() {return _miToConvolvePtr;}
        MaskedImagePtr getMiToNotConvolvePtr() {return _miToNotConvolvePtr;}

        typename ImageT::ConstPtr getImage() const;
        typename ImageT::Ptr copyImage() const;
        lsst::afw::math::Kernel::PtrT getKernel() const;
        Eigen::MatrixXd getM()  {return _M;}
        Eigen::VectorXd getB()  {return _B;}
        
        void setKernel(lsst::afw::math::Kernel::PtrT kernel) {_kernel = kernel; _haveKernel = true;}
        void setBackground(double background) {_background = background;}
        void setM(Eigen::MatrixXd M) {_M = M;}
        void setB(Eigen::VectorXd B) {_B = B;}
        
    private:
        MaskedImagePtr _miToConvolvePtr;                    ///< Subimage around which you build kernel
        MaskedImagePtr _miToNotConvolvePtr;                 ///< Subimage around which you build kernel

        lsst::afw::math::Kernel::PtrT _kernel;              ///< Derived single-object convolution kernel
        double _background;                                 ///< Derived differential background estimate
        Eigen::MatrixXd _M;                                 ///< Derived least squares matrix
        Eigen::VectorXd _B;                                 ///< Derived least squares vector

        bool mutable _haveImage;                            ///< do we have an Image to return?
        bool mutable _haveKernel;                           ///< do we have a Kernel to return?
    };

    /**
     * Return a KernelCandidate pointer of the right sort
     *
     */
    template <typename PixelT>
    typename KernelCandidate<PixelT>::Ptr
    makeKernelCandidate(float const xCenter,
                        float const yCenter, 
                        boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > const& miToConvolvePtr,
                        boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > const& miToNotConvolvePtr) {
        
        return typename KernelCandidate<PixelT>::Ptr(new KernelCandidate<PixelT>(xCenter, yCenter,
                                                                                 miToConvolvePtr,
                                                                                 miToNotConvolvePtr));
    }

    template<typename PixelT>
    std::pair<lsst::afw::math::LinearCombinationKernel::PtrT, std::vector<double> >
    createPcaBasisFromCandidates(lsst::afw::math::SpatialCellSet const& psfCells,
                                 lsst::pex::policy::Policy const& policy);

    template<typename PixelT>
    std::pair<bool, double>
    fitSpatialKernelFromCandidates(lsst::afw::math::Kernel *kernel,
                                   lsst::afw::math::SpatialCellSet const& psfCells,
                                   lsst::pex::policy::Policy const& policy);
    
    
}}}

#endif // LSST_IP_DIFFIM_SPATIALMODELKERNEL_H

