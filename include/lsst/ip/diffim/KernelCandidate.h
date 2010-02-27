// -*- lsst-c++ -*-
/**
 * @file KernelCandidate.h
 *
 * @brief Class used by SpatialModelCell for spatial Kernel fitting
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_KERNELCANDIDATE_H
#define LSST_IP_DIFFIM_KERNELCANDIDATE_H

#include "boost/shared_ptr.hpp"
#include "Eigen/Core"

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"

namespace lsst { 
namespace ip { 
namespace diffim {

    /** 
     * @brief Class stored in SpatialCells for spatial Kernel fitting
     * 
     * @note KernelCandidate is a single Kernel derived around a source.  We'll assign
     * them to sets of SpatialCells; these sets will then be used to fit a
     * spatial model to the Kernel.
     *
     * @ingroup ip_diffim
     */    
    template <typename _PixelT>
    class KernelCandidate : public lsst::afw::math::SpatialCellImageCandidate<
        lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel> 
        > {
    public: 
        typedef lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel> ImageT;
        typedef _PixelT PixelT;         // _after_ using lsst::afw::math::Kernel::Pixel

    private:
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::_image;

    public:
        typedef boost::shared_ptr<KernelCandidate> Ptr;
        typedef boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > MaskedImagePtr;

        /**
	 * @brief Constructor
         *
         * @param xCenter Col position of object
         * @param yCenter Row position of object
         * @param miToConvolvePtr  Pointer to template image
         * @param miToNotConvolvePtr  Pointer to science image
         * @param policy  Policy file
         */
        KernelCandidate(float const xCenter,
                        float const yCenter, 
                        MaskedImagePtr const& miToConvolvePtr,
                        MaskedImagePtr const& miToNotConvolvePtr,
                        lsst::pex::policy::Policy const& policy);        
        /// Destructor
        virtual ~KernelCandidate() {};

        /**
         * Return Cell rating
         * 
         * @note Required method for use by SpatialCell; e.g. total flux
         */
        double getCandidateRating() const { return _coreFlux; }

        MaskedImagePtr getMiToConvolvePtr() {return _miToConvolvePtr;}
        MaskedImagePtr getMiToNotConvolvePtr() {return _miToNotConvolvePtr;}

        /** 
         * @brief Calculate associated difference image using member _kernel/_background
         */
        lsst::afw::image::MaskedImage<PixelT> returnDifferenceImage();
	
        /** 
         * @brief Calculate associated difference image using provided kernel and background
         */
        lsst::afw::image::MaskedImage<PixelT> returnDifferenceImage(
            lsst::afw::math::Kernel::Ptr kernel,
            double background
            );


        typename ImageT::ConstPtr getImage() const;
        typename ImageT::Ptr copyImage() const;
        double getKsum() const;
        lsst::afw::math::Kernel::Ptr getKernel() const;
        double getBackground() const;
        boost::shared_ptr<Eigen::MatrixXd> const getM()  {return _mMat;}
        boost::shared_ptr<Eigen::VectorXd> const getB()  {return _bVec;}
        bool hasKernel() {return _haveKernel;}
        
        void setKernel(lsst::afw::math::Kernel::Ptr kernel);
        void setBackground(double background) {_background = background;}

        void setM(boost::shared_ptr<Eigen::MatrixXd> mMat) {_mMat = mMat;}
        void setB(boost::shared_ptr<Eigen::VectorXd> bVec) {_bVec = bVec;}

    private:
        MaskedImagePtr _miToConvolvePtr;                    ///< Subimage around which you build kernel
        MaskedImagePtr _miToNotConvolvePtr;                 ///< Subimage around which you build kernel
        lsst::pex::policy::Policy _policy;                  ///< Policy
        double _coreFlux;                                   ///< Mean S/N in the science image

        lsst::afw::math::Kernel::Ptr _kernel;               ///< Derived single-object convolution kernel
        double _kSum;                                       ///< Derived kernel sum
        double _background;                                 ///< Derived differential background estimate

        boost::shared_ptr<Eigen::MatrixXd> _mMat;           ///< Derived least squares matrix
        boost::shared_ptr<Eigen::VectorXd> _bVec;           ///< Derived least squares vector

        bool mutable _haveKernel;                           ///< do we have a Kernel to return?
    };

    /**
     * @brief Return a KernelCandidate pointer of the right sort
     *
     * @param xCenter  X-center of candidate
     * @param yCenter  Y-center of candidate
     * @param miToConvolvePtr  Template subimage 
     * @param miToNotConvolvePtr  Science image subimage
     *
     * @ingroup ip_diffim
     */
    template <typename PixelT>
    typename KernelCandidate<PixelT>::Ptr
    makeKernelCandidate(float const xCenter,
                        float const yCenter, 
                        boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > const& miToConvolvePtr,
                        boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > const& miToNotConvolvePtr,
                        lsst::pex::policy::Policy const& policy){
        
        return typename KernelCandidate<PixelT>::Ptr(new KernelCandidate<PixelT>(xCenter, yCenter,
                                                                                 miToConvolvePtr,
                                                                                 miToNotConvolvePtr,
                                                                                 policy));
    }


}}} // end of namespace lsst::ip::diffim

#endif
