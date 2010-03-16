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
#include "lsst/ip/diffim/KernelSolution.h"

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
    template <typename _PixelT, typename VarT=lsst::afw::image::VariancePixel>
    class KernelCandidate : public lsst::afw::math::SpatialCellImageCandidate<
        lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel> 
        > {
    public: 
        typedef lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel> ImageT;
        typedef _PixelT PixelT;         // _after_ lsst::afw::math::Kernel::Pixel

    private:
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::_image;

    public:
        typedef boost::shared_ptr<KernelCandidate> Ptr;
        typedef boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > MaskedImagePtr;
        typedef boost::shared_ptr<lsst::afw::image::Image<VarT> > VariancePtr;

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
        KernelCandidate(float const xCenter,
                        float const yCenter, 
                        MaskedImagePtr const& miToConvolvePtr,
                        MaskedImagePtr const& miToNotConvolvePtr,
                        VariancePtr varianceEstimate,
                        lsst::pex::policy::Policy const& policy);        
        /// Destructor
        virtual ~KernelCandidate() {};

        /**
         * @brief Return Candidate rating
         * 
         * @note Required method for use by SpatialCell; e.g. total flux
         */
        double getCandidateRating() const {return _coreFlux;}

        /**
         * @brief Return pointers to the image pixels used in kernel determination
         */
        MaskedImagePtr getMiToConvolvePtr() {return _miToConvolvePtr;}
        MaskedImagePtr getMiToNotConvolvePtr() {return _miToNotConvolvePtr;}

        /**
         * @brief Return results of kernel solution
         * 
         */
        lsst::afw::math::Kernel::Ptr getOrigKernel() const {return _kernelSolutionOrig->getKernel();}
        double getOrigBackground() const {return _kernelSolutionOrig->getBackground();}
        double getOrigKsum() const {return _kernelSolutionOrig->getKsum();}
        bool isInitialized() const {return _isInitialized;}

        lsst::afw::math::Kernel::Ptr getPcaKernel() const {return _kernelSolutionPca->getKernel();}
        double getPcaBackground() const {return _kernelSolutionPca->getBackground();}
        double getPcaKsum() const {return _kernelSolutionPca->getKsum();}

        void setVariance(VariancePtr var) {_varianceEstimate = var;}

        /**
         * @brief Return pointers to the image of the kernel.  Needed for Pca.
         */
        typename ImageT::ConstPtr getOrigImage() const;
        typename ImageT::Ptr copyOrigImage() const;
        typename ImageT::ConstPtr getPcaImage() const;
        typename ImageT::Ptr copyPcaImage() const;

        /** 
         * @brief Calculate associated difference image using internal solutions
         */
        lsst::afw::image::MaskedImage<PixelT> returnOrigDifferenceImage();
        lsst::afw::image::MaskedImage<PixelT> returnPcaDifferenceImage();

        /** 
         * @brief Calculate associated difference image using input kernel and background.
         * 
         * @note Useful for spatial modeling
         */
        lsst::afw::image::MaskedImage<PixelT> returnDifferenceImage(
            lsst::afw::math::Kernel::Ptr kernel,
            double background
            );

        /** 
         * @brief Core functionality of KernelCandidate, to build and fill a KernelSolution
         *
         * @note This is an expensive step involving matrix math, and one that
         * may be called multiple times per candidate.  Use cases are:
         *
         *  o _isInitialized = false.  This is a constructed but not initialized
         *  KernelCandidate.  When build() is called, M and B are derived from
         *  the MaskedImages and the basisList.  KernelCandidate owns the
         *  knowledge of how to fill this KernelSolution; the solution knows how
         *  to solve itself and how to turn that into an output kernel.  This
         *  solution ends up being _kernelSolution0.
         *
         *  o _isInitialized = true.  This is for when build() is re-called
         *  using a different basis list, e.g. one based on Pca.  We need to use
         *  M and B for the spatial modeling, but do *not* want to override the
         *  original KernelSolution.  This solution ends up as
         *  _kernelSolutionCurrent.
         */
        //void build(boost::shared_ptr<lsst::afw::math::KernelList> const& basisList);

        /** 
         * @brief Build KernelSolution matrices for M x = B with regularization matrix H 
         *
         * @note Modified equation is (Mt.M + lambda H) x = Mt.B with lambda a
         * degree of freedom describing the "strength" of the regularization.
         * The larger the value of lambda, the smoother the kernel, but the
         * larger the residuals in the difference image.
         *
         * @note A value of lambda = Trace(Mt.M) / Tr(H) will yield essentially
         * equivalent power in the kernel smoothness and in the diffim quality.
         * We scale this estimate by lambdaScaling to give more/less
         * consideration to the smoothness of the kernel.
         */
        void build(
            boost::shared_ptr<lsst::afw::math::KernelList> const& basisList,
            boost::shared_ptr<Eigen::MatrixXd> hMat = boost::shared_ptr<Eigen::MatrixXd>()
            );


    private:
        MaskedImagePtr _miToConvolvePtr;                    ///< Subimage around which you build kernel
        MaskedImagePtr _miToNotConvolvePtr;                 ///< Subimage around which you build kernel
        VariancePtr _varianceEstimate;                      ///< Estimate of the local variance
        lsst::pex::policy::Policy _policy;                  ///< Policy
        double _coreFlux;                                   ///< Mean S/N in the science image
        bool _isInitialized;                                ///< Has the kernel been built

        /* best single raw kernel */
        boost::shared_ptr<StaticKernelSolution> _kernelSolutionOrig;    ///< Original basis kernel solution

        /* with Pca basis */
        boost::shared_ptr<StaticKernelSolution> _kernelSolutionPca;     ///< Most recent kernel solution

    };


    /**
     * @brief Return a KernelCandidate pointer of the right sort
     *
     * @param xCenter  X-center of candidate
     * @param yCenter  Y-center of candidate
     * @param miToConvolvePtr  Template subimage 
     * @param miToNotConvolvePtr  Science image subimage
     * @param policy   Policy file for creation of rating
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
