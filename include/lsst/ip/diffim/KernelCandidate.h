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
    template <typename _PixelT>
    class KernelCandidate :
        public lsst::afw::math::SpatialCellImageCandidate<lsst::afw::math::Kernel::Pixel> {
    public: 
        typedef lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel> ImageT;
        typedef _PixelT PixelT;         // _after_ lsst::afw::math::Kernel::Pixel

    private:
        using lsst::afw::math::SpatialCellImageCandidate<lsst::afw::math::Kernel::Pixel>::_image;

    public:
        typedef boost::shared_ptr<KernelCandidate> Ptr;
        typedef boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > MaskedImagePtr;
        typedef boost::shared_ptr<lsst::afw::image::Image<lsst::afw::image::VariancePixel> > VariancePtr;

        enum CandidateSwitch {
            ORIG    = 0,
            PCA     = 1,
            RECENT  = 2
        };

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
        lsst::afw::math::Kernel::Ptr getKernel(CandidateSwitch cand) const;
        double getBackground(CandidateSwitch cand) const;
        double getKsum(CandidateSwitch cand) const;
        PTR(ImageT) getKernelImage(CandidateSwitch cand) const;
        CONST_PTR(ImageT) getImage() const; // For SpatialCellImageCandidate
        boost::shared_ptr<StaticKernelSolution<PixelT> > getKernelSolution(CandidateSwitch cand) const; 
        
        /** 
         * @brief Calculate associated difference image using internal solutions
         */
        lsst::afw::image::MaskedImage<PixelT> getDifferenceImage(CandidateSwitch cand);

        /** 
         * @brief Calculate associated difference image using input kernel and background.
         * 
         * @note Useful for spatial modeling
         */
        lsst::afw::image::MaskedImage<PixelT> getDifferenceImage(
            lsst::afw::math::Kernel::Ptr kernel,
            double background
            );
        
        
        bool isInitialized() const {return _isInitialized;}


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

        /*
         * @note This method uses an estimate of the variance which is the
         * straight difference of the 2 images.  If requested in the Policy
         * ("iterateSingleKernel"), the kernel will be rebuilt using the
         * variance of the difference image resulting from this first
         * approximate step.  This is particularly useful when convolving a
         * single-depth science image; the variance (and thus resulting kernel)
         * generally converges after 1 iteration.  If
         * "constantVarianceWeighting" is requested in the Policy, no iterations
         * will be performed even if requested.
         */

        void build(
            lsst::afw::math::KernelList const& basisList
            );
        void build(
            lsst::afw::math::KernelList const& basisList,
            boost::shared_ptr<Eigen::MatrixXd> hMat
            );

       

    private:
        MaskedImagePtr _miToConvolvePtr;                    ///< Subimage around which you build kernel
        MaskedImagePtr _miToNotConvolvePtr;                 ///< Subimage around which you build kernel
        VariancePtr _varianceEstimate;                      ///< Estimate of the local variance
        lsst::pex::policy::Policy _policy;                  ///< Policy
        double _coreFlux;                                   ///< Mean S/N in the science image
        bool _isInitialized;                                ///< Has the kernel been built
        bool _useRegularization;                            ///< Use regularization?              
        bool _fitForBackground;

        /* best single raw kernel */
        boost::shared_ptr<StaticKernelSolution<PixelT> > _kernelSolutionOrig; ///< Original basis solution

        /* with Pca basis */
        boost::shared_ptr<StaticKernelSolution<PixelT> > _kernelSolutionPca;  ///< Most recent  solution

        void _buildKernelSolution(lsst::afw::math::KernelList const& basisList,
                                  boost::shared_ptr<Eigen::MatrixXd> hMat);
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
    boost::shared_ptr<KernelCandidate<PixelT> >
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
