// -*- lsst-c++ -*-
/**
 * @file KernelSolution.h
 *
 * @brief Declaration of classes to store the solution for convolution kernels
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_KERNELSOLUTION_H
#define LSST_IP_DIFFIM_KERNELSOLUTION_H

#include "boost/shared_ptr.hpp"
#include "Eigen/Core"

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"

namespace lsst { 
namespace ip { 
namespace diffim {

    /* 
     * @brief Method used to solve for M and B
     */

    class KernelSolution {
    public:
        typedef boost::shared_ptr<KernelSolution> Ptr;
        typedef lsst::afw::math::Kernel::Pixel PixelT;
        typedef lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel> ImageT;

        enum KernelSolvedBy {
            NONE          = 0,
            CHOLESKY_LDLT = 1,
            CHOLESKY_LLT  = 2,
            LU            = 3,
            EIGENVECTOR   = 4
        };

        explicit KernelSolution(boost::shared_ptr<Eigen::MatrixXd> mMat,
                                boost::shared_ptr<Eigen::VectorXd> bVec,
                                bool fitForBackground);
        explicit KernelSolution(bool fitForBackground);
        explicit KernelSolution();

        virtual ~KernelSolution() {};
        virtual void solve();
        virtual void solve(Eigen::MatrixXd _mMat, 
                           Eigen::VectorXd _bVec);

        inline boost::shared_ptr<Eigen::MatrixXd> getM() {return _mMat;}
        inline boost::shared_ptr<Eigen::VectorXd> getB() {return _bVec;}
        inline int getId() const { return _id; }

    protected:
        int _id;                                                ///< Unique ID for object
        boost::shared_ptr<Eigen::MatrixXd> _mMat;               ///< Derived least squares M matrix
        boost::shared_ptr<Eigen::VectorXd> _bVec;               ///< Derived least squares B vector
        boost::shared_ptr<Eigen::VectorXd> _aVec;               ///< Derived least squares solution matrix
        KernelSolvedBy _solvedBy;                               ///< Type of algorithm used to make solution
        bool _fitForBackground;                                 ///< Background terms included in fit
        static int _SolutionId;                                 ///< Unique identifier for solution

    };

    template <typename InputT>
    class StaticKernelSolution : public KernelSolution {
    public:
        typedef boost::shared_ptr<StaticKernelSolution<InputT> > Ptr;

        StaticKernelSolution(lsst::afw::math::KernelList const& basisList,
                             bool fitForBackground);
        virtual ~StaticKernelSolution() {};

        /* Overrides KernelSolution */
        void solve();

        /* Used by RegularizedKernelSolution */
        virtual void build(lsst::afw::image::Image<InputT> const &imageToConvolve,
                           lsst::afw::image::Image<InputT> const &imageToNotConvolve,
                           lsst::afw::image::Image<lsst::afw::image::VariancePixel> const &varianceEstimate);
        virtual lsst::afw::math::Kernel::Ptr getKernel();
        virtual ImageT::Ptr makeKernelImage();
        virtual double getBackground();
        virtual double getKsum();
        virtual std::pair<boost::shared_ptr<lsst::afw::math::Kernel>, double> getSolutionPair();

    protected:
        boost::shared_ptr<Eigen::MatrixXd> _cMat;               ///< K_i x R
        boost::shared_ptr<Eigen::VectorXd> _iVec;               ///< Vectorized I
        boost::shared_ptr<Eigen::VectorXd> _vVec;               ///< Variance

        lsst::afw::math::Kernel::Ptr _kernel;                   ///< Derived single-object convolution kernel
        double _background;                                     ///< Derived differential background estimate
        double _kSum;                                           ///< Derived kernel sum

        void _setKernel();                                      ///< Set kernel after solution
        void _setKernelUncertainty();                           ///< Not implemented
    };


    template <typename InputT>
    class RegularizedKernelSolution : public StaticKernelSolution<InputT> {
    public:
        typedef boost::shared_ptr<RegularizedKernelSolution<InputT> > Ptr;

        RegularizedKernelSolution(lsst::afw::math::KernelList const& basisList,
                                  bool fitForBackground,
                                  boost::shared_ptr<Eigen::MatrixXd> hMat,
                                  lsst::pex::policy::Policy policy
                                  );
        virtual ~RegularizedKernelSolution() {};
        void solve();
        double estimateRisk();
        double estimateGcv();

    private:
        boost::shared_ptr<Eigen::MatrixXd> _hMat;               ///< Regularization weights
        double _lambda;                                         ///< Overall regularization strength
        lsst::pex::policy::Policy _policy;
    };





    class SpatialKernelSolution : public KernelSolution {
    public:
        typedef boost::shared_ptr<SpatialKernelSolution> Ptr;

        SpatialKernelSolution(lsst::afw::math::KernelList const& basisList,
                              lsst::pex::policy::Policy policy
            );

        virtual ~SpatialKernelSolution() {};
        
        void addConstraint(float xCenter, float yCenter,
                           boost::shared_ptr<Eigen::MatrixXd> qMat,
                           boost::shared_ptr<Eigen::VectorXd> wVec);

        void solve();
        ImageT::Ptr makeKernelImage();
        std::pair<lsst::afw::math::LinearCombinationKernel::Ptr,
                  lsst::afw::math::Kernel::SpatialFunctionPtr> getKernelSolution();
        std::pair<lsst::afw::math::LinearCombinationKernel::Ptr,
                  lsst::afw::math::Kernel::SpatialFunctionPtr> getKernelUncertainty();
    private:
        lsst::afw::math::Kernel::SpatialFunctionPtr _spatialKernelFunction; ///< Spatial function for Kernel
        lsst::afw::math::Kernel::SpatialFunctionPtr _spatialBgFunction;     ///< Spatial function for Bg
        bool _constantFirstTerm;                                            ///< Is the first term constant

        lsst::afw::math::LinearCombinationKernel::Ptr _kernel;   ///< Spatial convolution kernel
        lsst::afw::math::Kernel::SpatialFunctionPtr _background; ///< Spatial background model
        double _kSum;                                           ///< Derived kernel sum

        lsst::afw::math::LinearCombinationKernel::Ptr _kernelErr;   ///< Kernel uncertainty
        lsst::afw::math::Kernel::SpatialFunctionPtr _backgroundErr; ///< Background uncertainty
        bool _errCalculated;                                        ///< Has the uncertainty been calculated?
        
        lsst::pex::policy::Policy _policy;
        int _nbases;
        int _nkt;
        int _nbt;
        int _nt;

        void _setKernel();
        void _setKernelUncertainty();
    };

}}} // end of namespace lsst::ip::diffim

#endif
