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

#include <memory>
#include "Eigen/Core"

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/geom.h"
#include "lsst/daf/base.h"

namespace lsst {
namespace ip {
namespace diffim {

    /*
     * @brief Method used to solve for M and B
     */

    class KernelSolution {
    public:
        typedef std::shared_ptr<KernelSolution> Ptr;
        typedef lsst::afw::math::Kernel::Pixel PixelT;
        typedef lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel> ImageT;

        enum KernelSolvedBy {
            NONE          = 0,
            CHOLESKY_LDLT = 1,
            CHOLESKY_LLT  = 2,
            LU            = 3,
            EIGENVECTOR   = 4
        };

        enum ConditionNumberType {
            EIGENVALUE = 0,
            SVD        = 1
        };

        explicit KernelSolution(Eigen::MatrixXd mMat,
                                Eigen::VectorXd bVec,
                                bool fitForBackground);
        explicit KernelSolution(bool fitForBackground);
        explicit KernelSolution();

        virtual ~KernelSolution() {};
        virtual void solve();
        virtual void solve(Eigen::MatrixXd const& mMat,
                           Eigen::VectorXd const& bVec);
        KernelSolvedBy getSolvedBy() {return _solvedBy;}
        virtual double getConditionNumber(ConditionNumberType conditionType);
        virtual double getConditionNumber(Eigen::MatrixXd const& mMat, ConditionNumberType conditionType);

        inline Eigen::MatrixXd const& getM() {return _mMat;}
        inline Eigen::VectorXd const& getB() {return _bVec;}
        void printM() {std::cout << _mMat << std::endl;}
        void printB() {std::cout << _bVec << std::endl;}
        void printA() {std::cout << _aVec << std::endl;}
        inline int getId() const { return _id; }

    protected:
        int _id;                                                ///< Unique ID for object
        Eigen::MatrixXd _mMat;               ///< Derived least squares M matrix
        Eigen::VectorXd _bVec;               ///< Derived least squares B vector
        Eigen::VectorXd _aVec;               ///< Derived least squares solution matrix
        KernelSolvedBy _solvedBy;                               ///< Type of algorithm used to make solution
        bool _fitForBackground;                                 ///< Background terms included in fit
        static int _SolutionId;                                 ///< Unique identifier for solution

    };

    template <typename InputT>
    class StaticKernelSolution : public KernelSolution {
    public:
        typedef std::shared_ptr<StaticKernelSolution<InputT> > Ptr;

        StaticKernelSolution(lsst::afw::math::KernelList const& basisList,
                             bool fitForBackground);
        virtual ~StaticKernelSolution() {};

        /* Overrides KernelSolution */
        void solve();

        /* Used by RegularizedKernelSolution */
        virtual void build(lsst::afw::image::Image<InputT> const &templateImage,
                           lsst::afw::image::Image<InputT> const &scienceImage,
                           lsst::afw::image::Image<lsst::afw::image::VariancePixel> const &varianceEstimate);
        virtual std::shared_ptr<lsst::afw::math::Kernel> getKernel();
        virtual std::shared_ptr<lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel>> makeKernelImage();
        virtual double getBackground();
        virtual double getKsum();
        virtual std::pair<std::shared_ptr<lsst::afw::math::Kernel>, double> getSolutionPair();

    protected:
        Eigen::MatrixXd _cMat;               ///< K_i x R
        Eigen::VectorXd _iVec;               ///< Vectorized I
        Eigen::VectorXd _ivVec;              ///< Inverse variance

        std::shared_ptr<lsst::afw::math::Kernel> _kernel;                   ///< Derived single-object convolution kernel
        double _background;                                     ///< Derived differential background estimate
        double _kSum;                                           ///< Derived kernel sum

        void _setKernel();                                      ///< Set kernel after solution
        void _setKernelUncertainty();                           ///< Not implemented
    };


    template <typename InputT>
    class MaskedKernelSolution : public StaticKernelSolution<InputT> {
    public:
        typedef std::shared_ptr<MaskedKernelSolution<InputT> > Ptr;

        MaskedKernelSolution(lsst::afw::math::KernelList const& basisList,
                             bool fitForBackground);
        virtual ~MaskedKernelSolution() {};
        virtual void buildOrig(lsst::afw::image::Image<InputT> const &templateImage,
                               lsst::afw::image::Image<InputT> const &scienceImage,
                               lsst::afw::image::Image<lsst::afw::image::VariancePixel>
                               const &varianceEstimate,
                               lsst::afw::image::Mask<lsst::afw::image::MaskPixel> pixelMask);

        virtual void buildWithMask(lsst::afw::image::Image<InputT> const &templateImage,
                                   lsst::afw::image::Image<InputT> const &scienceImage,
                                   lsst::afw::image::Image<lsst::afw::image::VariancePixel>
                                   const &varianceEstimate,
                                   lsst::afw::image::Mask<lsst::afw::image::MaskPixel> const &pixelMask);

        virtual void buildSingleMaskOrig(lsst::afw::image::Image<InputT> const &templateImage,
                                         lsst::afw::image::Image<InputT> const &scienceImage,
                                         lsst::afw::image::Image<lsst::afw::image::VariancePixel>
                                         const &varianceEstimate,
                                         lsst::geom::Box2I maskBox);
    };



    template <typename InputT>
    class RegularizedKernelSolution : public StaticKernelSolution<InputT> {
    public:
        typedef std::shared_ptr<RegularizedKernelSolution<InputT> > Ptr;

        RegularizedKernelSolution(lsst::afw::math::KernelList const& basisList,
                                  bool fitForBackground,
                                  Eigen::MatrixXd const& hMat,
                                  lsst::daf::base::PropertySet const& ps
                                  );
        virtual ~RegularizedKernelSolution() {};
        void solve();
        double getLambda() {return _lambda;}
        double estimateRisk(double maxCond);

        /* Include additive term (_lambda * _hMat) in M matrix? */
        Eigen::MatrixXd getM(bool includeHmat = true);

    private:
        Eigen::MatrixXd const _hMat;               ///< Regularization weights
        double _lambda;                                         ///< Overall regularization strength
        lsst::daf::base::PropertySet::Ptr _ps;

        std::vector<double> _createLambdaSteps();
    };


    class SpatialKernelSolution : public KernelSolution {
    public:
        typedef std::shared_ptr<SpatialKernelSolution> Ptr;

        /* Creates a polynomial SpatialFunction */
        SpatialKernelSolution(lsst::afw::math::KernelList const& basisList,
                              lsst::afw::math::Kernel::SpatialFunctionPtr spatialKernelFunction,
                              lsst::afw::math::Kernel::SpatialFunctionPtr background,
                              lsst::daf::base::PropertySet const& ps
            );

        virtual ~SpatialKernelSolution() {};

        void addConstraint(float xCenter, float yCenter,
                           Eigen::MatrixXd const& qMat,
                           Eigen::VectorXd const& wVec);

        void solve();
        std::shared_ptr<lsst::afw::image::Image<lsst::afw::math::Kernel::Pixel>> makeKernelImage(lsst::geom::Point2D const& pos);
        std::pair<std::shared_ptr<lsst::afw::math::LinearCombinationKernel>,
                  lsst::afw::math::Kernel::SpatialFunctionPtr> getSolutionPair();

    private:
        lsst::afw::math::Kernel::SpatialFunctionPtr _spatialKernelFunction; ///< Spatial function for Kernel
        bool _constantFirstTerm;                                            ///< Is the first term constant

        std::shared_ptr<lsst::afw::math::LinearCombinationKernel> _kernel;   ///< Spatial convolution kernel
        lsst::afw::math::Kernel::SpatialFunctionPtr _background; ///< Spatial background model
        double _kSum;                                            ///< Derived kernel sum

        lsst::daf::base::PropertySet::Ptr _ps;                   ///< Config to control processing
        int _nbases;                                             ///< Number of basis functions
        int _nkt;                                                ///< Number of kernel terms
        int _nbt;                                                ///< Number of background terms
        int _nt;                                                 ///< Total number of terms

        void _setKernel();                                       ///< Set kernel after solution
        void _setKernelUncertainty();                            ///< Not implemented
    };

}}} // end of namespace lsst::ip::diffim

#endif
