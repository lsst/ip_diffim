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
#include <lsst/afw/math/SpatialCell.h>
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

    template <typename ImageT, typename PixelT>
    class KernelCandidate : public lsst::afw::math::SpatialCellImageCandidate<ImageT> {
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getXCenter;
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getYCenter;
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getWidth;
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getHeight;
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::_image;
    public: 
        typedef boost::shared_ptr<KernelCandidate> Ptr;
        typedef boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > MaskedImagePtr;
	
        /** Constructor
         *
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
        double getCandidateRating();

        MaskedImagePtr getMiToConvolvePtr() {return _miToConvolvePtr;}
        MaskedImagePtr getMiToNotConvolvePtr() {return _miToNotConvolvePtr;}

        typename ImageT::ConstPtr getImage() const;
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
    template <typename ImageT, typename PixelT>
    typename KernelCandidate<ImageT,PixelT>::Ptr
    makeKernelCandidate(boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > const& miToConvolvePtr,
			boost::shared_ptr<lsst::afw::image::MaskedImage<PixelT> > const& miToNotConvolvePtr) {
        
        return typename KernelCandidate<ImageT,PixelT>::Ptr(new KernelCandidate<ImageT,PixelT>(miToConvolvePtr,
                                                                                               miToNotConvolvePtr));
    }

    template<typename PixelT>
    std::pair<lsst::afw::math::LinearCombinationKernel::PtrT, std::vector<double> >
    createPcaBasisFromCandidates(lsst::afw::math::SpatialCellSet const& psfCells,
                                 int const nEigenComponents,
                                 int const spatialOrder,
                                 int const nStarPerCell=-1);

    /*
    template<typename PixelT>
    std::pair<lsst::afw::math::LinearCombinationKernel::PtrT, std::vector<double> >
    returnKernelFromCandidates(lsst::afw::math::SpatialCellSet const& psfCells,
                               int const spatialOrder,
                               int const nStarPerCell=-1);
    */

    template<typename PixelT>
    std::pair<bool, double>
    fitSpatialKernelNonlinear(lsst::afw::math::Kernel *kernel,
                              lsst::afw::math::SpatialCellSet const& psfCells, 
                              int const nStarPerCell=-1,
                              double const tolerance=1e-5);
    
    template<typename PixelT>
    std::pair<bool, double>
    fitSpatialKernelLinear(lsst::afw::math::Kernel *kernel,
                           lsst::afw::math::SpatialCellSet const& psfCells, 
                           int const nStarPerCell=-1,
                           double const tolerance=1e-5);
    
#if 0

    /* OLD CODE BELOW */

    /** 
     * 
     * @brief Class used by SpatialModelCell for spatial Kernel fitting
     * 
     * A Kernel model is built for a given Footprint within a MaskedImage.  An
     * ensemble of Kernels, distributed evenly across the image using
     * SpatialModelCell, is used to fit for a spatial function.  If this Kernel
     * is a poor fit to the spatial function, another member of SpatialModelCell
     * will be used instead.
     *
     * This class needs to know how to build itself, meaning it requires the
     * basis functions used to create the Kernel, as well as the input images
     * that it is to compare.
     *
     * @see lsst/ip/diffim/SpatialModelCell.h for required methods
     */    
    template <typename ImageT>
    class SpatialModelKernel {
    public: 
        typedef boost::shared_ptr<SpatialModelKernel<ImageT> > Ptr;
        typedef std::vector<Ptr> SpatialModelKernelPtrList;
        typedef boost::shared_ptr<lsst::afw::image::MaskedImage<ImageT> > MaskedImagePtr; 

        /** Constructor
         *
         * @param fpPtr  Pointer to footprint of pixels used to build Kernel
         * @param miToConvolvePtr  Pointer to template image
         * @param miToNotConvolvePtr  Pointer to science image
         * @param kernelFunctor  Functor to build the PSF Mathching Kernel
         * @param policy  Policy for operations
         * @param build  Build upon construction?  Default is false.
         */
        SpatialModelKernel(lsst::afw::detection::Footprint::Ptr const& fpPtr,
                           MaskedImagePtr const& miToConvolvePtr,
                           MaskedImagePtr const& miToNotConvolvePtr,
                           boost::shared_ptr<PsfMatchingFunctor<ImageT> > const& kernelFunctor,
                           lsst::pex::policy::Policy const& policy,
                           bool build=false);

        /** Destructor
         */
        virtual ~SpatialModelKernel() {};


        /****
           Methods required for SpatialModelCell 
        */

        /** Execute the time-consuming process of building the local model
         * 
         * @note Required method for use by SpatialModelCell
         */
        bool buildModel();

        /** Return Cell rating
         * 
         * @note Required method for use by SpatialModelCell
         */
        double returnCellRating();

        /** Get its build status
         * 
         * @note Required method for use by SpatialModelCell
         */
        bool isBuilt() {return _isBuilt;};


        /****
         Build status 
        */

        /** Set its build status
         *
         * @param built  Boolean status of build
         */
        void setBuildStatus(bool built) {_isBuilt = built;};

        /** Get its build status
         */
        bool getBuildStatus() {return _isBuilt;};


        /****
           Quality status
        */

        /** Set its quality status
         */
        void setStatus(bool status      //!< Boolean status of build
                      ) {_isGood = status;};

        /** Get its quality status
         */
        bool getStatus() {return _isGood;};

        /** Get its quailty status
         */
        bool isGood() {return _isGood;};


        /****
           Position on the image
        */

        /** Set col centroid of Model; range -1 to 1
         *
         * @param colc  Column center
         */
        void setColc(double colc) {_colc = colc;};

        /** Get col centroid of Model; range -1 to 1
         */
        double getColc() {return _colc;};

        /** Set row centroid of Model; range -1 to 1
         *
         * @param rowc  Row center
         */
        void setRowc(double rowc) {_rowc = rowc;};

        /** Get row centroid of Model; range -1 to 1
         */
        double getRowc() {return _rowc;};


        /****
           Getters/setters
        */

        /** Get Footprint pointer for the Kernel model
         */
        lsst::afw::detection::Footprint::Ptr const& getFootprintPtr() const {return _fpPtr;};

        /** Get template's MaskedImage pointer for the Kernel model
         */
        MaskedImagePtr const& getMiToConvolvePtr() const {return _miToConvolvePtr;};

        /** Get image's MaskedImage pointer for the Kernel model
         */
        MaskedImagePtr const& getMiToNotConvolvePtr() const {return _miToNotConvolvePtr;};

        /** Set Kernel pointer associated with the Footprint; the core of this Model
         *
         * @param kPtr  pointer to the Kernel
         */
        void setKernelPtr(lsst::afw::math::Kernel::PtrT kPtr) {_kPtr = kPtr;};
        /** Get Kernel pointer associated with the Footprint
         */
        lsst::afw::math::Kernel::PtrT const& getKernelPtr() const {return _kPtr;};

        /** Set pointer associated with the uncertainty in the Kernel
         *
         * @param kPtr  pointer to the Kernel uncertainty; represent as a Kernel
         */
        void setKernelErrPtr(lsst::afw::math::Kernel::PtrT kPtr) {_kErrPtr = kPtr;};
        /** Get pointer associated with the uncertainty in the Kernel
         */
        lsst::afw::math::Kernel::PtrT getKernelErrPtr() {return _kErrPtr;};

        /** Set Kernel sum
         *
         * @param kSum  Kernel sum
         */
        void setKernelSum(double kSum) {_kSum = kSum;};
        /** Get Kernel sum
         */
        double getKernelSum() {return _kSum;};

        /** Set differential background value associated with the Kernel
         *
         * @param bg  Background value
         */
        void setBg(double bg) {_bg = bg;};
        /** Get differential background value associated with the Kernel
         */
        double getBg() {return _bg;};

        /** Set uncertainty in the differential background determination
         *
         * @param bgErr  Uncertainty in background 
         */
        void setBgErr(double bgErr) {_bgErr = bgErr;};
        /** Get uncertainty in the differential background determination
         */
        double getBgErr() {return _bgErr;};

        /** Set differential background value associated with the Kernel
         *
         * @param bg  Background value
         */
        void setBackground(double bg) {_bg = bg;};
        /** Get differential background value associated with the Kernel
         */
        double getBackground() {return _bg;};

        /** Set uncertainty in the differential background determination
         *
         * @param bgErr  Uncertainty in background 
         */
        void setBackgroundErr(double bgErr) {_bgErr = bgErr;};
        /** Get uncertainty in the differential background determination
         */
        double getBackgroundErr() {return _bgErr;};

        /** Set class instance associated with residuals in the derived difference image
         *
         * @param kStats  Pointer to instance of ImageStatistics class
         *
         * @note Ideally will be replaced by Sdqa
         *
         * @note Has to be a pointer since there is no empty constructor of FootprintFunctor
         */
        void setStats(typename ImageStatistics<ImageT>::Ptr kStats) {_kStats = kStats;};
        /** Get class instance associated with residuals in the derived difference image
         */
        typename ImageStatistics<ImageT>::Ptr getStats() {return _kStats;};

    private: 
        /** Objects needed to build itself; only initializable during construction
         */
        lsst::afw::detection::Footprint::Ptr _fpPtr; ///< Footprint containing pixels used to build Kernel
        MaskedImagePtr _miToConvolvePtr;             ///< Subimage around which you build kernel
        MaskedImagePtr _miToNotConvolvePtr;          ///< Subimage around which you build kernel
        typename PsfMatchingFunctor<ImageT>::Ptr _kFunctor; ///< Functor to build PSF matching kernel
        lsst::pex::policy::Policy _policy;           ///< Policy file for operations

        double _colc;     ///< Effective col position of model in overall image
        double _rowc;     ///< Effective col position of model in overall image

        /** Results from single Kernel fit
         */
        lsst::afw::math::Kernel::PtrT _kPtr;    ///< Kernel
        lsst::afw::math::Kernel::PtrT _kErrPtr; ///< Uncertainty in Kernel
        double _kSum;                                        ///< Kernel sum
        double _bg;                                          ///< Differential background value
        double _bgErr;                                       ///< Uncertainty in background
        typename ImageStatistics<ImageT>::Ptr _kStats; ///< Home-grown statistics; placeholder for Sdqa

        /** Status of model
         */
        bool _isBuilt;    ///< Model has been built
        bool _isGood;     ///< Passes local and/or Sdqa requirments


    }; // end of class

#endif
    
}}}

#endif // LSST_IP_DIFFIM_SPATIALMODELKERNEL_H

