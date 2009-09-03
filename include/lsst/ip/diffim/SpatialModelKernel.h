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

#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/math/KernelFunctions.h>
#include <lsst/pex/policy/Policy.h>
#include <lsst/afw/detection/Footprint.h>
#include <lsst/sdqa/SdqaRating.h>

#include <lsst/ip/diffim/ImageSubtract.h>

namespace lsst {
namespace ip {
namespace diffim {

//    /** 
//     * 
//     * @brief Class stored in SpatialCells for spatial Kernel fitting
//     * 
//     * KernelCandidate is a single Kernel derived around a source.  We'll assign
//     * them to sets of SpatialCells; these sets will then be used to fit a
//     * spatial model to the PSF.
//     */    
//    template <typename ImageT>
//    class KernelCandidate : public lsst::afw::math::SpatialCellImageCandidate<ImageT> {
//        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getXCenter;
//        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getYCenter;
//        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getWidth;
//        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getHeight;
//    public: 
//        typedef boost::shared_ptr<KernelCandidate> Ptr;
//
//        /// Constructor
//        KernelCandidate(lsst::afw::detection::Source const& source, ///< The detected Source
//                        typename ImageT::ConstPtr parentImage       ///< The image wherein lie the Sources
//                        ) :
//            lsst::afw::math::SpatialCellImageCandidate<ImageT>(source.getXAstrom(), source.getYAstrom()),
//            _parentImage(parentImage),
//            _source(source),
//            _haveImage(false) {
//        }
//
//        /// Destructor
//        ~KernelCandidate() {};
//
//        /**
//         * Return Cell rating
//         * 
//         * @note Required method for use by SpatialCell
//         */
//        double getCandidateRating() const { return _source.getPsfFlux(); }
//
//        /// Return the original Source
//        lsst::afw::detection::Source const& getSource() const { return _source; }
//
//        typename ImageT::ConstPtr getImage() const;
//    private:
//        typename ImageT::ConstPtr _parentImage; // the %image that the Sources are found in
//        lsst::afw::detection::Source const _source; // the Source itself
//        bool mutable _haveImage;                    // do we have an Image to return?
//    };
//
//    /**
//     * Return a KernelCandidate of the right sort
//     *
//     * Cf. std::make_pair
//     */
//    template <typename ImageT>
//    typename KernelCandidate<ImageT>::Ptr
//    makeKernelCandidate(lsst::afw::detection::Source const& source, ///< The detected Source
//                        typename ImageT::ConstPtr image             ///< The image wherein lies the object
//                        ) {
//        
//        return typename KernelCandidate<ImageT>::Ptr(new KernelCandidate<ImageT>(source, image));
//    }
//
//    template<typename PixelT>
//    std::pair<lsst::afw::math::LinearCombinationKernel::PtrT, std::vector<double> >
//    createKernelFromCandidates(lsst::afw::math::SpatialCellSet const& kernelCells,
//                               int const nEigenComponents,
//                               int const spatialOrder,
//                               int const ksize,
//                               int const nPerCell=-1                                  
//                               );
//
//    template<typename PixelT>
//    std::pair<bool, double>
//    fitSpatialKernelFromCandidates(lsst::afw::math::Kernel *kernel,
//                                   lsst::afw::math::SpatialCellSet const& kernelCells, 
//                                   int const nStarPerCell=-1,
//                                   double const tolerance=1e-5);
    
    


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


}}}

#endif // LSST_IP_DIFFIM_SPATIALMODELKERNEL_H

