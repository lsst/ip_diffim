// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Implementation of SpatialModelKernel class
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */
#include <boost/timer.hpp> 

#include <lsst/afw/image/Image.h>
#include <lsst/afw/image/ImagePca.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/math/FunctionLibrary.h>
#include <lsst/afw/math/Statistics.h>
#include <lsst/afw/detection/Footprint.h>

#include <lsst/pex/exceptions/Runtime.h>
#include <lsst/pex/policy/Policy.h>
#include <lsst/pex/logging/Trace.h>

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/QR>

#include <lsst/ip/diffim/SpatialModelKernel.h>

#define DEBUG_MATRIX 0

namespace afwMath        = lsst::afw::math;
namespace afwImage       = lsst::afw::image;
namespace pexLogging     = lsst::pex::logging; 
namespace pexExcept      = lsst::pex::exceptions; 
namespace pexPolicy      = lsst::pex::policy; 

namespace lsst {
namespace ip {
namespace diffim {

template <typename PixelT>
KernelCandidate<PixelT>::ImageT::ConstPtr KernelCandidate<PixelT>::getImage() const {
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Kernel to make KernelCandidate Image from");
    }
    return _image;
}

template <typename PixelT>
KernelCandidate<PixelT>::ImageT::Ptr KernelCandidate<PixelT>::copyImage() const {
    return typename KernelCandidate<PixelT>::ImageT::Ptr(new typename KernelCandidate<PixelT>::ImageT(*getImage(), true));
    /*
    typename KernelCandidate<PixelT>::ImageT::Ptr imcopy(
        new typename KernelCandidate<PixelT>::ImageT(*_image, true)
        );
    return imcopy;
    */
}

  
template <typename PixelT>
void KernelCandidate<PixelT>::setKernel(lsst::afw::math::Kernel::Ptr kernel) {
    _kernel     = kernel; 
    _haveKernel = true;

    setWidth(_kernel->getWidth());
    setHeight(_kernel->getHeight());
    
    typename KernelCandidate<PixelT>::ImageT::Ptr image (
        new typename KernelCandidate<PixelT>::ImageT(_kernel->getDimensions())
        );
    _kSum  = _kernel->computeImage(*image, false);                    
    _image = image;
}

template <typename PixelT>
afwMath::Kernel::Ptr KernelCandidate<PixelT>::getKernel() const {
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Kernel for KernelCandidate");
    }
    return _kernel;
}

template <typename PixelT>
double KernelCandidate<PixelT>::getBackground() const {
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Kernel for KernelCandidate");
    }
    return _background;
}

template <typename PixelT>
double KernelCandidate<PixelT>::getKsum() const {
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Kernel for KernelCandidate");
    }
    return _kSum;
}

template <typename PixelT>
lsst::afw::image::MaskedImage<PixelT> KernelCandidate<PixelT>::returnDifferenceImage() {
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Kernel for KernelCandidate");
    }
    return returnDifferenceImage(_kernel, _background);
}

template <typename PixelT>
lsst::afw::image::MaskedImage<PixelT> KernelCandidate<PixelT>::returnDifferenceImage(
    lsst::afw::math::Kernel::Ptr kernel,
    double background
    ) {
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Kernel for KernelCandidate");
    }
    
    /* Make diffim and set chi2 from result */
    lsst::afw::image::MaskedImage<PixelT> diffim = convolveAndSubtract(*_miToConvolvePtr,
                                                                       *_miToNotConvolvePtr,
                                                                       *kernel,
                                                                       background);
    return diffim;

}

/* Currently implemented in an anonymous namespace; consider moving to ::detail
   namespace so that they can be unit tested */
namespace {
    template<typename PixelT>
    class KernelSumVisitor : public afwMath::CandidateVisitor {
        typedef afwImage::Image<lsst::afw::math::Kernel::Pixel> ImageT;
    public:
        enum Mode {AGGREGATE = 0, REJECT = 1};

        KernelSumVisitor(lsst::pex::policy::Policy const& policy) :
            afwMath::CandidateVisitor(),
            _mode(AGGREGATE),
            _kSums(std::vector<double>()),
            _kSumMean(0.),
            _kSumStd(0.),
            _dkSumMax(0.),
            _kSumNpts(0),
            _nRejected(0),
            _policy(policy) {}

        void setMode(Mode mode) {_mode = mode;}
        int getNRejected() {
            return _nRejected;
        }
        void reset() {
            _kSums.clear();
            _nRejected = 0;
        }
        void processCandidate(afwMath::SpatialCellCandidate *candidate) {
            KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
            if (kCandidate == NULL) {
                throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                                  "Failed to cast SpatialCellCandidate to KernelCandidate");
            }
            pexLogging::TTrace<6>("lsst.ip.diffim.KernelSumVisitor.processCandidate", 
                                  "Processing candidate %d, mode %d", kCandidate->getId(), _mode);

            /* Grab all kernel sums and look for outliers */
            if (_mode == AGGREGATE) {
                _kSums.push_back(kCandidate->getKsum());
            }
            else if (_mode == REJECT) {
                if ( fabs(kCandidate->getKsum() - _kSumMean) > (_dkSumMax) ) {
                    kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
                    pexLogging::TTrace<5>("lsst.ip.diffim.KernelSumVisitor.processCandidate", 
                                          "Rejecting candidate %d due to bad source kernel sum : (%.2f)",
                                          kCandidate->getId(),
                                          kCandidate->getKsum());
                    _nRejected += 1;
                }
            }
            
        }
        
        void processKsumDistribution() {
            afwMath::Statistics stats = afwMath::makeStatistics(_kSums, afwMath::NPOINT | afwMath::MEANCLIP | afwMath::STDEVCLIP); 
            _kSumMean = stats.getValue(afwMath::MEANCLIP);
            _kSumStd  = stats.getValue(afwMath::STDEVCLIP);
            _kSumNpts = static_cast<int>(stats.getValue(afwMath::NPOINT));
            _dkSumMax = _policy.getDouble("maxKsumSigma") * _kSumStd;
            pexLogging::TTrace<3>("lsst.ip.diffim.KernelSumVisitor.processCandidate", 
                                  "Kernel Sum Distribution : %.3f +/- %.3f (%d points)", _kSumMean, _kSumStd, _kSumNpts);
        }
 
        
    private:
        Mode _mode;
        std::vector<double> _kSums;
        double _kSumMean;
        double _kSumStd;
        double _dkSumMax;
        int    _kSumNpts;
        int    _nRejected;
        lsst::pex::policy::Policy _policy;
    };    

    template<typename PixelT>
    class SetPcaImageVisitor : public afwMath::CandidateVisitor {
        typedef afwImage::Image<lsst::afw::math::Kernel::Pixel> ImageT;
    public:

        SetPcaImageVisitor(afwImage::ImagePca<ImageT> *imagePca // Set of Images to initialise
            ) :
            afwMath::CandidateVisitor(),
            _imagePca(imagePca) {}
        
        // Called by SpatialCellSet::visitCandidates for each Candidate
        void processCandidate(afwMath::SpatialCellCandidate *candidate) {
            KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
            if (kCandidate == NULL) {
                throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                                  "Failed to cast SpatialCellCandidate to KernelCandidate");
            }
            
            try {
                /* 
                   We don't necessarily want images with larger kernel sums to
                   have more weight.  Each kernel should have constant weight in
                   the Pca.  For simplicity we will also scale them to have the
                   same kernel sum, 1.0, and send to ImagePca that the flux is
                   1.0.
                */
                ImageT::Ptr kImage = kCandidate->copyImage();
                *kImage           /= kCandidate->getKsum();
                
                _imagePca->addImage(kImage, 1.0);
            } catch(lsst::pex::exceptions::LengthErrorException &e) {
                return;
            }
        }
    private:
        afwImage::ImagePca<ImageT> *_imagePca; 
    };


    template<typename PixelT>
    class BuildSingleKernelVisitor : public afwMath::CandidateVisitor {
        typedef afwImage::MaskedImage<PixelT> MaskedImageT;
    public:
        BuildSingleKernelVisitor(PsfMatchingFunctor<PixelT> &kFunctor,
                                 lsst::pex::policy::Policy const& policy) :
            afwMath::CandidateVisitor(),
            _kFunctor(kFunctor),
            _policy(policy),
            _imstats(ImageStatistics<PixelT>()),
            _setCandidateKernel(true),
            _skipBuilt(true),
            _nRejected(0)
            {}

        /* 
         * This functionality allows the user to not set the "kernel" and thus
         * "image" values of the KernelCandidate.  When running a PCA fit on the
         * kernels, we want to keep the delta-function representation of the raw
         * kernel which are used to derive the eigenBases, while still being
         * able to modify the _M and _B matrices with linear fits to the
         * eigenBases themselves.
         */
        void setCandidateKernel(bool set) {_setCandidateKernel = set;}

        /* 
           Don't reprocess candidate if its already been build.  The use
           case for this functionality is : when iterating over all Cells
           and rejecting bad Kernels, we need to re-visit *all* Cells to
           build the next candidate in the list.  Without this flag we would
           unncessarily re-build all the good Kernels.
        */
        void setSkipBuilt(bool skip)      {_skipBuilt = skip;}

        /* 
           Since this is the base class that builds a kernel, we need to make
           sure that the current Kernel in the Cell is initialized.  To do this,
           if we set something afwMath::SpatialCellCandidate::BAD we have to go
           back over the Cells and build the next Candidate, until we are out of
           Candidates (in which case this will get called on no Cells and
           nRejected=0) or all current Candidates are
           afwMath::SpatialCellCandidate::GOOD.
        */
        void reset()          {_nRejected = 0;}
        int  getNRejected()   {return _nRejected;}


        void processCandidate(afwMath::SpatialCellCandidate *candidate) {
            KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
            if (kCandidate == NULL) {
                throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                                  "Failed to cast SpatialCellCandidate to KernelCandidate");
            }

            if (_skipBuilt and kCandidate->hasKernel()) {
                return;
            }

            pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                  "Processing candidate %d", kCandidate->getId());
            
            /* Estimate of the variance */
            MaskedImageT var = MaskedImageT(*(kCandidate->getMiToNotConvolvePtr()), true);
            if (_policy.getBool("constantVarianceWeighting")) {
                /* Constant variance weighting */
                *var.getVariance() = 1;
            }
            else {
                /* Variance estimate is the straight difference */
                var -= *(kCandidate->getMiToConvolvePtr());
            }
            
            /* Build its kernel here */
            try {
                _kFunctor.apply(*(kCandidate->getMiToConvolvePtr()->getImage()),
                                *(kCandidate->getMiToNotConvolvePtr()->getImage()),
                                *(var.getVariance()),
                                _policy);
            } catch (lsst::pex::exceptions::Exception &e) {
                LSST_EXCEPT_ADD(e, "Unable to calculate Kernel");
                throw e;
            }
            
            /* 
               Sometimes you do not want to override the kernel; e.g. on a
               second fitting loop after the results of the first fitting loop
               are used to define a PCA basis
            */
            std::pair<boost::shared_ptr<lsst::afw::math::Kernel>, double> KB = _kFunctor.getKernel();
            if (_setCandidateKernel) {
                kCandidate->setKernel(KB.first);
                kCandidate->setBackground(KB.second);
            }

            /* 
             * However you *always* need to reset M and B since these are used *
             * in the spatial fitting
             */
            std::pair<boost::shared_ptr<Eigen::MatrixXd>, boost::shared_ptr<Eigen::VectorXd> > MB = _kFunctor.getAndClearMB();
            kCandidate->setM(MB.first);
            kCandidate->setB(MB.second);

            /* 
             * Make diffim and set chi2 from result.  Note that you need to send
             * the newly-derived kernel and background in the case that
             * _setCandidateKernel = false.
            */
            MaskedImageT diffim = kCandidate->returnDifferenceImage(KB.first, KB.second);
            
            /* 
             * Remake the kernel using the first iteration difference image
             * variance as a better estimate of the true diffim variance
             */
            if (_policy.getBool("iterateSingleKernel")) {
                try {
                    _kFunctor.apply(*(kCandidate->getMiToConvolvePtr()->getImage()),
                                    *(kCandidate->getMiToNotConvolvePtr()->getImage()),
                                    *(diffim.getVariance()),
                                    _policy);
                } catch (lsst::pex::exceptions::Exception &e) {
                    LSST_EXCEPT_ADD(e, "Unable to recalculate Kernel");
                    throw e;
                }
                KB = _kFunctor.getKernel();
                if (_setCandidateKernel) {
                    kCandidate->setKernel(KB.first);
                    kCandidate->setBackground(KB.second);
                }
                MB = _kFunctor.getAndClearMB();
                kCandidate->setM(MB.first);
                kCandidate->setB(MB.second);
                diffim = kCandidate->returnDifferenceImage(KB.first, KB.second);                
            }

            _imstats.apply(diffim);
            kCandidate->setChi2(_imstats.getVariance());

            /* When using a Pca basis, we don't reset the kernel or background,
               so we need to evaluate these locally for the Trace */
            afwImage::Image<double> kImage(KB.first->getDimensions());
            double kSum = KB.first->computeImage(kImage, false);
            double background = KB.second;

            pexLogging::TTrace<4>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                  "Chi2 = %.2f", kCandidate->getChi2());
            pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate",
                                  "X = %.2f Y = %.2f",
                                  kCandidate->getXCenter(), 
                                  kCandidate->getYCenter());
            pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate",
                                  "Kernel Sum = %.3f", kSum);
            pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate",
                                  "Background = %.3f", background);
            pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate",
                                  "Diffim residuals = %.2f +/- %.2f sigma",
                                  _imstats.getMean(),
                                  _imstats.getRms());

            if (fabs(_imstats.getMean()) > _policy.getDouble("candidateResidualMeanMax")) {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
                pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                      "Rejecting due to bad source kernel mean residuals : |%.2f| > %.2f",
                                      _imstats.getMean(),
                                      _policy.getDouble("candidateResidualMeanMax"));
                _nRejected += 1;
            }
            else if (_imstats.getRms() > _policy.getDouble("candidateResidualStdMax")) {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
                pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                      "Rejecting due to bad source kernel residual rms : %.2f > %.2f",
                                      _imstats.getRms(),
                                      _policy.getDouble("candidateResidualStdMax"));
                _nRejected += 1;
            }
            else {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::GOOD);
                pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor.processCandidate", 
                                      "Source kernel OK");
            }
        }
    private:
        PsfMatchingFunctor<PixelT> _kFunctor;
        lsst::pex::policy::Policy _policy;
        ImageStatistics<PixelT> _imstats;
        bool _setCandidateKernel;
        bool _skipBuilt;
        int _nRejected;
    };


    template<typename PixelT>
    class AssessSpatialKernelVisitor : public afwMath::CandidateVisitor {
        typedef afwImage::MaskedImage<PixelT> MaskedImageT;
    public:
        AssessSpatialKernelVisitor(afwMath::LinearCombinationKernel::Ptr spatialKernel,
                                   afwMath::Kernel::SpatialFunctionPtr spatialBackground,
                                   lsst::pex::policy::Policy const& policy) :
            afwMath::CandidateVisitor(),
            _spatialKernel(spatialKernel),
            _spatialBackground(spatialBackground),
            _policy(policy),
            _imstats(ImageStatistics<PixelT>()),
            _nRejected(0) {}

        void processCandidate(afwMath::SpatialCellCandidate *candidate) {
            KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
            if (kCandidate == NULL) {
                throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                                  "Failed to cast SpatialCellCandidate to KernelCandidate");
            }
            if (!(kCandidate->hasKernel())) {
                pexLogging::TTrace<3>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                      "Cannot process candidate %d, continuing", kCandidate->getId());
                return;
            }
            
            pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                  "Processing candidate %d", kCandidate->getId());
            
            /* Make diffim and set chi2 from result */

            /* 
               Note - this is a hack until the Kernel API is upgraded by the
               Davis crew.  I need a "local" version of the spatially varying
               Kernel
            */
            afwImage::Image<double> kImage(_spatialKernel->getDimensions());
            double kSum = _spatialKernel->computeImage(kImage, false, 
                                                       kCandidate->getXCenter(),
                                                       kCandidate->getYCenter());
            boost::shared_ptr<afwMath::Kernel>
                kernelPtr( new afwMath::FixedKernel(kImage) );
            /* </hack> */

            double background = (*_spatialBackground)(kCandidate->getXCenter(), kCandidate->getYCenter());
            
            MaskedImageT diffim = kCandidate->returnDifferenceImage(kernelPtr, background);
            _imstats.apply(diffim);
            kCandidate->setChi2(_imstats.getVariance());
            
            pexLogging::TTrace<4>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                  "Chi2 = %.2f", kCandidate->getChi2());
            pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate",
                                  "X = %.2f Y = %.2f",
                                  kCandidate->getXCenter(), 
                                  kCandidate->getYCenter());
            pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate",
                                  "Kernel Sum = %.3f", kSum);
            pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate",
                                  "Background = %.3f", background);
            pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate",
                                  "Diffim residuals = %.2f +/- %.2f sigma",
                                  _imstats.getMean(),
                                  _imstats.getRms());
            
            if (fabs(_imstats.getMean()) > _policy.getDouble("candidateResidualMeanMax")) {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
                pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                      "Rejecting due to bad spatial kernel mean residuals : |%.2f| > %.2f",
                                      _imstats.getMean(),
                                      _policy.getDouble("candidateResidualMeanMax"));
                _nRejected += 1;
            }
            else if (_imstats.getRms() > _policy.getDouble("candidateResidualStdMax")) {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
                pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                      "Rejecting due to bad spatial kernel residual rms : %.2f > %.2f",
                                      _imstats.getRms(),
                                      _policy.getDouble("candidateResidualStdMax"));
                _nRejected += 1;
            }
            else {
                kCandidate->setStatus(afwMath::SpatialCellCandidate::GOOD);
                pexLogging::TTrace<5>("lsst.ip.diffim.AssessSpatialKernelVisitor.processCandidate", 
                                      "Spatial kernel OK");
            }
        }

        int getNRejected() {
            return _nRejected;
        }

    private:
        afwMath::LinearCombinationKernel::Ptr _spatialKernel;
        afwMath::Kernel::SpatialFunctionPtr _spatialBackground;
        lsst::pex::policy::Policy _policy;
        ImageStatistics<PixelT> _imstats;
        int _nRejected;
    };


    template<typename PixelT>
    class BuildSpatialKernelVisitor : public afwMath::CandidateVisitor {
    public:
        BuildSpatialKernelVisitor(
            afwMath::KernelList basisList, ///< Basis functions used in the fit
            int const spatialKernelOrder,   ///< Order of spatial kernel variation (cf. lsst::afw::math::PolynomialFunction2)
            int const spatialBgOrder,       ///< Order of spatial bg variation (cf. lsst::afw::math::PolynomialFunction2)
            lsst::pex::policy::Policy policy
            ) :
            afwMath::CandidateVisitor(),
            _basisList(basisList),
            _M(Eigen::MatrixXd()),
            _B(Eigen::VectorXd()),
            _Soln(Eigen::VectorXd()),
            _spatialKernelOrder(spatialKernelOrder),
            _spatialBgOrder(spatialBgOrder),
            _spatialKernelFunction( new afwMath::PolynomialFunction2<double>(spatialKernelOrder) ),
            _spatialBgFunction( new afwMath::PolynomialFunction2<double>(spatialBgOrder) ),
            _nbases(basisList.size()),
            _policy(policy),
            _constantFirstTerm(false){

            /* 
               NOTE : The variable _constantFirstTerm allows that the first
               component of the basisList has no spatial variation.  This is
               useful to conserve the kernel sum across the image.  There are 2
               options to enable this : we can make the input matrices/vectors
               smaller by (_nkt-1), or we can create _M and _B with empty values
               for the first component's spatial terms.  The latter should cause
               problems for the matrix math even though its probably more
               readable, so we go with the former.
            */
            if ( (_policy.getString("kernelBasisSet") == "alard-lupton") || _policy.getBool("usePcaForSpatialKernel")) {
	       _constantFirstTerm = true;
            }

            /* Bookeeping terms */
            _nkt = _spatialKernelFunction->getParameters().size();
            _nbt = _spatialBgFunction->getParameters().size();
            if (_constantFirstTerm) {
 	       _nt = (_nbases-1)*_nkt+1 + _nbt;
            } else {
	       _nt = _nbases*_nkt + _nbt;
            }
	    _M.resize(_nt, _nt);
	    _B.resize(_nt);
	    _M.setZero();
	    _B.setZero();

            pexLogging::TTrace<5>("lsst.ip.diffim.LinearSpatialFitVisitor", 
                                  "Initializing with size %d %d %d and constant first term = %s",
                                  _nkt, _nbt, _nt,
                                  _constantFirstTerm ? "true" : "false");
        }
        
        void processCandidate(afwMath::SpatialCellCandidate *candidate) {
            KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
            if (kCandidate == NULL) {
                throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                                  "Failed to cast SpatialCellCandidate to KernelCandidate");
            }
            if (!(kCandidate->hasKernel())) {
                pexLogging::TTrace<3>("lsst.ip.diffim.BuildSpatialKernelVisitor.processCandidate", 
                                      "Cannot process candidate %d, continuing", kCandidate->getId());
                return;
            }

            pexLogging::TTrace<6>("lsst.ip.diffim.BuildSpatialKernelVisitor.processCandidate", 
                                  "Processing candidate %d", kCandidate->getId());
            
            /* Calculate P matrices */
            /* Pure kernel terms */
            std::vector<double> paramsK = _spatialKernelFunction->getParameters();
            for (unsigned int idx = 0; idx < _nkt; idx++) { paramsK[idx] = 0.0; }
            Eigen::VectorXd Pk(_nkt);
            for (unsigned int idx = 0; idx < _nkt; idx++) {
                paramsK[idx] = 1.0;
                _spatialKernelFunction->setParameters(paramsK);
                Pk(idx) = (*_spatialKernelFunction)( kCandidate->getXCenter(), 
                                                     kCandidate->getYCenter() );
                paramsK[idx] = 0.0;
            }
            Eigen::MatrixXd PkPkt = (Pk * Pk.transpose());
            
            /* Pure background terms */
            std::vector<double> paramsB = _spatialBgFunction->getParameters();
            for (unsigned int idx = 0; idx < _nbt; idx++) { paramsB[idx] = 0.0; }
            Eigen::VectorXd Pb(_nbt);
            for (unsigned int idx = 0; idx < _nbt; idx++) {
                paramsB[idx] = 1.0;
                _spatialBgFunction->setParameters(paramsB);
                Pb(idx) = (*_spatialBgFunction)( kCandidate->getXCenter(), 
                                                 kCandidate->getYCenter() );
                paramsB[idx] = 0.0;
            }
            Eigen::MatrixXd PbPbt = (Pb * Pb.transpose());

            /* Cross terms */
            Eigen::MatrixXd PkPbt = (Pk * Pb.transpose());

            if (DEBUG_MATRIX) {
                std::cout << "Spatial weights" << std::endl;
                std::cout << "PkPkt " << PkPkt << std::endl;
                std::cout << "PbPbt " << PbPbt << std::endl;
                std::cout << "PkPbt " << PkPbt << std::endl;
            }
            
            /* Add each candidate to the M, B matrix */
            boost::shared_ptr<Eigen::MatrixXd> Q = kCandidate->getM();
            boost::shared_ptr<Eigen::VectorXd> W = kCandidate->getB();

            if (DEBUG_MATRIX) {
                std::cout << "Spatial matrix inputs" << std::endl;
                std::cout << "M " << (*Q) << std::endl;
                std::cout << "B " << (*W) << std::endl;
            }

	    /* first index to start the spatial blocks; default=0 for non-constant first term */
	    unsigned int m0 = 0; 
	    /* how many rows/cols to adjust the matrices/vectors; default=0 for non-constant first term */
	    unsigned int dm = 0; 
	    /* where to start the background terms; default always true */
	    unsigned int mb = _nt - _nbt;

            if (_constantFirstTerm) {
	       m0 = 1;       /* we need to manually fill in the first (non-spatial) terms below */
	       dm = _nkt-1;  /* need to shift terms due to lack of spatial variation in first term */

	       _M(0, 0) += (*Q)(0,0);
	       for(unsigned int m2 = 1; m2 < _nbases; m2++)  {
		  _M.block(0, m2*_nkt-dm, 1, _nkt) += (*Q)(0,m2) * Pk.transpose();
	       }
	       
	       _M.block(0, mb, 1, _nbt) += (*Q)(0,_nbases) * Pb.transpose();
	       _B(0) += (*W)(0);
	    }
	    
	    /* Fill in the spatial blocks */
	    for(unsigned int m1 = m0; m1 < _nbases; m1++)  {
	       /* Diagonal kernel-kernel term; only use upper triangular part of PkPkt */
	       _M.block(m1*_nkt-dm, m1*_nkt-dm, _nkt, _nkt) += (*Q)(m1,m1) * PkPkt.part<Eigen::UpperTriangular>();

	       /* Kernel-kernel terms */
	       for(unsigned int m2 = m1+1; m2 < _nbases; m2++)  {
		  _M.block(m1*_nkt-dm, m2*_nkt-dm, _nkt, _nkt) += (*Q)(m1,m2) * PkPkt;
	       }
	       
	       /* Kernel cross terms with background */
	       _M.block(m1*_nkt-dm, mb, _nkt, _nbt) += (*Q)(m1,_nbases) * PkPbt;
               
	       /* B vector */
	       _B.segment(m1*_nkt-dm, _nkt) += (*W)(m1) * Pk;
	    }
            
	    /* Background-background terms only */
	    _M.block(mb, mb, _nbt, _nbt) += (*Q)(_nbases,_nbases) * PbPbt.part<Eigen::UpperTriangular>();
	    _B.segment(mb, _nbt)         += (*W)(_nbases) * Pb;

#if 0
            if (_constantFirstTerm) {
	       /* Fill in matrices */
	       unsigned int dm = _nkt-1;
	       unsigned int mb = (_nbases-1)*_nkt+1;
	       
	       /* Zeroth term, m1=0 */

	       /* The rest of the terms */
	       for(unsigned int m1 = 1; m1 < _nbases; m1++)  {
		  /* Diagonal kernel-kernel term; only use upper triangular part of PkPkt */
		  _M.block(m1*_nkt-dm, m1*_nkt-dm, _nkt, _nkt) += (*Q)(m1,m1) * PkPkt.part<Eigen::UpperTriangular>();
		  
		  /* Kernel-kernel terms */
		  for(unsigned int m2 = m1+1; m2 < _nbases; m2++)  {
		     _M.block(m1*_nkt-dm, m2*_nkt-dm, _nkt, _nkt) += (*Q)(m1,m2) * PkPkt;
		  }
		  
		  /* Kernel cross terms with background */
		  _M.block(m1*_nkt-dm, mb, _nkt, _nbt) += (*Q)(m1,_nbases) * PkPbt;

		  /* B vector */
		  _B.segment(m1*_nkt-dm, _nkt) += (*W)(m1) * Pk;
	       } 

	       /* Background-background terms only */
	       _M.block(mb, mb, _nbt, _nbt) += (*Q)(_nbases,_nbases) * PbPbt.part<Eigen::UpperTriangular>();
	       _B.segment(mb, _nbt)         += (*W)(_nbases) * Pb;
            }
            else {
                /* Fill in matrices */
	        unsigned int mb = _nbases*_nkt;
                for(unsigned int m1 = 0; m1 < _nbases; m1++)  {
		    /* Diagonal kernel-kernel term; only use upper triangular part of PkPkt */
 		    _M.block(m1*_nkt, m1*_nkt, _nkt, _nkt) += (*Q)(m1,m1) * PkPkt.part<Eigen::UpperTriangular>();
                    
                    /* Kernel-kernel terms */
                    for(unsigned int m2 = m1+1; m2 < _nbases; m2++)  {
		       _M.block(m1*_nkt, m2*_nkt, _nkt, _nkt) += (*Q)(m1,m2) * PkPkt;
                    }
	            
                    /* Kernel cross terms with background */
                    _M.block(m1*_nkt, mb, _nkt, _nbt) += (*Q)(m1,_nbases) * PkPbt;
                    
                    /* B vector */
                    _B.segment(m1*_nkt, _nkt) += (*W)(m1) * Pk;
                }
                
                /* Background-background terms only */
                _M.block(mb, mb, _nbt, _nbt) += (*Q)(_nbases,_nbases) * PbPbt.part<Eigen::UpperTriangular>();
                _B.segment(mb, _nbt)         += (*W)(_nbases) * Pb;
            }
#endif
            if (DEBUG_MATRIX) {
                std::cout << "Spatial matrix outputs" << std::endl;
                std::cout << "_M " << _M << std::endl;
                std::cout << "_B " << _B << std::endl;
            }

        }
        
        void solveLinearEquation() {
            boost::timer t;
            t.restart();

            /* Fill in the other half of _M */
            for (unsigned int i = 0; i < _nt; i++) {
                for (unsigned int j = i+1; j < _nt; j++) {
                    _M(j,i) = _M(i,j);
                }
            }
            _Soln = Eigen::VectorXd::Zero(_nt);
            
            if (DEBUG_MATRIX) {
                std::cout << "Solving for _M:" << std::endl;
                std::cout << _M << std::endl;
                std::cout << _B << std::endl;
            }

            if (!( _M.ldlt().solve(_B, &_Soln) )) {
                pexLogging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.solveLinearEquation", 
                                      "Unable to determine kernel via Cholesky LDL^T");
                if (!( _M.llt().solve(_B, &_Soln) )) {
                    pexLogging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.solveLinearEquation", 
                                          "Unable to determine kernel via Cholesky LL^T");
                    if (!( _M.lu().solve(_B, &_Soln) )) {
                        pexLogging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.solveLinearEquation", 
                                              "Unable to determine kernel via LU");
                        // LAST RESORT
                        try {
                            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eVecValues(_M);
                            Eigen::MatrixXd const& R = eVecValues.eigenvectors();
                            Eigen::VectorXd eValues  = eVecValues.eigenvalues();
                            
                            for (int i = 0; i != eValues.rows(); ++i) {
                                if (eValues(i) != 0.0) {
                                    eValues(i) = 1.0/eValues(i);
                                }
                            }
                            
                            _Soln = R*eValues.asDiagonal()*R.transpose()*_B;
                        } catch (pexExcept::Exception& e) {
                            pexLogging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.solveLinearEquation", 
                                                  "Unable to determine kernel via eigen-values");
                            
                            throw LSST_EXCEPT(pexExcept::Exception, 
                                              "Unable to determine kernel solution in SpatialModelKernel::solveLinearEquation");
                        }
                    }
                }
            }

            if (DEBUG_MATRIX) {
                std::cout << "Solution:" << std::endl;
                std::cout << _Soln << std::endl;
            }

            double time = t.elapsed();
            pexLogging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.solveLinearEquation", 
                                  "Compute time to do spatial matrix math : %.2f s", time);
        }

        
        Eigen::VectorXd getSolution() {
            return _Soln; 
        }
        
        std::pair<afwMath::LinearCombinationKernel::Ptr, afwMath::Kernel::SpatialFunctionPtr> getSpatialModel() {
            /* Set up kernel */
            std::vector<afwMath::Kernel::SpatialFunctionPtr> spatialFunctionList;
            for (unsigned int i = 0; i < _nbases; i++) {
                afwMath::Kernel::SpatialFunctionPtr spatialFunction(_spatialKernelFunction->copy());
                spatialFunctionList.push_back(spatialFunction);
            }
            afwMath::LinearCombinationKernel::Ptr spatialKernel(new afwMath::LinearCombinationKernel(_basisList, spatialFunctionList));
            
            /* Set up background */
            afwMath::Kernel::SpatialFunctionPtr bgFunction(_spatialBgFunction->copy());
            
            /* Set the kernel coefficients */
            std::vector<std::vector<double> > kCoeffs;
            kCoeffs.reserve(_nbases);
            for (unsigned int i = 0, idx = 0; i < _nbases; i++) {
                kCoeffs.push_back(std::vector<double>(_nkt));

		/* Deal with the possibility the first term doesn't vary spatially */
		if ( (i == 0) && (_constantFirstTerm) ) {
		   kCoeffs[i][0] = _Soln[idx++];
		}
		else {
		   for (unsigned int j = 0; j < _nkt; j++) {
		      kCoeffs[i][j] = _Soln[idx++];
		   }
		}
            }
            
            /* Set the background coefficients */
            std::vector<double> bgCoeffs(_nbt);
            for (unsigned int i = 0; i < _nbt; i++) {
                bgCoeffs[i] = _Soln[_nt - _nbt + i];
            }
            
            spatialKernel->setSpatialParameters(kCoeffs);
            bgFunction->setParameters(bgCoeffs);
            
            return std::make_pair(spatialKernel, bgFunction);
        }

    private:
        afwMath::KernelList _basisList;
        Eigen::MatrixXd _M;       ///< Least squares matrix
        Eigen::VectorXd _B;       ///< Least squares vector
        Eigen::VectorXd _Soln;    ///< Least squares solution
        int const _spatialKernelOrder;
        int const _spatialBgOrder;
        afwMath::Kernel::SpatialFunctionPtr _spatialKernelFunction;
        afwMath::Kernel::SpatialFunctionPtr _spatialBgFunction;
        unsigned int _nbases;     ///< Number of bases being fit for
        unsigned int _nkt;        ///< Number of kernel terms in spatial model
        unsigned int _nbt;        ///< Number of backgruond terms in spatial model
        unsigned int _nt;         ///< Total number of terms in the solution; also dimensions of matrices
        lsst::pex::policy::Policy _policy;
        bool _constantFirstTerm;  ///< Is the first term spatially invariant?
    };
}

/************************************************************************************************************/

template<typename PixelT>
std::pair<afwMath::LinearCombinationKernel::Ptr, afwMath::Kernel::SpatialFunctionPtr>
fitSpatialKernelFromCandidates(
    PsfMatchingFunctor<PixelT> &kFunctor,       ///< kFunctor used to build the kernels
    afwMath::SpatialCellSet const& kernelCells, ///< A SpatialCellSet containing KernelCandidates
    pexPolicy::Policy const& policy             ///< Policy to control the processing
                                 ) {
    typedef typename afwImage::Image<afwMath::Kernel::Pixel> ImageT;

    /* There are a variety of recipes for creating a spatial kernel which I will
     * outline here :
     *
     * 1a) Using unregularized delta function kernels, run a full spatial model
           where effectively each delta function basis varys spatially
           individually.  While this is the most general solution and may be
           fast due to the specialization of delta-function convolution, it has
           also been shown to lead to noisy kernels.  This is not recommended.

     * 1b) Using unregularized delta function kernels, do a PCA of the returned
           Kernels, and use these to create a new basis set.  This requires a
           first call to singleKernelFitter, then an instance of
           SetPcaImageVisitor() to do the PCA, creation of a new kFunctor with
           the eigenBases, a new call to singleKernelFitter using these new
           bases then a call to spatialKernelFitter.  It appaears that the
           kernels are not self-similar enough to make this a viable solution.

     * 2a) Using regularized delta function kernels, run a full spatial model
           where effectively each delta function basis varys spatially
           individually.  This merely requires repeated calls to
           singleKernelFitter and spatialKernelFitter with the supplied
           kFunctor, same as option 1a) and 3).  While this is general and may
           be fast due to the specialized delta-function convolution, we cannot
           enforce that the kernel sum does not vary spatially.

     * 2b) Using regularized delta function kernels, do a PCA of the returned
           Kernels, and use these to create a new basis set.  This requires a
           first call to singleKernelFitter, then an instance of
           SetPcaImageVisitor() to do the PCA, creation of a new kFunctor with
           the eigenBases, a new call to singleKernelFitter using these new
           bases then a call to spatialKernelFitter.  While this seems somewhat
           circuitous, we should be able to use many fewer basis functions,
           making the final image convolution faster.  We can also enforce that
           the kernel sum does not vary spatially by modifying the eigenBases.

     * 3)  Use Alard Lupton basis set.  This merely requires repeated calls to
           singleKernelFitter and spatialKernelFitter with the supplied
           kFunctor.  With these we can enforce that the kernel sum does not
           vary spatially.
     * 
     */
    

    int const maxKsumIterations       = policy.getInt("maxKsumIterations");
    int const maxSpatialIterations    = policy.getInt("maxSpatialIterations");
    int const nStarPerCell            = policy.getInt("nStarPerCell");
    int const spatialKernelOrder      = policy.getInt("spatialKernelOrder");
    int const spatialBgOrder          = policy.getInt("spatialBgOrder");
    bool const usePcaForSpatialKernel = policy.getBool("usePcaForSpatialKernel");

    boost::timer t;
    t.restart();
    
    /* The basis used for the spatial fit may change if we run Pca */
    boost::shared_ptr<afwMath::KernelList> basisListToUse;

    afwMath::LinearCombinationKernel::Ptr spatialKernel;
    afwMath::Kernel::SpatialFunctionPtr spatialBackground;

    /* Visitor for the single kernel fit */
    BuildSingleKernelVisitor<PixelT> singleKernelFitter(kFunctor, policy);

    /* Visitor for the kernel sum rejection */
    KernelSumVisitor<PixelT> kernelSumVisitor(policy);
    
    for (int i=0, nRejected=-1; i < maxSpatialIterations; i++) {
        /* Make sure there are no uninitialized candidates as current occupant of Cell */
        while (nRejected != 0) {
            pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor", 
                                  "Building single kernels A...");
            singleKernelFitter.reset();
            kernelCells.visitCandidates(&singleKernelFitter, nStarPerCell);
            nRejected = singleKernelFitter.getNRejected();
        }
        
        /* Reject outliers in kernel sum */
        for (int j=0; j < maxKsumIterations; j++) {
            kernelSumVisitor.reset();
            kernelSumVisitor.setMode(KernelSumVisitor<PixelT>::AGGREGATE);
            kernelCells.visitCandidates(&kernelSumVisitor, nStarPerCell);
            kernelSumVisitor.processKsumDistribution();
            kernelSumVisitor.setMode(KernelSumVisitor<PixelT>::REJECT);
            kernelCells.visitCandidates(&kernelSumVisitor, nStarPerCell);
            nRejected = kernelSumVisitor.getNRejected();
            pexLogging::TTrace<3>("lsst.ip.diffim.KernelSumVisitor", 
                                  "Ksum Iteration %d, rejected %d Kernels", j, nRejected);
            if (nRejected == 0) {
                break;
            }
            else {
                /* Make sure there are no uninitialized candidates as current occupant of Cell */
                while (nRejected != 0) {
                    pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor", 
                                          "Building single kernels B...");
                    singleKernelFitter.reset();
                    kernelCells.visitCandidates(&singleKernelFitter, nStarPerCell);
                    nRejected = singleKernelFitter.getNRejected();
                }
            }
        }

        /* 
           At this stage we can either apply the spatial fit to the kernels, or
           we run a PCA, use these as a *new* basis set with lower
           dimensionality, and then apply the spatial fit to these kernels 
        */
        if (usePcaForSpatialKernel) {
            int const nEigenComponents = policy.getInt("numPrincipalComponents");

            pexLogging::TTrace<5>("lsst.ip.diffim.SetPcaImageVisitor", 
                                  "Building Pca Basis");
            afwImage::ImagePca<ImageT> imagePca;
            SetPcaImageVisitor<PixelT> importStarVisitor(&imagePca);
            kernelCells.visitCandidates(&importStarVisitor, nStarPerCell);
            imagePca.analyze();
            std::vector<typename ImageT::Ptr> eigenImages = imagePca.getEigenImages();
            std::vector<double> eigenValues = imagePca.getEigenValues();
            int const nEigen = static_cast<int>(eigenValues.size());
            int const ncomp  = (nEigenComponents <= 0 || nEigen < nEigenComponents) ? nEigen : nEigenComponents;
            //
            afwMath::KernelList kernelListRaw;
            for (int i = 0; i != ncomp; ++i) {
                kernelListRaw.push_back(afwMath::Kernel::Ptr(
                                            new afwMath::FixedKernel(afwImage::Image<afwMath::Kernel::Pixel>(*eigenImages[i], true))));
            }
            /* Put all the power in the first kernel, which will not vary spatially */
            afwMath::KernelList kernelListPca = diffim::renormalizeKernelList(kernelListRaw);

            /* New PsfMatchingFunctor and Kernel visitor for this new basis list */
            diffim::PsfMatchingFunctor<PixelT> kFunctorPca(kernelListPca);
            BuildSingleKernelVisitor<PixelT> singleKernelFitterPca(kFunctorPca, policy);

            /* Always true for Pca kernel; leave original kernel alone for future Pca iterations */
            singleKernelFitterPca.setCandidateKernel(false);

            pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor", 
                                  "Rebuilding kernels using Pca basis");

            /* Only true for the first visit so we rebuild each good kernel with
             * its PcaBasis representation */
            singleKernelFitterPca.setSkipBuilt(false);
            kernelCells.visitCandidates(&singleKernelFitterPca, nStarPerCell);

            /* Once they are built we don't have to revisit */
            singleKernelFitterPca.setSkipBuilt(true);

            /* Were any rejected?  If so revisit the bad Cells */
            nRejected = singleKernelFitterPca.getNRejected();
            while (nRejected != 0) {
                pexLogging::TTrace<5>("lsst.ip.diffim.BuildSingleKernelVisitor", 
                                      "Building single kernels C...");
                singleKernelFitterPca.reset();
                kernelCells.visitCandidates(&singleKernelFitterPca, nStarPerCell);
                nRejected = singleKernelFitterPca.getNRejected();
            }
            basisListToUse.reset(new afwMath::KernelList(kFunctorPca.getBasisList()));
        }
        else {
            basisListToUse.reset(new afwMath::KernelList(kFunctor.getBasisList()));
        }
        
        /* Visitor for the spatial kernel fit */
        BuildSpatialKernelVisitor<PixelT> spatialKernelFitter(*basisListToUse, spatialKernelOrder, spatialBgOrder, policy);
        kernelCells.visitCandidates(&spatialKernelFitter, nStarPerCell);
        spatialKernelFitter.solveLinearEquation();
        std::pair<afwMath::LinearCombinationKernel::Ptr, 
            afwMath::Kernel::SpatialFunctionPtr> KB = spatialKernelFitter.getSpatialModel();
        spatialKernel     = KB.first;
        spatialBackground = KB.second;
        
        /* Visitor for the spatial kernel result */
        AssessSpatialKernelVisitor<PixelT> spatialKernelAssessor(spatialKernel, spatialBackground, policy);
        kernelCells.visitCandidates(&spatialKernelAssessor, nStarPerCell);
        nRejected = spatialKernelAssessor.getNRejected();
        pexLogging::TTrace<3>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                              "Spatial Kernel iteration %d, rejected %d Kernels", i, nRejected);
        if (nRejected == 0) {
            break;
        }
        /* 
         * Don't need to call kernelCells.visitCandidates() here since its done
         * at the begnning of the loop
        */
        
    }

    double time = t.elapsed();
    pexLogging::TTrace<3>("lsst.ip.diffim.fitSpatialKernelFromCandidates", 
                          "Total time to compute the spatial kernel : %.2f s", time);

    return std::make_pair(spatialKernel, spatialBackground);
}

/************************************************************************************************************/
//
// Explicit instantiations
//
/// \cond
    typedef float PixelT;
    template class KernelCandidate<PixelT>;

    template
    std::pair<afwMath::LinearCombinationKernel::Ptr, afwMath::Kernel::SpatialFunctionPtr>
    fitSpatialKernelFromCandidates<PixelT>(PsfMatchingFunctor<PixelT> &,
                                           lsst::afw::math::SpatialCellSet const&,
                                           lsst::pex::policy::Policy const&);
    
/// \endcond

}}} // end of namespace lsst::ip::diffim

