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

#include <lsst/afw/image/Image.h>
#include <lsst/afw/image/ImagePca.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/math/FunctionLibrary.h>
#include <lsst/afw/detection/Footprint.h>

#include <lsst/pex/exceptions/Runtime.h>
#include <lsst/pex/policy/Policy.h>
#include <lsst/pex/logging/Trace.h>

#include <lsst/ip/diffim/SpatialModelKernel.h>

namespace afwMath        = lsst::afw::math;
namespace afwImage       = lsst::afw::image;

namespace lsst {
namespace ip {
namespace diffim {

template <typename ImageT, typename PixelT>
typename ImageT::ConstPtr lsst::ip::diffim::KernelCandidate<ImageT, PixelT>::getImage() const {
    int const width = getWidth() == 0 ? 15 : getWidth();
    int const height = getHeight() == 0 ? 15 : getHeight();
    
    if (_haveImage && (width != _image->getWidth() || height != _image->getHeight())) {
        _haveImage = false;
    }
    
    if (!_haveKernel) {
        lsst::pex::exceptions::NotFoundException e;
        LSST_EXCEPT_ADD(e, "No Kernel to make KernelCandidate Image from");
        throw e;
    }

    if (!_haveImage) {
        /* Calculate it from the Kernel */
        afwImage::Image<afwMath::Kernel::PixelT> _image(_kernel->getDimensions());
        (void)_kernel->computeImage(_image, false);                    
    }
    
    return _image;
}

template <typename ImageT, typename PixelT>
lsst::afw::math::Kernel::PtrT lsst::ip::diffim::KernelCandidate<ImageT, PixelT>::getKernel() const {
    if (!_haveKernel) {
        lsst::pex::exceptions::NotFoundException e;
        LSST_EXCEPT_ADD(e, "No Kernel for KernelCandidate");
        throw e;
    }
    return _kernel;
}

namespace {
    /* Lets assume this steps over the bad footprints */
    template<typename PixelT>
    class BuildKernelVisitor : public afwMath::CandidateVisitor {
        typedef afwImage::MaskedImage<PixelT> MaskedImageT;
    public:
        BuildKernelVisitor(typename PsfMatchingFunctor<PixelT>::Ptr const& kFunctor,
                           lsst::pex::policy::Policy const& policy) :
            afwMath::CandidateVisitor(),
            _kFunctor(kFunctor),
            _policy(policy),
            _imstats( ImageStatistics<PixelT>() ){}
        
        void processCandidate(afwMath::SpatialCellCandidate *candidate) {
            KernelCandidate<MaskedImageT, PixelT> *kCandidate = dynamic_cast<KernelCandidate<MaskedImageT, PixelT> *>(candidate);
            if (kCandidate == NULL) {
                throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                                  "Failed to cast SpatialCellCandidate to KernelCandidate");
            }

            /* Build its kernel here */
            MaskedImageT var = MaskedImageT(kCandidate->getMiToNotConvolvePtr(), true);
            var             -= kCandidate->getMiToConvolvePtr();

            _kFunctor.apply(kCandidate->getMiToConvolvePtr->getImage(),
                            kCandidate->getMiToNotConvolvePtr->getImage(),
                            var.getVariance(),
                            _policy);

            kCandidate->setKernel(_kFunctor.getKernel());
            kCandidate->setBackground(_kFunctor.getBackground());

            /* Make diffim and set chi2 from result */
            MaskedImageT diffim = convolveAndSubtract(kCandidate->getMiToConvolvePtr(),
                                                      kCandidate->getMiToNotConvolvePtr(),
                                                      _kFunctor.getKernel(),
                                                      _kFunctor.getBackground());
            _imstats.apply(diffim);
            kCandidate->setChi2(_imstats.getVariance());

        }
    private:
        typename PsfMatchingFunctor<PixelT>::Ptr _kFunctor;
        lsst::pex::policy::Policy _policy;
        ImageStatistics<PixelT> _imstats;
    };
}

namespace {
    template<typename PixelT>
    class SetPcaImageVisitor : public afwMath::CandidateVisitor {
        typedef afwImage::Image<PixelT> ImageT;
        typedef afwImage::MaskedImage<PixelT> MaskedImageT;
    public:
        SetPcaImageVisitor(afwImage::ImagePca<ImageT> *imagePca // Set of Images to initialise
                          ) :
            afwMath::CandidateVisitor(),
            _imagePca(imagePca) {}
        
        // Called by SpatialCellSet::visitCandidates for each Candidate
        void processCandidate(afwMath::SpatialCellCandidate *candidate) {
            KernelCandidate<MaskedImageT, PixelT> *kCandidate = dynamic_cast<KernelCandidate<MaskedImageT, PixelT> *>(candidate);
            if (kCandidate == NULL) {
                throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                                  "Failed to cast SpatialCellCandidate to KernelCandidate");
            }

            try {
                _imagePca->addImage(kCandidate->getImage()->getImage(), kCandidate->getSource().getPsfFlux());
            } catch(lsst::pex::exceptions::LengthErrorException &e) {
                return;
            }
        }
    private:
        afwImage::ImagePca<ImageT> *_imagePca; 
    };
}

namespace {
    template<typename PixelT>
    class LinearSpatialFitVisitor : public afwMath::CandidateVisitor {
        typedef afwImage::Image<PixelT> ImageT;
        typedef afwImage::MaskedImage<PixelT> MaskedImageT;
    public:
        LinearSpatialFitVisitor() :
            afwMath::CandidateVisitor() {}
        
        void processCandidate(afwMath::SpatialCellCandidate *candidate) {
            KernelCandidate<MaskedImageT, PixelT> *kCandidate = dynamic_cast<KernelCandidate<MaskedImageT, PixelT> *>(candidate);
            if (kCandidate == NULL) {
                throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                                  "Failed to cast SpatialCellCandidate to KernelCandidate");
            }
            
            
        }
    private:
        afwImage::ImagePca<ImageT> *_imagePca; 
    };
}

template<typename PixelT>
std::pair<afwMath::LinearCombinationKernel::PtrT, std::vector<double> > createPcaBasisFromCandidates(
        afwMath::SpatialCellSet const& psfCells, ///< A SpatialCellSet containing PsfCandidates
        int const nEigenComponents,              ///< number of eigen components to keep; <= 0 => infty
        int const spatialOrder,                  ///< Order of spatial variation (cf. lsst::afw::math::PolynomialFunction2)
        int const nStarPerCell                   ///< max no. of stars per cell; <= 0 => infty
    ) {
    typedef typename afwImage::Image<PixelT> ImageT;

    afwImage::ImagePca<ImageT> imagePca;
    SetPcaImageVisitor<PixelT> importStarVisitor(&imagePca);
    psfCells.visitCandidates(&importStarVisitor, nStarPerCell);
    imagePca.analyze();

    std::vector<typename ImageT::Ptr> eigenImages = imagePca.getEigenImages();
    std::vector<double> eigenValues               = imagePca.getEigenValues();
    int const nEigen = static_cast<int>(eigenValues.size());
    int const ncomp  = (nEigenComponents <= 0 || nEigen < nEigenComponents) ? nEigen : nEigenComponents;

    //
    // Now build our LinearCombinationKernel; build the lists of basis functions
    // and spatial variation, then assemble the Kernel
    //
    afwMath::KernelList<afwMath::Kernel> kernelList;
    std::vector<afwMath::Kernel::SpatialFunctionPtr> spatialFunctionList;
    
    for (int i = 0; i != ncomp; ++i) {
        kernelList.push_back(afwMath::Kernel::PtrT(
                new afwMath::FixedKernel(afwImage::Image<afwMath::Kernel::PixelT>(*eigenImages[i], true)))
            );

        afwMath::Kernel::SpatialFunctionPtr spatialFunction(new afwMath::PolynomialFunction2<double>(spatialOrder));
        if (i == 0) 
            spatialFunction->setParameter(0, 1.0); // the constant term = mean kernel; all others are 0
        spatialFunctionList.push_back(spatialFunction);
    }

    afwMath::LinearCombinationKernel::PtrT kernel(new afwMath::LinearCombinationKernel(kernelList, spatialFunctionList));
    return std::make_pair(kernel, eigenValues);
}



































#if 0
template<typename PixelT>
std::pair<afwMath::LinearCombinationKernel::PtrT, std::vector<double> > returnKernelFromCandidates(
        afwMath::SpatialCellSet const& psfCells, ///< A SpatialCellSet containing PsfCandidates
        int const spatialOrder,                  ///< Order of spatial variation (cf. lsst::afw::math::PolynomialFunction2)
        int const nStarPerCell                   ///< max no. of stars per cell; <= 0 => infty
    ) {
    typedef typename afwImage::Image<PixelT> ImageT;

    /* 
       These are linear combination kernels by definition; we want to keep the
       basis functions and not do PCA
    */
    afwMath::KernelList<afwMath::LinearCombinationKernel> kernelList;
    SetKernelVisitor<PixelT> importKernelVisitor(&kernelList);
    psfCells.visitCandidates(&importKernelVisitor, nStarPerCell);

    afwMath::KernelList<afwMath::Kernel> kernelBasisList;
    std::vector<afwMath::Kernel::SpatialFunctionPtr> spatialFunctionList;
    
    for (int i = 0; i != ncomp; ++i) {
        if (i == 0) {
            /* Grab the basis functions first time through */
            //kernelList.push_back(afwMath::Kernel::PtrT(
            //new afwMath::FixedKernel(afwImage::Image<afwMath::Kernel::PixelT>(*eigenImages[i], true)))
            //);
        }

        /* Grab the coefficients */
        //afwMath::Kernel::SpatialFunctionPtr spatialFunction(new afwMath::PolynomialFunction2<double>(spatialOrder));
    }

    afwMath::LinearCombinationKernel::PtrT kernel(new afwMath::LinearCombinationKernel(kernelList, spatialFunctionList));
    return std::make_pair(kernel, eigenValues);
}
#endif













#if 0
    
template <typename ImageT>
SpatialModelKernel<ImageT>::SpatialModelKernel(
    lsst::afw::detection::Footprint::Ptr const &fpPtr,
    MaskedImagePtr const &miToConvolvePtr,
    MaskedImagePtr const &miToNotConvolvePtr,
    boost::shared_ptr<PsfMatchingFunctor<ImageT> > const &kFunctor,
    lsst::pex::policy::Policy const &policy,
    bool build
    ) :
    _fpPtr(fpPtr),
    _miToConvolvePtr(miToConvolvePtr),
    _miToNotConvolvePtr(miToNotConvolvePtr),
    _kFunctor(kFunctor),
    _policy(policy),
    _colc(0.),
    _rowc(0.),
    _kPtr(),
    _kErrPtr(),
    _kSum(0.),
    _bg(0.),
    _bgErr(0.),
    _kStats(),
    _isBuilt(false),
    _isGood(false)
{

    if (build == true) {
        buildModel();
    }
}

template <typename ImageT>
bool SpatialModelKernel<ImageT>::buildModel() {

    if (isBuilt() == true) {
        return false;
    }

    // fill in information on position in the image
    image::BBox fpBBox = _fpPtr->getBBox();

    // NOTE : since we can't remap pixel range to go from -1 to 1 in convolve(),
    // we have to use the actual pixel value here.  Not optimal.
    // setColc(float(fpBBox.getX0() + fpBBox.getX1()) / _miToConvolveParentPtr->getWidth() - 1.0);
    // setRowc(float(fpBBox.getY0() + fpBBox.getY1()) / _miToConvolveParentPtr->getHeight() - 1.0);
    setColc(0.5 * float(fpBBox.getX0() + fpBBox.getX1()));
    setRowc(0.5 * float(fpBBox.getY0() + fpBBox.getY1()));

    logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                       "Footprint = %d,%d -> %d,%d (center %d,%d)",
                       fpBBox.getX0(), fpBBox.getY0(),
                       fpBBox.getX1(), fpBBox.getY1(),
		       int(getColc()), int(getRowc()));

    // Estimate of the variance for first kernel pass
    // True argument is for a deep copy, so -= does not modify the original pixels
    image::MaskedImage<ImageT> varEstimate = image::MaskedImage<ImageT>(*_miToNotConvolvePtr, true);
    varEstimate -= *_miToConvolvePtr;

    try {
        _kFunctor->apply(*_miToConvolvePtr->getImage(), 
                         *_miToNotConvolvePtr->getImage(),
                         *varEstimate.getVariance(), 
                         _policy);
    } catch (exceptions::Exception& e) {
        setStatus(false);
        logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                           "Exception caught from kFunctor.apply"); 
        logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                           e.what());
        return false;
    }
    math::Kernel::PtrT kernelPtr       = _kFunctor->getKernel();
    math::Kernel::PtrT kernelErrorPtr  = _kFunctor->getKernelError();
    double background =      _kFunctor->getBackground();
    double backgroundError = _kFunctor->getBackgroundError();

    // Compute kernel sum
    double kSum = 0.;
    unsigned int kCols = _policy.getInt("kernelCols");
    unsigned int kRows = _policy.getInt("kernelRows");
    image::Image<double> kImage(kCols, kRows);
    kSum = kernelPtr->computeImage(kImage, false);

    // Create difference image and calculate associated statistics
    image::MaskedImage<ImageT> diffIm = convolveAndSubtract(*_miToConvolvePtr, *_miToNotConvolvePtr,
                                                            *kernelPtr, background);

    // Find statistics of difference image 
    typename diffim::ImageStatistics<ImageT>::Ptr kStats(new diffim::ImageStatistics<ImageT>());
    (*kStats).apply(diffIm);
    logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                       "Kernel pass 1 : Kernel Sum = %.3f; Background = %.3f +/- %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                       kSum, background, backgroundError,
                       (*kStats).getMean(),
                       (*kStats).getRms());

    // A second pass with a better variance estimate from first difference image
    // An issue here is the boundary pixels are contaminated
    bool iterateKernel = _policy.getBool("iterateKernel");
    if (iterateKernel) {
        try {
            try {
                _kFunctor->apply(*_miToConvolvePtr->getImage(), 
                                 *_miToNotConvolvePtr->getImage(),
                                 *diffIm.getVariance(), 
                                 _policy);
            } catch (exceptions::Exception& e) {
                throw;
            }
            kernelPtr       = _kFunctor->getKernel();
            kernelErrorPtr  = _kFunctor->getKernelError();
            background      = _kFunctor->getBackground();
            backgroundError = _kFunctor->getBackgroundError();
            
            kSum    = 0.;
            kSum    = kernelPtr->computeImage(kImage, false);
            diffIm  = convolveAndSubtract(*_miToConvolvePtr, *_miToNotConvolvePtr, *kernelPtr, background);

            // Reset the image its looking at
            kStats->apply(diffIm);

            logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                               "Kernel pass 2 : Kernel Sum = %.3f; Background = %.3f +/- %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                               kSum, background, backgroundError, kStats->getMean(), kStats->getRms());
            
        } catch (exceptions::Exception& e) {
            logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                               "Exception caught from kFunctor.apply, reverting to first solution");
            logging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                               e.what());
        }
    }
        
    // Updates for derived class
    _kPtr    = kernelPtr;
    _kErrPtr = kernelErrorPtr;
    _kSum    = kSum;
    _bg      = background;
    _bgErr   = backgroundError;
    _kStats  = kStats;
    // Updates for base class
    setStatus(kStats->evaluateQuality(_policy));
    setBuildStatus(true);

    logging::TTrace<4>("lsst.ip.diffim.SpatialModelKernel.buildModel",
                       "Kernel : Kernel Sum = %.3f; Background = %.3f +/- %.3f; Diffim residuals = %.2f +/- %.2f sigma",
                       _kSum, background, backgroundError, _kStats->getMean(), _kStats->getRms());

    // Return quality of the kernel
    return getStatus();
}

template <typename ImageT>
double SpatialModelKernel<ImageT>::returnCellRating() {
    // Currently, just check the total flux in the template image
    FindCounts<ImageT> counter;
    counter.apply(*_miToConvolvePtr);
    return counter.getCounts();
}

#endif

}}} // end of namespace lsst::ip::diffim

