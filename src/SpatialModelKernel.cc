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

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/QR>

#include <lsst/ip/diffim/SpatialModelKernel.h>

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
    int const width = getWidth() == 0 ? 15 : getWidth();
    int const height = getHeight() == 0 ? 15 : getHeight();
    
    if (_haveImage && (width != _image->getWidth() || height != _image->getHeight())) {
        _haveImage = false;
    }
    
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Kernel to make KernelCandidate Image from");
    }
    
    if (!_haveImage) {
        // Calculate it from the Kernel 
        (void)_kernel->computeImage(*_image, false);                    
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
afwMath::Kernel::PtrT KernelCandidate<PixelT>::getKernel() const {
    if (!_haveKernel) {
        throw LSST_EXCEPT(pexExcept::Exception, "No Kernel for KernelCandidate");
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
            KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
            if (kCandidate == NULL) {
                throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                                  "Failed to cast SpatialCellCandidate to KernelCandidate");
            }
            
            /* Build its kernel here */
            MaskedImageT var = MaskedImageT(kCandidate->getMiToNotConvolvePtr(), true);
            var             -= kCandidate->getMiToConvolvePtr();
            
            try {
                _kFunctor.apply(kCandidate->getMiToConvolvePtr->getImage(),
                                kCandidate->getMiToNotConvolvePtr->getImage(),
                                var.getVariance(),
                                _policy);
            } catch (lsst::pex::exceptions::Exception &e) {
                LSST_EXCEPT_ADD(e, "Unable to calculate Kernel");
                throw e;
            }
            
            /* Update the candidate with derived products */
            kCandidate->setM(_kFunctor.getM());
            kCandidate->setB(_kFunctor.getB());
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
        typedef afwImage::Image<lsst::afw::math::Kernel::PixelT> ImageT;
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
                _imagePca->addImage(kCandidate->copyImage(), kCandidate->getCandidateRating());
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
    public:
        LinearSpatialFitVisitor(
            typename PsfMatchingFunctor<PixelT>::Ptr const& kFunctor, ///< Basis functions used in the fit
            int const spatialKernelOrder,  ///< Order of spatial kernel variation (cf. lsst::afw::math::PolynomialFunction2)
            int const spatialBgOrder       ///< Order of spatial bg variation (cf. lsst::afw::math::PolynomialFunction2)
            ):
            afwMath::CandidateVisitor(),
            _M(Eigen::MatrixXd()),
            _B(Eigen::VectorXd()),
            _Soln(Eigen::VectorXd()),
            _spatialKernelOrder(spatialKernelOrder),
            _spatialBgOrder(spatialBgOrder),
            _spatialKernelFunction( new afwMath::PolynomialFunction2<double>(spatialKernelOrder) ),
            _spatialBgFunction( new afwMath::PolynomialFunction2<double>(spatialBgOrder) ),
            _nbases(kFunctor->getBasisList().size()) {
            
            /* Bookeeping terms */
            _nkt = _spatialKernelFunction->getParameters().size();
            _nbt = _spatialBgFunction->getParameters().size();

            /* Nbases + 1 term for background */
            _M.resize((_nbases+1)*(_nkt + _nbt), (_nbases+1)*(_nkt + _nbt));
            _B.resize((_nbases+1)*(_nkt + _nbt));
            _M.setZero();
            _B.setZero();
        }
        
        void processCandidate(afwMath::SpatialCellCandidate *candidate) {
            KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
            if (kCandidate == NULL) {
                throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                                  "Failed to cast SpatialCellCandidate to KernelCandidate");
            }
            
            /* Calculate P matrices */
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
            
            
            /* Add 'em to the M, B matrix */
            Eigen::MatrixXd Q   = kCandidate->getM();
            Eigen::VectorXd W   = kCandidate->getB();
            
            for(unsigned int m1 = 0; m1 < _nkt; m1++)  {
                for(unsigned int m2 = m1; m2 < _nkt; m2++)  {
                    _M.block(m1*_nkt, m2*_nkt, _nkt, _nkt) += Q(m1,m2) * PkPkt;
                    _M.block(m2*_nkt, m1*_nkt, _nkt, _nkt) += Q(m2,m1) * PkPkt;
                }
                _B.segment(m1*_nkt,_nkt)                   += W(m1) * Pk;
            } 
            
            /* shift things by m0 so as to not mix the background and kernel terms */
            for(unsigned int m1 = 0, m0 = _nkt*_nkt; m1 < _nbt; m1++)  {
                for(unsigned int m2 = m1; m2 < _nbt; m2++)  {
                    _M.block(m0+m1*_nbt, m0+m2*_nbt, _nbt, _nbt) += Q(m1+_nkt,m2+_nkt) * PbPbt;
                    _M.block(m0+m2*_nbt, m0+m1*_nbt, _nbt, _nbt) += Q(m2+_nkt,m1+_nkt) * PbPbt;
                }
                _B.segment(m0+m1*_nbt, _nbt)                     += W(m1+_nkt) * Pb;
            } 
            
            
        }
        
        void solveLinearEquation() {
            _Soln = Eigen::VectorXd::Zero(_nkt + _nbt);
            
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
        }

        
        Eigen::VectorXd getSolution() {
            return _Soln; 
        }

    private:
        Eigen::MatrixXd _M;    ///< Least squares matrix
        Eigen::VectorXd _B;    ///< Least squares vector
        Eigen::VectorXd _Soln; ///< Least squares solution
        int const _spatialKernelOrder;
        int const _spatialBgOrder;
        afwMath::Kernel::SpatialFunctionPtr _spatialKernelFunction;
        afwMath::Kernel::SpatialFunctionPtr _spatialBgFunction;
        unsigned int _nbases;  ///< Number of bases being fit for
        unsigned int _nkt;     ///< Number of kernel terms in spatial model
        unsigned int _nbt;     ///< Number of backgruond terms in spatial model
    };
}

/************************************************************************************************************/

template<typename PixelT>
std::pair<afwMath::LinearCombinationKernel::PtrT, afwMath::Kernel::SpatialFunctionPtr>
fitSpatialKernelFromCandidates(
    typename PsfMatchingFunctor<PixelT>::Ptr const& kFunctor, ///< kFunctor used to build the kernels
    afwMath::SpatialCellSet const& psfCells,                  ///< A SpatialCellSet containing PsfCandidates
    pexPolicy::Policy const& policy                           ///< Policy to control the processing
                                 ) {
    
    int const nStarPerCell       = policy.getInt("nStarPerCell");
    int const spatialKernelOrder = policy.getInt("spatialKernelOrder");
    int const spatialBgOrder     = policy.getInt("spatialBgOrder");
    
    /* Do the linear fit */
    LinearSpatialFitVisitor<PixelT> linearFitter(kFunctor, spatialKernelOrder, spatialBgOrder);
    psfCells.visitCandidates(&linearFitter, nStarPerCell);
    linearFitter.solveLinearEquation();
    Eigen::VectorXd solution = linearFitter.getSolution();

    /* Set up kernel */
    std::vector<afwMath::Kernel::SpatialFunctionPtr> spatialFunctionList;
    afwMath::KernelList<afwMath::Kernel> kernelList = kFunctor->getBasisList();
    unsigned int const nComponents                  = kernelList.size();
    for (unsigned int i = 0; i < nComponents; i++) {
        afwMath::Kernel::SpatialFunctionPtr spatialFunction(new afwMath::PolynomialFunction2<double>(spatialKernelOrder));
        spatialFunctionList.push_back(spatialFunction);
    }
    afwMath::LinearCombinationKernel::PtrT spatialKernel(new afwMath::LinearCombinationKernel(kernelList, spatialFunctionList));

    /* Set up background */
    afwMath::Kernel::SpatialFunctionPtr bgFunction(new afwMath::PolynomialFunction2<double>(spatialBgOrder));
    
    /* Set the coefficients */
    unsigned int nk = spatialFunctionList[0]->getParameters().size();
    unsigned int nb = bgFunction->getParameters().size();
    
    std::vector<std::vector<double> > kCoeffs;
    kCoeffs.reserve(nComponents);
    for (unsigned int i = 0, idx = 0; i < nComponents; i++, idx++) {
        kCoeffs.push_back(std::vector<double>(nk));
        for (unsigned int j = 0; j < nk; j++) {
            kCoeffs[i][j] = solution[idx];
        }
    }

    std::vector<double> bgCoeffs;
    for (unsigned int i = 0; i < nb; i++) {
        bgCoeffs[i] = solution[i + nComponents*nk];
    }
    
    spatialKernel->setSpatialParameters(kCoeffs);
    bgFunction->setParameters(bgCoeffs);
    
    return std::make_pair(spatialKernel, bgFunction);
}

/************************************************************************************************************/

template<typename PixelT>
std::pair<afwMath::LinearCombinationKernel::PtrT, std::vector<double> > createPcaBasisFromCandidates(
    afwMath::SpatialCellSet const& psfCells, ///< A SpatialCellSet containing PsfCandidates
    pexPolicy::Policy const& policy  ///< Policy to control the processing
    ) {
    typedef typename afwImage::Image<lsst::afw::math::Kernel::PixelT> ImageT;

    int const nEigenComponents   = policy.getInt("nEigenComponents");   // number of eigen components to keep; <= 0 => infty
    int const nStarPerCell       = policy.getInt("nStarPerCell");       // order of spatial variation
    int const spatialKernelOrder = policy.getInt("spatialKernelOrder"); // max no. of stars per cell; <= 0 => infty
    
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
        
        afwMath::Kernel::SpatialFunctionPtr spatialFunction(new afwMath::PolynomialFunction2<double>(spatialKernelOrder));
        if (i == 0) 
            spatialFunction->setParameter(0, 1.0); // the constant term = mean kernel; all others are 0
        spatialFunctionList.push_back(spatialFunction);
    }
    
    afwMath::LinearCombinationKernel::PtrT kernel(new afwMath::LinearCombinationKernel(kernelList, spatialFunctionList));
    return std::make_pair(kernel, eigenValues);
}

/************************************************************************************************************/
//
// Explicit instantiations
//
/// \cond
    typedef float PixelT;
    template class KernelCandidate<PixelT>;

    template
    std::pair<lsst::afw::math::LinearCombinationKernel::PtrT, std::vector<double> >
    createPcaBasisFromCandidates<PixelT>(lsst::afw::math::SpatialCellSet const&,
                                         lsst::pex::policy::Policy const&);

    template
    std::pair<afwMath::LinearCombinationKernel::PtrT, afwMath::Kernel::SpatialFunctionPtr>
    fitSpatialKernelFromCandidates<PixelT>(PsfMatchingFunctor<PixelT>::Ptr const&,
                                           lsst::afw::math::SpatialCellSet const&,
                                           lsst::pex::policy::Policy const&);
    
/// \endcond

}}} // end of namespace lsst::ip::diffim

