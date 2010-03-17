// -*- lsst-c++ -*-
/**
 * @file KernelSolution.cc
 *
 * @brief Implementation of KernelSolution class
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include "boost/timer.hpp" 

#include "Eigen/Core"
#include "Eigen/QR"
#include "Eigen/LU"

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/ip/diffim/ImageSubtract.h"
#include "lsst/ip/diffim/KernelSolution.h"

namespace afwMath        = lsst::afw::math;
namespace afwImage       = lsst::afw::image;
namespace pexLog         = lsst::pex::logging; 
namespace pexExcept      = lsst::pex::exceptions; 

namespace lsst { 
namespace ip { 
namespace diffim {
    
    KernelSolution::KernelSolution(
        boost::shared_ptr<Eigen::MatrixXd> mMat,
        boost::shared_ptr<Eigen::VectorXd> bVec
        ) :
        _mMat(mMat),
        _bVec(bVec),
        _sVec(),
        _solvedBy(KernelSolution::NONE)
    {};

    void KernelSolution::solve() {
        Eigen::VectorXd sVec = Eigen::VectorXd::Zero(_bVec->size());

        boost::timer t;
        t.restart();
        
        pexLog::TTrace<2>("lsst.ip.diffim.KernelSolution.solve", 
                          "Solving for kernel");
        
        _solvedBy = KernelSolution::CHOLESKY_LDLT;
        if (!(_mMat->ldlt().solve(*_bVec, &sVec))) {
            pexLog::TTrace<5>("lsst.ip.diffim.KernelSolution.solve", 
                              "Unable to determine kernel via Cholesky LDL^T");
            
            _solvedBy = KernelSolution::CHOLESKY_LLT;
            if (!(_mMat->llt().solve(*_bVec, &sVec))) {
                pexLog::TTrace<5>("lsst.ip.diffim.KernelSolution.solve", 
                                  "Unable to determine kernel via Cholesky LL^T");
                
                _solvedBy = KernelSolution::LU;
                if (!(_mMat->lu().solve(*_bVec, &sVec))) {
                    pexLog::TTrace<5>("lsst.ip.diffim.KernelSolution.solve", 
                                      "Unable to determine kernel via LU");
                    /* LAST RESORT */
                    try {
                        
                        _solvedBy = KernelSolution::EIGENVECTOR;
                        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eVecValues(_mMat);
                        Eigen::MatrixXd const& rMat = eVecValues.eigenvectors();
                        Eigen::VectorXd eValues = eVecValues.eigenvalues();
                        
                        for (int i = 0; i != eValues.rows(); ++i) {
                            if (eValues(i) != 0.0) {
                                eValues(i) = 1.0/eValues(i);
                            }
                        }
                        
                        sVec = rMat*eValues.asDiagonal()*rMat.transpose()*(*_bVec);
                    } catch (pexExcept::Exception& e) {
                        
                        _solvedBy = KernelSolution::NONE;
                        pexLog::TTrace<5>("lsst.ip.diffim.KernelSolution.solve", 
                                          "Unable to determine kernel via eigen-values");
                        
                        throw LSST_EXCEPT(pexExcept::Exception, "Unable to determine kernel solution");
                    }
                }
            }
        }

        double time = t.elapsed();
        pexLog::TTrace<3>("lsst.ip.diffim.KernelSolution.solve", 
                          "Compute time for matrix math : %.2f s", time);

        boost::shared_ptr<Eigen::VectorXd> _sVec (new Eigen::VectorXd(sVec));
    }

    /* 
     *
     */

    StaticKernelSolution::StaticKernelSolution(
        boost::shared_ptr<Eigen::MatrixXd> mMat,
        boost::shared_ptr<Eigen::VectorXd> bVec,
        boost::shared_ptr<lsst::afw::math::KernelList> const& basisList
        ) 
        :
        KernelSolution(mMat, bVec),
        _basisList(basisList),
        _kernel(),
        _background(0.0),
        _kSum(0.0),
        _kernelErr(),
        _backgroundErr(0.0),
        _errCalculated(false)
    {};

    lsst::afw::math::Kernel::Ptr StaticKernelSolution::getKernel() {
        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return solution");
        }
        return _kernel;
    }

    lsst::afw::image::Image<lsst::afw::image::VariancePixel>::Ptr StaticKernelSolution::makeKernelImage() {
        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return image");
        }
        lsst::afw::image::Image<lsst::afw::image::VariancePixel>::Ptr image (
            new lsst::afw::image::Image<lsst::afw::image::VariancePixel>::Image(_kernel->getDimensions())
            );
        _kSum  = _kernel->computeImage(*image, false);              
        return image;
    }

    double StaticKernelSolution::getBackground() {
        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return background");
        }
        return _background;
    }

    double StaticKernelSolution::getKsum() {
        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return ksum");
        }
        return _kSum;
    }

    void StaticKernelSolution::solve(bool calculateUncertainties) {
        try {
            KernelSolution::solve();
        } catch (pexExcept::Exception &e) {
            LSST_EXCEPT_ADD(e, "Unable to solve static kernel matrix");
            throw e;
        }
        
        std::pair<boost::shared_ptr<lsst::afw::math::Kernel>, double> kb = getKernelSolution();
        _kernel = kb.first;
        _background = kb.second;

        /* set kernel sum */
        (void)makeKernelImage();

        if (calculateError) {
            std::pair<boost::shared_ptr<lsst::afw::math::Kernel>, double> dkb = getKernelUncertainty();
            _kernelErr = dkb.first;
            _backgroundErr = kdb.second;
        }
    }

    std::pair<boost::shared_ptr<lsst::afw::math::Kernel>, double>
    StaticKernelSolution::getKernelSolution() {

        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return solution");
        }

        unsigned int const nParameters           = _sVec->size();
        unsigned int const nBackgroundParameters = 1;
        unsigned int const nKernelParameters     = _basisList->size();
        if (nParameters != (nKernelParameters + nBackgroundParameters)) 
            throw LSST_EXCEPT(pexExcept::Exception, "Mismatched sizes in kernel solution");

        /* Fill in the kernel results */
        std::vector<double> kValues(nKernelParameters);
        for (unsigned int idx = 0; idx < nKernelParameters; idx++) {
            if (std::isnan((*_sVec)(idx))) {
                throw LSST_EXCEPT(pexExcept::Exception, 
                                  str(boost::format("Unable to determine kernel solution %d (nan)") % idx));
            }
            kValues[idx] = (*_sVec)(idx);
        }
        boost::shared_ptr<afwMath::Kernel> kernel( 
            new afwMath::LinearCombinationKernel(_basisList, kValues) 
            );
        
        if (std::isnan((*_sVec)(nParameters-1))) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              str(boost::format("Unable to determine background solution %d (nan)") % 
                                  (nParameters-1)));
        }
        double background = (*_sVec)(nParameters-1);
        
        return std::make_pair(kernel, background);
    }        


    std::pair<boost::shared_ptr<lsst::afw::math::Kernel>, double>
    StaticKernelSolution::getSolutionUncertainty() {

        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return uncertainty");
        }

        if (_errCalculated) {
            return std::make_pair(_kernelErr, _backgroundErr);
        }

        unsigned int const nParameters           = _sVec->size();
        unsigned int const nBackgroundParameters = 1;
        unsigned int const nKernelParameters     = _basisList->size();
        if (nParameters != (nKernelParameters + nBackgroundParameters)) 
            throw LSST_EXCEPT(pexExcept::Exception, "Mismatched sizes in kernel solution");
        
        
        /* Estimate of parameter uncertainties comes from the inverse of the
         * covariance matrix (noise spectrum).  
         * N.R. 15.4.8 to 15.4.15
         * 
         * Since this is a linear problem no need to use Fisher matrix
         * N.R. 15.5.8
         *
         * Although I might be able to take advantage of the solution above.
         * Since this now works and is not the rate limiting step, keep as-is for DC3a.
         *
         * Use Cholesky decomposition again.
         * Cholkesy:
         * Cov       =  L L^t
         * Cov^(-1)  = (L L^t)^(-1)
         *           = (L^T)^-1 L^(-1)
         */
        Eigen::MatrixXd             Cov    = (*_mMat).transpose() * (*_mMat);
        Eigen::LLT<Eigen::MatrixXd> llt    = Cov.llt();
        Eigen::MatrixXd             Error2 = llt.matrixL().transpose().inverse() * llt.matrixL().inverse();
        
        std::vector<double> kErrValues(nKernelParameters);
        for (unsigned int idx = 0; idx < nKernelParameters; idx++) {
            // Insanity checking
            if (std::isnan(Error2(idx, idx))) {
                throw LSST_EXCEPT(pexExcept::Exception, 
                                  str(boost::format("Unable to determine kernel err %d (nan)") % idx));
            }
            if (Error2(idx, idx) < 0.0) {
                throw LSST_EXCEPT(pexExcept::Exception,
                                  str(boost::format("Unable to determine kernel err %d (%.3e)") % 
                                      idx % Error2(idx, idx)));
            }
            kErrValues[idx] = std::sqrt(Error2(idx, idx));
        }
        boost::shared_ptr<afwMath::Kernel> kernelErr( 
            new afwMath::LinearCombinationKernel(_basisList, kErrValues) 
            );
        
        /* Estimate of Background and Background Error */
        if (std::isnan(Error2(nParameters-1, nParameters-1))) {
            throw LSST_EXCEPT(pexExcept::Exception, "Unable to determine background err (nan)");
        }
        if (Error2(nParameters-1, nParameters-1) < 0.0) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              str(boost::format("Unable to determine background err, negative var (%.3e)") % 
                                  Error2(nParameters-1, nParameters-1) 
                                  ));
        }
        double backgroundErr = std::sqrt(Error2(nParameters-1, nParameters-1));

        _errCalculated = true;
        return std::make_pair(kernelErr, backgroundErr);
    }

    /* 
     *
     */

    SpatialKernelSolution::SpatialKernelSolution(
        boost::shared_ptr<Eigen::MatrixXd> mMat,
        boost::shared_ptr<Eigen::VectorXd> bVec,
        boost::shared_ptr<lsst::afw::math::KernelList> const& basisList,
        lsst::afw::math::Kernel::SpatialFunctionPtr spatialKernelFunction,
        lsst::afw::math::Kernel::SpatialFunctionPtr spatialBgFunction,
        bool constantFirstTerm
        ) :
        KernelSolution(mMat, bVec),
        _basisList(basisList),
        _spatialKernelFunction(spatialKernelFunction),
        _spatialBgFunction(spatialBgFunction),
        _constantFirstTerm(constantFirstTerm),
        _kernel(),
        _background(),
        _kernelErr(),
        _backgroundErr(),
        _errCalculated(false)
    {};

    void SpatialKernelSolution::solve(bool calculateUncertainties) {
        try {
            //KernelSolution::solve();
            this->solve();
        } catch (pexExcept::Exception &e) {
            LSST_EXCEPT_ADD(e, "Unable to solve spatial kernel matrix");
            throw e;
        }
            
        std::pair<afwMath::LinearCombinationKernel::Ptr, 
            afwMath::Kernel::SpatialFunctionPtr> kb = getKernelSolution();
        _kernel = kb.first;
        _background = kb.second;
        
        if (calculateUncertainties) {
            std::pair<afwMath::LinearCombinationKernel::Ptr, 
                afwMath::Kernel::SpatialFunctionPtr> dkb = getKernelUncertainty();
            _kernelErr = dkb.first;
            _backgroundErr = dkb.second;
        }
    }

    std::pair<afwMath::LinearCombinationKernel::Ptr, 
              afwMath::Kernel::SpatialFunctionPtr> SpatialKernelSolution::getKernelSolution() {

        unsigned int nkt    = _spatialKernelFunction->getParameters().size();
        unsigned int nbt    = _spatialBgFunction->getParameters().size();
        unsigned int nt     = _sVec->size();
        unsigned int nbases = _basisList->size();

        /* Set up kernel */
        afwMath::LinearCombinationKernel::Ptr spatialKernel(
            new afwMath::LinearCombinationKernel(_basisList, *_spatialKernelFunction)
            );
        
        if (nkt == 1) {
            /* Not spatially varying; this fork is a specialization
             * for convolution speed--up */
            
            /* Make sure the coefficients look right */
            if (static_cast<int>(nbases) != (nt - nbt)) {
                throw LSST_EXCEPT(pexExcept::Exception, 
                                  "Wrong number of terms for non spatially varying kernel");
            }
            
            /* Set the basis coefficients */
            std::vector<double> kCoeffs(nbases);
            for (unsigned int i = 0; i < nbases; i++) {
                kCoeffs[i] = (*_sVec)(i);
            }
            spatialKernel.reset(
                new afwMath::LinearCombinationKernel(_basisList, kCoeffs)
                );
        }
        else {
            
            /* Set the kernel coefficients */
            std::vector<std::vector<double> > kCoeffs;
            kCoeffs.reserve(nbases);
            for (unsigned int i = 0, idx = 0; i < nbases; i++) {
                kCoeffs.push_back(std::vector<double>(nkt));
                
                /* Deal with the possibility the first term doesn't vary spatially */
                if ((i == 0) && (_constantFirstTerm)) {
                    kCoeffs[i][0] = (*_sVec)(idx++);
                }
                else {
                    for (unsigned int j = 0; j < nkt; j++) {
                        kCoeffs[i][j] = (*_sVec)(idx++);
                    }
                }
            }
            spatialKernel->setSpatialParameters(kCoeffs);
        }
        
        /* Set up background */
        afwMath::Kernel::SpatialFunctionPtr bgFunction(_spatialBgFunction->clone());
        
        /* Set the background coefficients */
        std::vector<double> bgCoeffs(nbt);
        for (unsigned int i = 0; i < nbt; i++) {
            bgCoeffs[i] = (*_sVec)(nt - nbt + i);
        }
        bgFunction->setParameters(bgCoeffs);
        
        return std::make_pair(spatialKernel, bgFunction);
    }


    std::pair<afwMath::LinearCombinationKernel::Ptr, 
              afwMath::Kernel::SpatialFunctionPtr> SpatialKernelSolution::getKernelUncertainty() {
        throw LSST_EXCEPT(pexExcept::Exception, 
                          "SpatialKernelSolution::getKernelUncertainty not implemented");
    }


}}} // end of namespace lsst::ip::diffim
