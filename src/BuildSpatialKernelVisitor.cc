// -*- lsst-c++ -*-
/**
 * @file BuildSpatialKernelVisitor.h
 *
 * @brief Implementation of BuildSpatialKernelVisitor 
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#include "boost/shared_ptr.hpp" 
#include "boost/timer.hpp" 

#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/LU"
#include "Eigen/QR"

#include "lsst/afw/math.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/ip/diffim/KernelCandidate.h"
#include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"

#define DEBUG_MATRIX 0

namespace afwMath        = lsst::afw::math;
namespace pexLogging     = lsst::pex::logging; 
namespace pexPolicy      = lsst::pex::policy; 
namespace pexExcept      = lsst::pex::exceptions; 

namespace lsst { 
namespace ip { 
namespace diffim {
namespace detail {
    /**
     * @class BuildSpatialKernelVisitor
     * @ingroup ip_diffim
     *
     * @brief Creates a spatial kernel and background from a list of candidates
     *
     * @code
        Policy::Ptr policy(new Policy);
        policy->set("spatialKernelOrder", spatialKernelOrder);
        policy->set("spatialBgOrder", spatialBgOrder);
        policy->set("kernelBasisSet", "delta-function");
        policy->set("usePcaForSpatialKernel", true);
    
        detail::BuildSpatialKernelVisitor<PixelT> spatialKernelFitter(*basisListToUse, 
                                                                      *policy);
        kernelCells.visitCandidates(&spatialKernelFitter, nStarPerCell);
        spatialKernelFitter.solveLinearEquation();
        std::pair<afwMath::LinearCombinationKernel::Ptr, 
            afwMath::Kernel::SpatialFunctionPtr> kb = spatialKernelFitter.getSpatialModel();
        spatialKernel     = kb.first;
        spatialBackground = kb.second; 
     * @endcode
     *
     * @note After visiting all candidates, solveLinearEquation() must be called to
     * trigger the matrix math.
     *
     * @note The user has the option to enfore conservation of the kernel sum across
     * the image through the policy.  In this case, all terms but the first are fit
     * for spatial variation.  This requires a little extra code to make sure the
     * matrices are the correct size, and that it is accessing the appropriate terms
     * in the matrices when creating the spatial models.
     * 
     */
    template<typename PixelT>
    BuildSpatialKernelVisitor<PixelT>::BuildSpatialKernelVisitor(
        lsst::afw::math::KernelList basisList,   ///< Basis functions used in the fit
        lsst::pex::policy::Policy policy         ///< Policy file directing behavior
        ) :
        afwMath::CandidateVisitor(),
        _basisList(basisList),
        _mMat(Eigen::MatrixXd()),
        _bVec(Eigen::VectorXd()),
        _sVec(Eigen::VectorXd()),
        _spatialKernelOrder(policy.getInt("spatialKernelOrder")),
        _spatialBgOrder(policy.getInt("spatialBgOrder")),
        _spatialKernelFunction(new afwMath::PolynomialFunction2<double>(_spatialKernelOrder)),
        _spatialBgFunction(new afwMath::PolynomialFunction2<double>(_spatialBgOrder)),
        _nbases(basisList.size()),
        _policy(policy),
        _constantFirstTerm(false),
        _nCandidates(0) {
        
        /* 
           NOTE : The variable _constantFirstTerm allows that the first
           component of the basisList has no spatial variation.  This is useful
           to conserve the kernel sum across the image.  There are 2 ways to
           implement this : we can make the input matrices/vectors smaller by
           (_nkt-1), or we can create _M and _B with empty values for the first
           component's spatial terms.  The latter could cause problems for the
           matrix math even though its probably more readable, so we go with the
           former.
        */
        bool isAlardLupton = _policy.getString("kernelBasisSet") == "alard-lupton";
        bool usePca        = _policy.getBool("usePcaForSpatialKernel");
        if (isAlardLupton || usePca) {
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
        _mMat.resize(_nt, _nt);
        _bVec.resize(_nt);
        _mMat.setZero();
        _bVec.setZero();
        
        pexLogging::TTrace<5>("lsst.ip.diffim.LinearSpatialFitVisitor", 
                              "Initializing with size %d %d %d and constant first term = %s",
                              _nkt, _nbt, _nt,
                              _constantFirstTerm ? "true" : "false");
    }
    
    
    template<typename PixelT>
    void BuildSpatialKernelVisitor<PixelT>::processCandidate(
        lsst::afw::math::SpatialCellCandidate *candidate
        ) {
        KernelCandidate<PixelT> *kCandidate = dynamic_cast<KernelCandidate<PixelT> *>(candidate);
        if (kCandidate == NULL) {
            throw LSST_EXCEPT(pexExcept::LogicErrorException,
                              "Failed to cast SpatialCellCandidate to KernelCandidate");
        }
        if (!(kCandidate->hasKernel())) {
            kCandidate->setStatus(afwMath::SpatialCellCandidate::BAD);
            pexLogging::TTrace<3>("lsst.ip.diffim.BuildSpatialKernelVisitor.processCandidate", 
                                  "Cannot process candidate %d, continuing", kCandidate->getId());
            return;
        }
        
        pexLogging::TTrace<6>("lsst.ip.diffim.BuildSpatialKernelVisitor.processCandidate", 
                              "Processing candidate %d", kCandidate->getId());
        _nCandidates += 1;

        /* Calculate P matrices */
        /* Pure kernel terms */
        std::vector<double> paramsK = _spatialKernelFunction->getParameters();
        for (unsigned int idx = 0; idx < _nkt; idx++) { paramsK[idx] = 0.0; }
        Eigen::VectorXd pK(_nkt);
        for (unsigned int idx = 0; idx < _nkt; idx++) {
            paramsK[idx] = 1.0;
            _spatialKernelFunction->setParameters(paramsK);
            pK(idx) = (*_spatialKernelFunction)(kCandidate->getXCenter(), 
                                                kCandidate->getYCenter());
            paramsK[idx] = 0.0;
        }
        Eigen::MatrixXd pKpKt = (pK * pK.transpose());
        
        /* Pure background terms */
        std::vector<double> paramsB = _spatialBgFunction->getParameters();
        for (unsigned int idx = 0; idx < _nbt; idx++) { paramsB[idx] = 0.0; }
        Eigen::VectorXd pB(_nbt);
        for (unsigned int idx = 0; idx < _nbt; idx++) {
            paramsB[idx] = 1.0;
            _spatialBgFunction->setParameters(paramsB);
            pB(idx) = (*_spatialBgFunction)(kCandidate->getXCenter(), 
                                            kCandidate->getYCenter());
            paramsB[idx] = 0.0;
        }
        Eigen::MatrixXd pBpBt = (pB * pB.transpose());
        
        /* Cross terms */
        Eigen::MatrixXd pKpBt = (pK * pB.transpose());
        
        if (DEBUG_MATRIX) {
            std::cout << "Spatial weights" << std::endl;
            std::cout << "pKpKt " << pKpKt << std::endl;
            std::cout << "pBpBt " << pBpBt << std::endl;
            std::cout << "pKpBt " << pKpBt << std::endl;
        }
        
        /* Add each candidate to the M, B matrix */
        boost::shared_ptr<Eigen::MatrixXd> qMat = kCandidate->getM();
        boost::shared_ptr<Eigen::VectorXd> wVec = kCandidate->getB();
        
        if (DEBUG_MATRIX) {
            std::cout << "Spatial matrix inputs" << std::endl;
            std::cout << "M " << (*qMat) << std::endl;
            std::cout << "B " << (*wVec) << std::endl;
        }
        
        /* first index to start the spatial blocks; default=0 for non-constant first term */
        unsigned int m0 = 0; 
        /* how many rows/cols to adjust the matrices/vectors; default=0 for non-constant first term */
        unsigned int dm = 0; 
        /* where to start the background terms; this is always true */
        unsigned int mb = _nt - _nbt;
        
        if (_constantFirstTerm) {
            m0 = 1;       /* we need to manually fill in the first (non-spatial) terms below */
            dm = _nkt-1;  /* need to shift terms due to lack of spatial variation in first term */
            
            _mMat(0, 0) += (*qMat)(0,0);
            for(unsigned int m2 = 1; m2 < _nbases; m2++)  {
                _mMat.block(0, m2*_nkt-dm, 1, _nkt) += (*qMat)(0,m2) * pK.transpose();
            }
            
            _mMat.block(0, mb, 1, _nbt) += (*qMat)(0,_nbases) * pB.transpose();
            _bVec(0) += (*wVec)(0);
        }
        
        /* Fill in the spatial blocks */
        for(unsigned int m1 = m0; m1 < _nbases; m1++)  {
            /* Diagonal kernel-kernel term; only use upper triangular part of pKpKt */
            _mMat.block(m1*_nkt-dm, m1*_nkt-dm, _nkt, _nkt) += (*qMat)(m1,m1) * 
                pKpKt.part<Eigen::UpperTriangular>();
            
            /* Kernel-kernel terms */
            for(unsigned int m2 = m1+1; m2 < _nbases; m2++)  {
                _mMat.block(m1*_nkt-dm, m2*_nkt-dm, _nkt, _nkt) += (*qMat)(m1,m2) * pKpKt;
            }
            
            /* Kernel cross terms with background */
            _mMat.block(m1*_nkt-dm, mb, _nkt, _nbt) += (*qMat)(m1,_nbases) * pKpBt;
            
            /* B vector */
            _bVec.segment(m1*_nkt-dm, _nkt) += (*wVec)(m1) * pK;
        }
        
        /* Background-background terms only */
        _mMat.block(mb, mb, _nbt, _nbt) += (*qMat)(_nbases,_nbases) * 
            pBpBt.part<Eigen::UpperTriangular>();
        _bVec.segment(mb, _nbt)         += (*wVec)(_nbases) * pB;
        
        if (DEBUG_MATRIX) {
            std::cout << "Spatial matrix outputs" << std::endl;
            std::cout << "_mMat " << _mMat << std::endl;
            std::cout << "_bVec " << _bVec << std::endl;
        }
        
    }

    template<typename PixelT>
    void BuildSpatialKernelVisitor<PixelT>::solveLinearEquation() {
        boost::timer t;
        t.restart();

        pexLogging::TTrace<2>("lsst.ip.diffim.SpatialModelKernel.solveLinearEquation", 
                              "Solving for spatial model");
        
        /* Fill in the other half of _mMat */
        for (unsigned int i = 0; i < _nt; i++) {
            for (unsigned int j = i+1; j < _nt; j++) {
                _mMat(j,i) = _mMat(i,j);
            }
        }
        _sVec = Eigen::VectorXd::Zero(_nt);
        
        if (DEBUG_MATRIX) {
            std::cout << "Solving for _mMat:" << std::endl;
            std::cout << _mMat << std::endl;
            std::cout << _bVec << std::endl;
        }
        
        if (!(_mMat.ldlt().solve(_bVec, &_sVec))) {
            pexLogging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.solveLinearEquation", 
                                  "Unable to determine kernel via Cholesky LDL^T");
            if (!(_mMat.llt().solve(_bVec, &_sVec))) {
                pexLogging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.solveLinearEquation", 
                                      "Unable to determine kernel via Cholesky LL^T");
                if (!(_mMat.lu().solve(_bVec, &_sVec))) {
                    pexLogging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.solveLinearEquation", 
                                          "Unable to determine kernel via LU");
                    // LAST RESORT
                    try {
                        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eVecValues(_mMat);
                        Eigen::MatrixXd const& R = eVecValues.eigenvectors();
                        Eigen::VectorXd eValues  = eVecValues.eigenvalues();
                        
                        for (int i = 0; i != eValues.rows(); ++i) {
                            if (eValues(i) != 0.0) {
                                eValues(i) = 1.0/eValues(i);
                            }
                        }
                        
                        _sVec = R*eValues.asDiagonal()*R.transpose()*_bVec;
                    } catch (pexExcept::Exception& e) {
                        pexLogging::TTrace<5>("lsst.ip.diffim.SpatialModelKernel.solveLinearEquation", 
                                              "Unable to determine kernel via eigen-values");
                        
                        throw LSST_EXCEPT(pexExcept::Exception, 
                                          "Unable to determine kernel solution");
                    }
                }
            }
        }
        
        if (DEBUG_MATRIX) {
            std::cout << "Solution:" << std::endl;
            std::cout << _sVec << std::endl;
        }
        
        double time = t.elapsed();
        pexLogging::TTrace<3>("lsst.ip.diffim.SpatialModelKernel.solveLinearEquation", 
                              "Compute time to do spatial matrix math : %.2f s", time);
    }

    template<typename PixelT>
    std::pair<afwMath::LinearCombinationKernel::Ptr, afwMath::Kernel::SpatialFunctionPtr>
    BuildSpatialKernelVisitor<PixelT>::getSpatialModel() {
        /* Set up kernel */
        afwMath::LinearCombinationKernel::Ptr spatialKernel(
            new afwMath::LinearCombinationKernel(_basisList, *_spatialKernelFunction)
            );
        
        if (_spatialKernelOrder == 0) {
            /* Not spatially varying; specialization for convolution speed--up */
            
            /* Make sure the coefficients look right */
            if (static_cast<int>(_nbases) != (_sVec.size()-_nbt)) {
                throw LSST_EXCEPT(pexExcept::Exception, 
                                  "Wrong number of terms for non spatially varying kernel");
            }
            
            /* Set the basis coefficients */
            std::vector<double> kCoeffs(_nbases);
            for (unsigned int i = 0; i < _nbases; i++) {
                kCoeffs[i] = _sVec[i];
            }
            spatialKernel.reset(
                new afwMath::LinearCombinationKernel(_basisList, kCoeffs)
                );
            
        }
        else {
            /* Spatially varying */
            
            /* Set the kernel coefficients */
            std::vector<std::vector<double> > kCoeffs;
            kCoeffs.reserve(_nbases);
            for (unsigned int i = 0, idx = 0; i < _nbases; i++) {
                kCoeffs.push_back(std::vector<double>(_nkt));
                
                /* Deal with the possibility the first term doesn't vary spatially */
                if ((i == 0) && (_constantFirstTerm)) {
                    kCoeffs[i][0] = _sVec[idx++];
                }
                else {
                    for (unsigned int j = 0; j < _nkt; j++) {
                        kCoeffs[i][j] = _sVec[idx++];
                    }
                }
            }
            spatialKernel->setSpatialParameters(kCoeffs);
        }
        
        /* Set up background */
        afwMath::Kernel::SpatialFunctionPtr bgFunction(_spatialBgFunction->clone());
        
        /* Set the background coefficients */
        std::vector<double> bgCoeffs(_nbt);
        for (unsigned int i = 0; i < _nbt; i++) {
            bgCoeffs[i] = _sVec[_nt - _nbt + i];
        }
        bgFunction->setParameters(bgCoeffs);
        
        return std::make_pair(spatialKernel, bgFunction);
    }

    typedef float PixelT;
    template class BuildSpatialKernelVisitor<PixelT>;
  
}}}} // end of namespace lsst::ip::diffim::detail
    
