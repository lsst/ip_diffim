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
#include "lsst/ip/diffim/KernelSolution.h"
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
            afwMath::Kernel::SpatialFunctionPtr> kb = spatialKernelFitter.getKernelSolution();
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
        boost::shared_ptr<lsst::afw::math::KernelList> const& basisList, ///< Basis functions used in the fit
        lsst::pex::policy::Policy policy         ///< Policy file directing behavior
        ) :
        afwMath::CandidateVisitor(),
        _basisList(basisList),
        _spatialKernelOrder(policy.getInt("spatialKernelOrder")),
        _spatialBgOrder(policy.getInt("spatialBgOrder")),
        _spatialKernelFunction(new afwMath::PolynomialFunction2<double>(_spatialKernelOrder)),
        _spatialBgFunction(new afwMath::PolynomialFunction2<double>(_spatialBgOrder)),
        _nbases(_basisList->size()),
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
            _nt = (_nbases-1)*_nkt + 1 + _nbt;
        } else {
            _nt = _nbases*_nkt + _nbt;
        }

        /* This visitor creates the correct sized matrix for the solution */
        boost::shared_ptr<Eigen::MatrixXd> mMat (new Eigen::MatrixXd(_nt, _nt));
        boost::shared_ptr<Eigen::VectorXd> bVec (new Eigen::VectorXd(_nt));
        (*mMat).setZero();
        (*bVec).setZero();

        boost::shared_ptr<SpatialKernelSolution> _kernelSolution (
            new SpatialKernelSolution(mMat, bVec, 
                                      _basisList, 
                                      _spatialKernelFunction, 
                                      _spatialBgFunction,
                                      _constantFirstTerm)
            );
        
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
        if (!(kCandidate->isInitialized())) {
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
        for (int idx = 0; idx < _nkt; idx++) { paramsK[idx] = 0.0; }
        Eigen::VectorXd pK(_nkt);
        for (int idx = 0; idx < _nkt; idx++) {
            paramsK[idx] = 1.0;
            _spatialKernelFunction->setParameters(paramsK);
            pK(idx) = (*_spatialKernelFunction)(kCandidate->getXCenter(), 
                                                kCandidate->getYCenter());
            paramsK[idx] = 0.0;
        }
        Eigen::MatrixXd pKpKt = (pK * pK.transpose());
        
        /* Pure background terms */
        std::vector<double> paramsB = _spatialBgFunction->getParameters();
        for (int idx = 0; idx < _nbt; idx++) { paramsB[idx] = 0.0; }
        Eigen::VectorXd pB(_nbt);
        for (int idx = 0; idx < _nbt; idx++) {
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
        boost::shared_ptr<Eigen::MatrixXd> qMat = 
            kCandidate->getKernelSolution(KernelCandidate<PixelT>::RECENT)->getM();
        boost::shared_ptr<Eigen::VectorXd> wVec = 
            kCandidate->getKernelSolution(KernelCandidate<PixelT>::RECENT)->getB();
        
        if (DEBUG_MATRIX) {
            std::cout << "Spatial matrix inputs" << std::endl;
            std::cout << "M " << (*qMat) << std::endl;
            std::cout << "B " << (*wVec) << std::endl;
        }

        /* M and B owned by the KernelSolution */
        boost::shared_ptr<Eigen::MatrixXd> mMat = _kernelSolution->getM();
        boost::shared_ptr<Eigen::VectorXd> bVec = _kernelSolution->getB();
        
        /* first index to start the spatial blocks; default=0 for non-constant first term */
        int m0 = 0; 
        /* how many rows/cols to adjust the matrices/vectors; default=0 for non-constant first term */
        int dm = 0; 
        /* where to start the background terms; this is always true */
        int mb = _nt - _nbt;
        
        if (_constantFirstTerm) {
            m0 = 1;       /* we need to manually fill in the first (non-spatial) terms below */
            dm = _nkt-1;  /* need to shift terms due to lack of spatial variation in first term */
            
            (*mMat)(0, 0) += (*qMat)(0,0);
            for(int m2 = 1; m2 < _nbases; m2++)  {
                (*mMat).block(0, m2*_nkt-dm, 1, _nkt) += (*qMat)(0,m2) * pK.transpose();
            }
            
            (*mMat).block(0, mb, 1, _nbt) += (*qMat)(0,_nbases) * pB.transpose();
            (*bVec)(0) += (*wVec)(0);
        }
        
        /* Fill in the spatial blocks */
        for(int m1 = m0; m1 < _nbases; m1++)  {
            /* Diagonal kernel-kernel term; only use upper triangular part of pKpKt */
            (*mMat).block(m1*_nkt-dm, m1*_nkt-dm, _nkt, _nkt) += (*qMat)(m1,m1) * 
                pKpKt.part<Eigen::UpperTriangular>();
            
            /* Kernel-kernel terms */
            for(int m2 = m1+1; m2 < _nbases; m2++)  {
                (*mMat).block(m1*_nkt-dm, m2*_nkt-dm, _nkt, _nkt) += (*qMat)(m1,m2) * pKpKt;
            }
            
            /* Kernel cross terms with background */
            (*mMat).block(m1*_nkt-dm, mb, _nkt, _nbt) += (*qMat)(m1,_nbases) * pKpBt;
            
            /* B vector */
            (*bVec).segment(m1*_nkt-dm, _nkt) += (*wVec)(m1) * pK;
        }
        
        /* Background-background terms only */
        (*mMat).block(mb, mb, _nbt, _nbt) += (*qMat)(_nbases,_nbases) * 
            pBpBt.part<Eigen::UpperTriangular>();
        (*bVec).segment(mb, _nbt)         += (*wVec)(_nbases) * pB;

        /* Fill in the other half of mMat */
        for (int i = 0; i < _nt; i++) {
            for (int j = i+1; j < _nt; j++) {
                (*mMat)(j,i) = (*mMat)(i,j);
            }
        }
        
        if (DEBUG_MATRIX) {
            std::cout << "Spatial matrix outputs" << std::endl;
            std::cout << "mMat " << (*mMat) << std::endl;
            std::cout << "bVec " << (*bVec) << std::endl;
        }
        
    }

    template<typename PixelT>
    void BuildSpatialKernelVisitor<PixelT>::solveLinearEquation() {
        _kernelSolution->solve();
    }

    template<typename PixelT>
    std::pair<afwMath::LinearCombinationKernel::Ptr, afwMath::Kernel::SpatialFunctionPtr>
    BuildSpatialKernelVisitor<PixelT>::getKernelSolution() {
        return _kernelSolution->getKernelSolution();
    }

    typedef float PixelT;
    template class BuildSpatialKernelVisitor<PixelT>;
  
}}}} // end of namespace lsst::ip::diffim::detail
    
