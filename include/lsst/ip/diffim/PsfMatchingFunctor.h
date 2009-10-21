// -*- lsst-c++ -*-
/**
 * @file PsfMatchingFunctor.h
 *
 * @brief Class owning and implementing the core functionality of building a single Psf matching kernel
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_PSFMATCHINGFUNCTOR_H
#define LSST_IP_DIFFIM_PSFMATCHINGFUNCTOR_H

#include "Eigen/Core"

#include "boost/shared_ptr.hpp"

#include "lsst/pex/policy/Policy.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/image/MaskedImage.h"

namespace lsst { 
namespace ip { 
namespace diffim {
   
   /**
     * @brief Functor to create PSF Matching Kernel
     *
     * @note This class owns the functionality to make a single difference
     * imaging kernel around one object realized in 2 different images.  If
     * constructed with a regularization matrix, will use it by default.  This
     * creates the M and B vectors that are used to solve for the kernel
     * parameters 'x' as in Mx = B.  This creates a single kernel around a
     * single object, and operates in tandem with the KernelCandidate +
     * BuildSingleKernelVisitor classes for the spatial modeling.
     * 
     * @ingroup ip_diffim
     */
    template <typename PixelT, typename VarT=lsst::afw::image::VariancePixel>
    class PsfMatchingFunctor {
    public:
        typedef boost::shared_ptr<PsfMatchingFunctor> Ptr;
        typedef typename lsst::afw::image::MaskedImage<PixelT>::xy_locator xy_locator;
        typedef typename lsst::afw::image::Image<VarT>::xy_locator         xyi_locator;

        PsfMatchingFunctor(
            lsst::afw::math::KernelList const& basisList
            );
        PsfMatchingFunctor(
            lsst::afw::math::KernelList const& basisList,
            boost::shared_ptr<Eigen::MatrixXd> const& _hMat
            );
        virtual ~PsfMatchingFunctor() {};

        /* Shallow copy only; shared matrix product uninitialized */
        PsfMatchingFunctor(const PsfMatchingFunctor<PixelT,VarT> &rhs);

        std::pair<boost::shared_ptr<lsst::afw::math::Kernel>, double> getSolution();
        std::pair<boost::shared_ptr<lsst::afw::math::Kernel>, double> getSolutionUncertainty();
        
        /** Access to least squares info
         */
        std::pair<boost::shared_ptr<Eigen::MatrixXd>, boost::shared_ptr<Eigen::VectorXd> > getAndClearMB();

        /** Access to basis list
         */
        lsst::afw::math::KernelList getBasisList() const { return _basisList; }

        /* Create PSF matching kernel */
        void apply(lsst::afw::image::Image<PixelT> const& imageToConvolve,
                   lsst::afw::image::Image<PixelT> const& imageToNotConvolve,
                   lsst::afw::image::Image<VarT>   const& varEstimate,
                   lsst::pex::policy::Policy       const& policy
            );

    protected:
        lsst::afw::math::KernelList const _basisList;            ///< List of Kernel basis functions
        boost::shared_ptr<Eigen::MatrixXd> _mMat;                ///< Least squares matrix
        boost::shared_ptr<Eigen::VectorXd> _bVec;                ///< Least squares vector
        boost::shared_ptr<Eigen::VectorXd> _sVec;                ///< Least square solution
        boost::shared_ptr<Eigen::MatrixXd> const _hMat;          ///< Regularization matrix
        bool _initialized;                                       ///< Has been solved for
        bool _regularize;                                        ///< Has a _hMat matrix
    };
    
    /**
     * @brief Helper method to return a pointer to a PsfMatchingFunctor()
     *
     * @param basisList  Input set of basis kernels to use for Psf matching
     *
     * @ingroup ip_diffim
     */
    template <typename PixelT>
    typename PsfMatchingFunctor<PixelT>::Ptr
    makePsfMatchingFunctor(lsst::afw::math::KernelList const& basisList) {
        return typename PsfMatchingFunctor<PixelT>::Ptr(new PsfMatchingFunctor<PixelT>(basisList));
    }

    /**
     * @brief Helper method to return a pointer to a PsfMatchingFunctor() with regularization
     *
     * @param basisList  Input set of basis kernels to use for Psf matching
     * @param _hMat  Regularization matrix (for delta-function bases)
     *
     * @ingroup ip_diffim
     */
    template <typename PixelT>
    typename PsfMatchingFunctor<PixelT>::Ptr
    makePsfMatchingFunctor(lsst::afw::math::KernelList const& basisList,
                           boost::shared_ptr<Eigen::MatrixXd> const _hMat) {
        return typename PsfMatchingFunctor<PixelT>::Ptr(new PsfMatchingFunctor<PixelT>(basisList, _hMat));
    }

}}} // end of namespace lsst::ip::diffim

#endif
